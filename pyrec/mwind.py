"""Python version of mwind.f (YREC development, 09c3b9b)."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.amcalc import amcalc
from pyrec.wind import wind

# Fortran SAVEs its locals. TCZ is set only when LROSSBY is true, but read
# whenever LBSCALE is true, so its value can carry over from an earlier
# call; Fortran leaves its first value undefined (0.0 here).
_tcz = 0.0


def mwind(bl, delts, hsbot, hstop, jstart, jend, ljdot, omegas,
          smass, teffl, hicz, hjm):
    """Remove angular momentum from the surface convection zone by a wind.

    Uses the magnetized-wind law of Matt et al. (2012) with their
    centrifugal correction and optional Rossby scaling, over one time step
    (split into substeps if needed); hands over to WIND (the older law)
    if LMWIND is false. The loss of specific angular momentum is
    subtracted from every shell of the zone.

    Translates: mwind.f#MWIND (09c3b9b)

    Parameters
    ----------
    bl : float
        log10 L/Lsun.
    delts : float
        Time step (s).
    hsbot, hstop : float
        Masses at the bottom and top of the convection zone.
    jstart, jend : int
        First and last shells of the zone.
    ljdot : bool
        Passed on to WIND only.
    omegas : float
        Surface angular velocity at the start of the step.
    smass : float
        Mass in solar units.
    teffl : float
        log10 Teff.
    hicz : float
        Moment of inertia of the zone.
    hjm : ndarray, shape (JSON + 1,)
        Specific angular momentum, changed in place for shells
        ``jstart..jend``.

    Notes
    -----
    Sets ``nmlpar_mod.wnew`` (the new surface angular velocity), as the
    Fortran does through the module.
    """
    global _tcz
    c = constants_mod
    nml = nmlpar_mod

    # If desired, remove angular momentum using a Weber-Davis wind model;
    # otherwise use the old law.
    if not nml.lmwind:
        wind(bl, delts, hsbot, hstop, jstart, jend, ljdot, omegas,
             smass, teffl, hicz, hjm)
        return  # label 9999

    # Rossby scaling, if desired; tau_cz interpolated within the step.
    if nml.lrossby:
        _tcz = (ovrtrn_mod.taucz0
                + ovrtrn_mod.fracstep * (ovrtrn_mod.taucz - ovrtrn_mod.taucz0))
        if nml.lbscale:
            wnow = omegas * _tcz / nml.pmmsoltau
            wsat = nml.wmax
        else:
            wnow = omegas
            wsat = nml.wmax * nml.pmmsoltau / _tcz
    else:
        wnow = omegas
        wsat = nml.wmax

    # The structure factor of the Matt et al. (2012) law.
    structfactor = amcalc(smass, bl, teffl, 0.0)

    # Centrifugal reduction term (Matt et al. 2012, ApJ 754, L26),
    # relative to the Sun. C_2 = 0.0506.
    fsun = (0.5 * nml.pmmsolw ** 2 * c.crsun ** 3
            / math.exp(c.cln * c.cgl) / c.cmsun)
    # Radius.
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    fcorr = (0.5 * omegas ** 2 * math.exp(c.cln * (3.0 * rl - c.cgl))
             / smass / c.cmsun)
    fcen = ((c.c_2 ** 2 + fsun) / (c.c_2 ** 2 + fcorr)) ** nml.excen
    # First guess of the change in Omega; if it exceeds Omega, use
    # several smaller steps so that Omega stays positive.
    dwtest = ((delts / hicz) * nml.constfactor * structfactor * omegas
              * min(wnow, wsat) ** (nml.exw - 1.0) * fcen)
    if dwtest > omegas:
        nstep = int(dwtest / omegas) + 1
        dt = delts / float(nstep)
    else:
        nstep = 1
        dt = delts

    ws = omegas
    for n in range(1, nstep + 1):  # label 100
        # Question (see the page): the comment here on PARMIN and ALFA
        # describes the old law.
        # Iterative solution: first guess uses the starting Omega; then
        # the average of the starting and the previous ending Omega.
        k = 0
        w = ws
        wnew0 = ws
        while True:  # label 5
            # Question (see the page): _tcz is undefined here when
            # LBSCALE is true but LROSSBY is false.
            if nml.lbscale:
                wnow = w * _tcz / nml.pmmsoltau
            else:
                wnow = w
            k = k + 1
            fcorr = (0.5 * w ** 2 * math.exp(c.cln * (3.0 * rl - c.cgl))
                     / smass / c.cmsun)
            fcen = ((c.c_2 ** 2 + fsun) / (c.c_2 ** 2 + fcorr)) ** nml.excen
            # Question (see the page): WNEW is the namelist input in
            # nmlpar_mod, overwritten here.
            nml.wnew = (ws - (dt / hicz) * nml.constfactor * structfactor
                        * w * min(wnow, wsat) ** (nml.exw - 1.0) * fcen)
            dw = 2.0 * abs((wnew0 - nml.wnew) / (wnew0 + nml.wnew))
            if dw > 1.0e-6:
                w = 0.5 * (ws + nml.wnew)
                wnew0 = nml.wnew
                # Question (see the page): after 21 passes the loop ends
                # without a warning.
                if k <= 20:
                    continue
            break
        ws = nml.wnew

    # DM is the total mass in the convection zone; subtract the change in
    # angular momentum per unit mass from every shell of the zone.
    dm = hstop - hsbot
    fjdom = (omegas - nml.wnew) * hicz / dm
    for ii in range(jstart, jend + 1):  # label 10
        hjm[ii] = hjm[ii] - fjdom
