"""Python version of wind.f (YREC 09c3b9b): the older magnetized-wind loss
of angular momentum from the surface convection zone."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod


def wind(bl, delts, hsbot, hstop, jstart, jend, ljdot, omegas,
         smass, teffl, hicz, hjm):
    """Remove angular momentum from the surface convection zone over one
    step by the Weber-Davis wind law
    I dOmega/dt = -K (Mdot/1e-14)^a (R/Rsun)^c M^d Omega min(Omega, Omega_sat)^(b-1),
    integrated in substeps by iteration, and subtract the loss of specific
    angular momentum from every shell of the zone.

    Translates: wind.f#WIND (09c3b9b)

    Arguments: log L bl; the step delts (s); the zone's mass range hsbot,
    hstop and shells jstart..jend; ljdot, the switch; the surface rate
    omegas; smass; log Teff teffl; the zone's moment of inertia hicz; hjm,
    the specific angular momentum (index 0 unused), changed in place.
    None is returned (but nmlpar_mod.wnew is overwritten).
    """
    c = constants_mod
    p = nmlpar_mod
    if not p.linstb:
        wsat = p.wmax
    elif p.wmax > 1.0:
        # Question (see the page): wmax > 1 switches its meaning.
        if ovrtrn_mod.taucz > 1.0:
            wsat = p.wmax / ovrtrn_mod.taucz
        else:
            msg = (f'ERROR IN WIND - TAUCZ NOT DEFINED '
                   f'{p.wmax:12.3E}{ovrtrn_mod.taucz:12.3E}STOPPED')
            print(msg)                     # WRITE(*,911)
            raise SystemExit(msg)
    else:
        wsat = p.wmax
    if ljdot:
        cgrav = math.exp(c.cln * c.cgl)    # noqa: F841 (unused, as in Fortran)
        rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
        rtot = math.exp(c.cln * rl)
        # Question (see the page): a fixed mass-loss rate.
        dmdot = 2.0e-14
        dwtest = ((delts / hicz) * p.constfactor * (dmdot / 1.0e-14)**p.exmd
                  * omegas * (rtot / c.crsun)**p.exr * smass**p.exm
                  * min(omegas, wsat)**(p.exw - 1.0))
        if dwtest > omegas:
            nstep = int(dwtest / omegas) + 1
            dt = delts / float(nstep)
        else:
            nstep = 1
            dt = delts
        ws = omegas
        for n in range(1, nstep + 1):      # loop 100
            k = 0
            w = ws
            wnew0 = ws
            while True:                    # label 5
                k = k + 1
                # Question (see the page): wnew is nmlpar_mod's solid-body
                # rate, overwritten here.
                p.wnew = (ws - (dt / hicz) * p.constfactor
                          * (dmdot / 1.0e-14)**p.exmd
                          * w * (rtot / c.crsun)**p.exr * smass**p.exm
                          * min(w, wsat)**(p.exw - 1.0))
                dw = 2.0 * abs((wnew0 - p.wnew) / (wnew0 + p.wnew))
                if dw > 1.0e-6:
                    w = 0.5 * (ws + p.wnew)
                    wnew0 = p.wnew
                    if k <= 20:
                        continue
                    # Question (see the page): no message when it fails.
                break
            ws = p.wnew
        dm = hstop - hsbot
        fjdom = (omegas - p.wnew) * hicz / dm
        for ii in range(jstart, jend + 1):  # loop 10
            hjm[ii] = hjm[ii] - fjdom
