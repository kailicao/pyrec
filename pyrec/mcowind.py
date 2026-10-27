"""Python version of mcowind.f (YREC 09c3b9b): the change in the surface
convection zone's angular velocity from the newer magnetized-wind law."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.amcalc import amcalc
from pyrec.cowind import cowind

# GL is a SAVEd local that is never assigned (the stellar terms use CGL).
# Question (see the page): undefined G in the solar term; 0.0 here, as a
# compiler that zeroes static storage would give, so 10**GL = 1.
GL = 0.0


def mcowind(bl, dt, hicz, it, omegas, smass, teffl, wold, wind1, wind2):
    """Return the change in the surface convection zone's angular velocity
    over a diffusion step from magnetized-wind loss: a saturated power law
    in Omega, with optional Rossby scaling of the saturation and the
    centrifugal correction of Matt et al. (2012). With LMWIND false it
    hands over to the older law, COWIND.

    Translates: mcowind.f#MCOWIND (09c3b9b)

    Arguments: bl (log L/Lsun), dt (the diffusion step, s), hicz (moment of
    inertia of the convection zone), it (the iteration), omegas (surface
    Omega now), smass (M/Msun), teffl (log Teff), wold (surface Omega at
    the start of the step); wind1 is output only; wind2 is the running
    mean from the previous iteration (used for it > 1).

    Returns (wind1, wind2): the changes in Omega computed with Omega at the
    start and at the end of the step.
    """
    c = constants_mod
    n = nmlpar_mod
    if not n.lmwind:
        wind1, wind2 = cowind(bl, dt, hicz, it, omegas, smass,
                              teffl, wold, wind1, wind2)
        return wind1, wind2  # label 9999
    wfirst = wold
    wnow = omegas
    if n.lrossby:
        # Turnover time interpolated within the step; saturation scaled.
        tcz = (ovrtrn_mod.taucz0
               + ovrtrn_mod.fracstep * (ovrtrn_mod.taucz - ovrtrn_mod.taucz0))
        wsat = n.wmax * n.pmmsoltau / tcz
        if n.lbscale:
            # Question (see the page): only the rates in the min() are scaled.
            wfirst = wfirst * tcz / n.pmmsoltau
            wnow = wnow * tcz / n.pmmsoltau
            wsat = wsat * tcz / n.pmmsoltau
    else:
        wsat = n.wmax
    # New calculation: AMCALC sets the structure factor.
    structfactor = amcalc(smass, bl, teffl, 0.0)
    # The new wind coefficient.
    cc = (dt / hicz) * n.constfactor * structfactor
    # Centrifugal correction (Matt et al. 2012).
    # Question (see the page): GL is never set; CGL was probably meant.
    fsun = 0.5 * n.pmmsolw**2 * c.crsun**3 / math.exp(c.cln * GL) / c.cmsun
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    fcorr1 = 0.5 * wold**2 * math.exp(c.cln * (3.0 * rl - c.cgl)) / smass / c.cmsun
    fcorr2 = 0.5 * omegas**2 * math.exp(c.cln * (3.0 * rl - c.cgl)) / smass / c.cmsun
    fcen1 = ((c.c_2**2 + fsun) / (c.c_2**2 + fcorr1))**n.excen
    fcen2 = ((c.c_2**2 + fsun) / (c.c_2**2 + fcorr2))**n.excen
    wp = min(wfirst, wsat)
    wn = min(wnow, wsat)
    wind1 = cc * wp**(n.exw - 1.0) * wold * fcen1
    temp = cc * wn**(n.exw - 1.0) * omegas * fcen2
    if it == 1:
        wind2 = temp
    else:
        # Question (see the page): a running half-and-half mean.
        wind2 = 0.5 * (wind2 + temp)
    return wind1, wind2
