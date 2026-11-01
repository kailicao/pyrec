"""Python version of deutrate.f (YREC a8a5950): the rate of deuterium
burning out of equilibrium at one point."""

import math

from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import ovrtrn_mod

# DATA C21/5.240358D-8/: seconds per Gyr over Avogadro's number, probably.
C21 = 5.240358e-8
# DATA AMDEUT/2.013553D0/: the bare deuteron mass in amu, the same value as
# ANUC(3) in engeb.f (keep the two in step).
AMDEUT = 2.013553


def deutrate(dl, tl, x, i, itlvl):
    """Rate of 2H(p,gamma)3He at point i, without screening, per Gyr,
    capped in the surface convection zone so that deuterium burns no
    faster than the convective overturn time.

    Translates: deutrate.f#DEUTRATE (a8a5950)

    Arguments: dl, tl (log rho, log T); x, the hydrogen mass fraction;
    i, the point; itlvl, the level of iteration.

    Returns nothing: sets deuter_mod.drate[i], and at the first level of
    iteration also deuter_mod.drate0[i].
    """
    c = constants_mod
    rho = math.exp(c.cln * dl)
    t9 = math.exp(c.cln * (tl - 9.0))
    t9p13 = t9**c.cc13
    t9p23 = t9p13**2
    t9m13 = 1.0 / t9p13
    t9m23 = t9m13**2
    t9m1 = 1.0 / t9
    z = -3.72 * t9m13
    tfacdeut = 1.0 + 0.112 * t9p13 + 3.38 * t9p23 + 2.65 * t9
    # Avogadro's number over the deuteron's mass in amu.
    rdeut = (rho * 2.240e3 * t9m23 * math.exp(z) * tfacdeut * c.cavog
             / AMDEUT)
    # Cap in the surface convection zone by the overturn time.
    if i >= deuter_mod.jcz and ovrtrn_mod.taucz > 1.0:
        rdeutmax = c.cavog / AMDEUT / ovrtrn_mod.taucz
        rdeut2 = rdeut * x
        if rdeut2 > rdeutmax:
            if x > 1.0e-6:
                rdeut = rdeutmax / x
    deuter_mod.drate[i] = x * rdeut * C21
    if itlvl == 1:
        deuter_mod.drate0[i] = deuter_mod.drate[i]
