"""Python version of cowind.f (YREC 09c3b9b): the change of the surface
convection zone's angular velocity by the older (Kawaler-type) wind law."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod


def cowind(bl, dt, hicz, it, omegas, smass, teffl, wold, wind1, wind2):
    """The change of the surface convection zone's Omega over the step by
    magnetized-wind angular-momentum loss, with a saturation threshold:
    delta Omega = (dt/I_cz) K (Mdot/1e-14)**a (R/R_sun)**b (M/M_sun)**c
    min(Omega, Omega_sat)**(n-1) Omega, for Omega at the start of the
    step (wind1) and at its end (wind2).

    Translates: cowind.f#COWIND (09c3b9b)

    Arguments: bl, log L/L_sun; dt, the time step; hicz, the moment of
    inertia of the surface convection zone; it, the iteration; omegas,
    the surface Omega after the previous iteration; smass, M/M_sun;
    teffl, log Teff; wold, Omega at the start of the step; wind1, passed
    in (with any value); wind2, its value from the previous iteration.

    Returns (wind1, wind2).
    """
    c = constants_mod
    n = nmlpar_mod
    taucz = ovrtrn_mod.taucz
    # MHP 3/09: WMAX > 1 is read as WMAX(sun)*TAUCZ(sun), and the
    # saturation threshold is WMAX/TAUCZ(star).
    # Question (see the page): the test WMAX > 1 works as a switch.
    if n.wmax > 1.0:
        if taucz > 1.0:
            wsat = n.wmax / taucz
        else:
            msg = (f"ERROR IN WIND - TAUCZ NOT DEFINED "
                   f"{n.wmax:12.3E}{taucz:12.3E}STOPPED")
            print(msg)
            raise SystemExit(msg)
    else:
        wsat = n.wmax
    # The total radius of the star (cm).
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    rtot = math.exp(c.cln * rl)
    # The mass-loss rate, in solar masses per year.
    # Question (see the page): fixed, not the model's.
    dmdot = 2.0e-14
    # cc is the Fortran's C (c here names constants_mod).
    cc = (dt / hicz * n.constfactor * (dmdot / 1.0e-14) ** n.exmd
          * (rtot / c.crsun) ** n.exr * smass ** n.exm)
    wp = min(wold, wsat)
    wn = min(omegas, wsat)
    wind1 = cc * wp ** (n.exw - 1.0) * wold
    temp = cc * wn ** (n.exw - 1.0) * omegas
    # Question (see the page): the mean with the previous WIND2.
    if it == 1:
        wind2 = temp
    else:
        wind2 = 0.5 * (wind2 + temp)
    return wind1, wind2
