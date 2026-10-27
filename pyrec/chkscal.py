"""Python version of chkscal.f (YREC 09c3b9b): calibrating a star's track
to a target radius and luminosity."""

import math

from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod

# CHKSCAL has SAVE (formerly part of COMMON/CALSTAR/): the previous model's
# L and age, and L and age at R*; BLRP is set in one run and read in the
# next.
bli = 0.0
agei = 0.0
blr = 0.0
ager = 0.0
blrp = 0.0


def chkscal(bl, teffl, dage, nk, alri, lstar, lpassr):
    """Calibrate X (thus Y) so that a track passes through a target
    (R*, L*): detect the crossing of R*, interpolate L and the age there,
    and set X for the next run by a secant step, or set the final run's
    stopping age.

    Translates: chkscal.f#CHKSCAL (09c3b9b)

    Arguments: bl, log L/L_sun; teffl, log Teff; dage, the age (Gyr); nk,
    the kind card; alri, the previous model's R/R_sun (despite the comments,
    not a log); lstar, the star found; lpassr, R* just passed.

    Returns (alri, lstar, lpassr). The new X goes to runloop_mod.rescal and
    nmlpar_mod.xenv0a; the stopping age to runloop_mod.lendag and
    nmlpar_mod.endage.
    """
    global bli, agei, blr, ager, blrp
    c = constants_mod
    nml = nmlpar_mod
    rescal = runloop_mod.rescal
    itrack = iounits_mod.itrack

    # Has the star passed R*? If not, store L and age and return.
    # Question (see the page): linear or log radius; and the previous
    # model across runs.
    lpassr = False
    teff = 10.0**teffl
    alr = math.sqrt((10.0**bl) * c.clsun / (c.c4pi * c.csig)) / (
        teff * teff * c.crsun)
    if alr > alri:
        if not (alr > nml.sr and alri <= nml.sr):
            bli = bl
            agei = dage
            alri = alr
            return alri, lstar, lpassr
    else:
        if not (alr < nml.sr and alri >= nml.sr):
            bli = bl
            agei = dage
            alri = alr
            return alri, lstar, lpassr

    # Question (see the page): the first crossing only.
    print("  Just passed R*")
    iounits_mod.write(itrack, " #Just passed R*")
    lpassr = True

    # Interpolate L and the age at R*, linearly in R.
    dldr = (bl - bli) / (alr - alri)
    blr = bl + dldr * (nml.sr - alr)
    dadr = (dage - agei) / (alr - alri)
    ager = dage + dadr * (nml.sr - alr)
    print(f"  X, LogL/Lsun at R* = {rescal[2, nk - 1]} {blr}")
    iounits_mod.write(itrack, f" #X, LogL/Lsun at R* = {rescal[2, nk - 1]}"
                              f" {blr}")
    # Question (see the page): two tolerances.
    if abs(10.0**blr - nml.xls) <= nml.xlstol:
        # Found: one more run, with the same X, stopping at that age.
        lstar = True
        runloop_mod.lendag[nk + 1] = True
        nml.endage[nk + 1] = ager * 1.0e9
        runloop_mod.lendag[nk + 2] = True
        nml.endage[nk + 2] = ager * 1.0e9
        rescal[2, nk + 1] = rescal[2, nk - 1]
        nml.xenv0a[nk + 1] = rescal[2, nk + 1]
        nml.xenv0a[nk + 2] = nml.xenv0a[nk + 1]
        print(f"  Have hit R* & L*, prepare final run to age: {ager}")
        iounits_mod.write(itrack, " #Have hit R* & L*, prepare final run "
                                  f"to age: {ager}")
        return alri, lstar, lpassr
    if nk == 2:
        # First time through: save L at R*, lower X by 0.01.
        blrp = blr
        xp = rescal[2, nk - 1]
        x = xp - 0.01
        rescal[2, nk + 1] = x
        nml.xenv0a[nk + 1] = x
        nml.xenv0a[nk + 2] = x
        print(f"  NK=2, Y=Y+0.01, Setup next run, X= {x}")
        iounits_mod.write(itrack, f" #NK=2, Y=Y+0.01, Setup next run, X= {x}")
        return alri, lstar, lpassr
    # NK = 4, 6, 8, ...: a secant step in X against log L.
    x = rescal[2, nk - 1]
    xp = rescal[2, nk - 3]
    dxdl = (x - xp) / (blr - blrp)
    x = dxdl * (math.log10(nml.xls) - blr) + x
    print(f"  Setup next run, NK, X = {nk} {x}")
    iounits_mod.write(itrack, f"  Setup next run, NK, X = {nk} {x}")
    blrp = blr
    rescal[2, nk + 1] = x
    nml.xenv0a[nk + 1] = x
    nml.xenv0a[nk + 2] = x
    return alri, lstar, lpassr
