"""Python version of microdiff_setup.f (YREC 09c3b9b): the settling
region and the conversion to Bahcall and Loeb's units, for MICRODIFF."""

import math

from pyrec import confac_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod


def microdiff_setup(dt, hqpr, hr, hd, hs1, ht, lc, m, stot, hcomp, hru,
                    htu, ibegin, iend, lallcz, oden, otem):
    """Find the region that takes part in settling (IBEGIN to IEND) and
    convert the model to Bahcall and Loeb's units; skip settling
    (lallcz = True) for a fully convective model, or one with X <= XMIN or
    Y <= YMIN everywhere.

    Translates: microdiff_setup.f#MICRODIFF_SETUP (09c3b9b)

    Arguments: the step dt (s); the model's hqpr (dlnP/dr), hr (log r),
    hd (log rho), hs1 (mass, g), ht (log T), lc (convective flags) and
    hcomp[k, j], with m shells; the total mass stot (g). ibegin, iend and
    lallcz are passed in with any value. Arrays have index 0 unused.

    Sets in place: hru, htu (unlogged r and T in the new units), oden, otem
    (unlogged rho and T, cgs), and converts hs1 and hqpr; writes the
    conversion factors to confac_mod.

    Returns (dt, stot, ibegin, iend, lallcz): dt and stot in the new units
    (unchanged when settling is skipped), the region, and the skip flag.
    """
    cln = constants_mod.cln
    ishort = iounits_mod.ishort
    crsun_bah = 6.9598e10
    csecyr_bah = 3.1558e7
    lallcz = False
    # Check for a convective core.
    if lc[1]:
        i = 2
        while i <= m:                     # loop 10
            if not lc[i]:
                break                     # GOTO 20
            i += 1
        else:
            lallcz = True
            iounits_mod.write(ishort, "  FULLY CONVECTIVE MODEL - NO SETTLING")
            return dt, stot, ibegin, iend, lallcz      # GOTO 9999
        # Compute overshoot (to be added).
        ibegin = i - 1
    else:
        ibegin = 1
    # MHP 6/90: check for a hydrogen-exhausted core.
    for i in range(ibegin, m + 1):        # loop 23
        if hcomp[1, i] > nmlpar_mod.xmin:
            break                         # GOTO 25
    else:
        # Question (see the page): "HYDROGEN-FREE MODEL" comment, both exits.
        iounits_mod.write(ishort, f" X BELOW {nmlpar_mod.xmin:9.6f}"
                          " IN WHOLE MODEL-NO SETTLING")
        lallcz = True
        return dt, stot, ibegin, iend, lallcz
    ibegin = i                            # label 25
    # Check for a convective envelope.
    if lc[m]:
        i = m - 1
        while i >= 2:                     # loop 30
            if not lc[i]:
                break                     # GOTO 40
            i -= 1
        # On completion the Fortran's I is 1.
        iend = i + 1                      # label 40
    else:
        iend = m
    # Check for a helium-exhausted surface: the outer point is the last
    # one with Y > YMIN.
    for i in range(iend, 0, -1):          # loop 45
        if hcomp[2, i] > nmlpar_mod.ymin:
            break                         # GOTO 47
    else:
        iounits_mod.write(ishort, f" Y BELOW {nmlpar_mod.ymin:9.6f}"
                          " IN WHOLE MODEL-NO SETTLING")
        lallcz = True
        return dt, stot, ibegin, iend, lallcz
    iend = i                              # label 47
    # Question (see the page): no check that iend > ibegin.
    # Conversion factors (the time unit includes a factor of 2.2 from
    # ln Lambda).
    confac_mod.con_rad = 1.0 / crsun_bah
    confac_mod.con_mass = 1.0e-2 * confac_mod.con_rad**3
    confac_mod.con_temp = 1.0e-7
    confac_mod.con_time = 2.7e13 * csecyr_bah
    for i in range(1, m + 1):             # loop 50
        hru[i] = math.exp(cln * hr[i]) * confac_mod.con_rad
        htu[i] = math.exp(cln * ht[i]) * confac_mod.con_temp
        hs1[i] = hs1[i] * confac_mod.con_mass
        hqpr[i] = hqpr[i] / confac_mod.con_rad
    dt = dt / confac_mod.con_time
    stot = stot * confac_mod.con_mass
    # Rho and T for the later transformation to the equal grid.
    for i in range(1, m + 1):
        oden[i] = math.exp(cln * hd[i])
        otem[i] = math.exp(cln * ht[i])
    return dt, stot, ibegin, iend, lallcz
