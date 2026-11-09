"""Python version of rotmix.f (YREC 7be4133): nuclear burning and
gravitational settling for one step of the rotational-mixing loop."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod
from pyrec.grsett import grsett
from pyrec.kemcom import kemcom
from pyrec.microdiff import microdiff
from pyrec.rates import rates

JSON = dimensions_mod.json

# Local arrays (SAVE): the rates from rates, the settling's
# pressure gradient and flags, and the saved SDEL(2,.).
hr1 = np.zeros(JSON + 1)
hr2 = np.zeros(JSON + 1)
hr3 = np.zeros(JSON + 1)
hr4 = np.zeros(JSON + 1)
hr5 = np.zeros(JSON + 1)
hr6 = np.zeros(JSON + 1)
hr7 = np.zeros(JSON + 1)
hr8 = np.zeros(JSON + 1)
hr9 = np.zeros(JSON + 1)
hr10 = np.zeros(JSON + 1)
hr11 = np.zeros(JSON + 1)
hr12 = np.zeros(JSON + 1)
hr13 = np.zeros(JSON + 1)
hf1 = np.zeros(JSON + 1)
hf2 = np.zeros(JSON + 1)
hqpr = np.zeros(JSON + 1)
temp = np.zeros(JSON + 1)
lcz = np.zeros(JSON + 1, dtype=bool)


def rotmix(delts: float, hcomp: np.ndarray, hs2: np.ndarray, ht: np.ndarray, m: int, mrzone: np.ndarray,
           mxzone: np.ndarray, nrzone: int, nzone: int, hstot: float, hd: np.ndarray, hs: np.ndarray, hr: np.ndarray,
           hp: np.ndarray, lc: np.ndarray, hs1: np.ndarray) -> None:
    """Burn each hot radiative shell and each convection zone with kemcom
    (which also homogenizes the zones), apply the settling of helium
    (microdiff or grsett, in substeps limited by the settling time at the
    outer point), and clip and renormalize the composition.

    Translates: rotmix.f#ROTMIX (7be4133)

    Arguments: the step delts (s); the composition hcomp[species, shell]
    (changed in place); shell masses hs2; log T ht; m shells; the zone
    lists mrzone[k, 1..2] (radiative) and mxzone[k, 1..2] (convective),
    with their counts nrzone, nzone; for the settling, log total mass
    hstot, hd, hs, hr, hp (logs), the convection flags lc and the masses
    hs1 (g). None is set (nmlpar_mod.ldify may be switched off).

    """

    c = constants_mod
    n = nmlpar_mod

    # NSPEC is the number of species being tracked (Fortran lines 74–78).
    if n.lexcom:
        nspec = 15
    else:
        nspec = 11

    # The step in years (Fortran line 80).
    ddage = delts / c.csecyr

    # The burning rates on the model passed in (hd, ht and the sub-step's
    # own composition), as mix does (MHP via Claude, 10/29/26) (Fortran lines 87–110).
    for ii in range(1, m + 1):
        if ht[ii] <= n.tcut[1]:
            for arr in (hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10,
                        hr11, hr12, hr13, hf1, hf2):
                arr[ii] = 0.0
        else:
            rates(hd[ii], ht[ii], hcomp[1, ii], hcomp[2, ii], hcomp[4, ii],
                  hcomp[5, ii], hcomp[6, ii], hcomp[7, ii], hcomp[9, ii],
                  hcomp[11, ii], ii, hr1, hr2, hr3, hr4, hr5, hr6, hr7,
                  hr8, hr9, hr10, hr11, hr12, hr13, hf1, hf2)

    # Radiative zones, shell by shell; stop at the first shell below the
    # nuclear cutoff (GOTO 45 leaves both loops) (Fortran lines 118–130).
    # Question (see the page): any later radiative zone is skipped.
    cool = False
    for k in range(1, nrzone + 1):
        for j in range(mrzone[k, 1], mrzone[k, 2] + 1):
            if ht[j] <= n.tcut[1]:
                cool = True
                break
            ibegin = j
            iend = j
            kemcom(ht, ibegin, iend, hr1, hr2, hr3, hr4, hr5, hr6, hr7,
                   hr8, hr10, hr11, hr12, hf1, hs2, hcomp, ddage)
        if cool:
            break

    # label 45. Convection zones, each as a unit (kemcom homogenizes it) (Fortran lines 135–142).
    for k in range(1, nzone + 1):
        ibegin = mxzone[k, 1]
        iend = mxzone[k, 2]
        kemcom(ht, ibegin, iend, hr1, hr2, hr3, hr4, hr5, hr6, hr7,
               hr8, hr10, hr11, hr12, hf1, hs2, hcomp, ddage)

    # Microscopic diffusion of helium (Fortran lines 151–237).
    if n.ldify:
        _settle(delts, hcomp, ht, m, mrzone, mxzone, nrzone, nzone,
                hstot, hd, hs, hr, hp, lc, hs1)
    # End of gravitational settling (label 170).

    # Renormalize the composition (small negative abundances...) (Fortran lines 242–249).
    for i in range(1, m + 1):
        for j in range(1, nspec + 1):
            hcomp[j, i] = max(hcomp[j, i], 0.0)
            hcomp[j, i] = min(hcomp[j, i], 1.0)
        hcomp[3, i] = min(hcomp[3, i], 1.0 - hcomp[1, i] - hcomp[4, i])
        hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]


def _settle(delts: float, hcomp: np.ndarray, ht: np.ndarray, m: int, mrzone: np.ndarray, mxzone: np.ndarray,
            nrzone: int, nzone: int, hstot: float, hd: np.ndarray, hs: np.ndarray, hr: np.ndarray, hp: np.ndarray,
            lc: np.ndarray, hs1: np.ndarray) -> None:
    """The settling block inside IF(LDIFY); each early return is a
    GOTO 170.

    """

    c = constants_mod
    n = nmlpar_mod
    sdel = scrtch_mod.sdel

    # Settling stops for good once the central X is below XMIN (Fortran lines 152–155).
    if hcomp[1, 1] < n.xmin:
        n.ldify = False
        return

    # Locate the outer boundary (Fortran lines 163–180).
    if not lc[m]:
        text = ("NO SURFACE CZ - DIFFUSION NOT MEANINGFUL\n"
                "STOPPED IN SUBROUTINE MIX")   # format 911
        # Suspected bug (see the page): the message names MIX.
        print(" " + text)
        iounits_mod.write(iounits_mod.ishort, " " + text)
        raise SystemExit(1)    # STOP 1, the run-time error code
    # Exit if overshoot makes the star fully convective (MHP 7/92).
    if mxzone[nzone, 1] <= 2 and mxzone[nzone, 2] == m:
        return
    jmax = mrzone[nrzone, 2]
    for j in range(jmax, 0, -1):
        if hcomp[2, j] > n.ymin:
            break                                   # go to 150
    else:
        return                  # Y < YMIN for the whole star
    # label 150. Suspected bug (see the page): the loop above moves J, not
    # JMAX, so JMAX stays the top of the last radiative zone and the first
    # shell with Y > YMIN is never used.

    # The pressure gradient; SDEL(2,.) holds DELM meanwhile (Fortran lines 181–186).
    stot = math.exp(c.cln * hstot)
    for i in range(1, m + 1):
        temp[i] = sdel[2, i]
        sdel[2, i] = mdphy_mod.delm[i]
        hqpr[i] = -math.exp(c.cln * (hd[i] + c.cgl + hs[i]
                                     - 2.0 * hr[i] - hp[i]))

    # Only check for timestep cutting if JMAX is large (BC 6/92) (Fortran lines 189–213).
    if jmax >= 2:
        # FM is the mass fraction above the outer point.
        fm = (stot - hs1[jmax]) / stot
        # Settling time of helium there (Michaud et al. 1984).
        tscale = 4.348e21 * c.csecyr * fm / math.exp(c.cln * 1.5 * ht[jmax])
        dtmax = n.dt_gs * tscale
        nstep = int(delts / dtmax)
        if math.fmod(delts, dtmax) != 0.0 or nstep == 0:
            nstep = nstep + 1
        dt = delts / float(nstep)
    else:
        nstep = 1
        dt = delts

    # The convection flags from the zone lists (Fortran lines 214–223).
    for j in range(1, nrzone + 1):
        for i in range(mrzone[j, 1], mrzone[j, 2] + 1):
            lcz[i] = False
    for j in range(1, nzone + 1):
        for i in range(mxzone[j, 1], mxzone[j, 2] + 1):
            lcz[i] = True

    # The settling substeps (Fortran lines 224–232).
    for ns in range(1, nstep + 1):
        # LNEWDIF: the new routines in microdiff; else grsett.
        if n.lnewdif:
            # Since 7be4133 microdiff and grsett work on copies of dt
            # and stot and leave them unchanged.
            microdiff(dt, hcomp, hqpr, hr, hd, hs1, ht, lcz, m, stot)
        else:
            grsett(dt, hcomp, hqpr, hr, hd, hs1, ht, lcz, m, stot)

    # SDEL(2,.) restored (Fortran lines 233–235).
    for i in range(1, m + 1):
        sdel[2, i] = temp[i]
