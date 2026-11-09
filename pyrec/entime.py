"""Python version of entime.f (YREC 7be4133): the time step allowed by the
envelope triangle."""

import math

import numpy as np

from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# SAVE: the work array HMAX(4), filled and never read (Fortran lines
# 47–48).
hmax = np.zeros(4 + 1)


def entime(delts: float, hl: np.ndarray, teffl: float, m: int, edelt: float) -> float:
    """Scale the last step so that the changes of log Teff and log L from
    the previous model stay within the envelope triangle, by a factor
    bounded by atime[13].

    Translates: entime.f#ENTIME (7be4133)

    Arguments: delts, the previous step; hl, the luminosities; teffl,
    log Teff; m, the number of shells; edelt (set only).

    Returns edelt, the step this criterion allows.

    """

    # The triangle's sizes, copied into HMAX, which nothing reads (Fortran
    # lines 73–75).
    hmax[1] = nmlpar_mod.tridt
    hmax[2] = nmlpar_mod.tridl

    # The changes of log Teff and log L since the previous model (Fortran
    # lines 77–82).
    testt = abs(oldmod_mod.tefflo - teffl)
    testl = abs(math.log10(oldmod_mod.hlsuro) - math.log10(hl[m]))

    # The factor, within [1/atime[13], atime[13]], and the step (Fortran
    # lines 86–95).
    # Suspected bug (see the page, htimer.f#b1): a previous step of 0 gives
    # edelt = 0, which HTIMER does not catch.
    fac = testt / nmlpar_mod.tridt
    if testl / nmlpar_mod.tridl > fac:
        fac = testl / nmlpar_mod.tridl
    if fac == 0.0:
        fac = 1.0
    temp = nmlpar_mod.atime[13]
    if fac > temp:
        fac = temp
    if fac < 1.0 / temp:
        fac = 1.0 / temp
    edelt = delts / fac

    return edelt
