"""Python version of entime.f (YREC 09c3b9b): the time step allowed by the
envelope triangle."""

import math

import numpy as np

from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# SAVE: the work array HMAX(4), filled and never read.
hmax = np.zeros(4 + 1)


def entime(delts, hl, teffl, m, edelt):
    """Scale the last step so that the changes of log Teff and log L from
    the previous model stay within the envelope triangle, by a factor
    bounded by atime[13].

    Translates: entime.f#ENTIME (09c3b9b)

    Arguments: delts, the previous step; hl, the luminosities; teffl,
    log Teff; m, the number of shells; edelt (set only).

    Returns edelt, the step this criterion allows.
    """
    # Question (see the page): hmax is filled but never used.
    hmax[1] = nmlpar_mod.tridt
    hmax[2] = nmlpar_mod.tridl

    testt = abs(oldmod_mod.tefflo - teffl)
    testl = abs(math.log10(oldmod_mod.hlo[oldmod_mod.mo]) - math.log10(hl[m]))

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
