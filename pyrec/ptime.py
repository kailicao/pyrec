"""Python version of ptime.f (YREC 7be4133): the time step limited by the
structure's change since the previous model."""

import numpy as np

from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# Locals kept between calls (SAVE) (Fortran lines 46–47).
hmax = np.zeros(5)
nhmax = np.zeros(5, dtype=int)


def ptime(delts: float, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray, ht: np.ndarray, m: int, pdelt: float) -> float:
    """Propose a time step from the largest changes, over all shells, of
    log P, log T, log r and the relative luminosity since the previous
    model (OLDMOD_MOD), against the limits ATIME(8)-ATIME(11), with the
    factor bounded by ATIME(13).

    Translates: ptime.f#PTIME (7be4133)

    Arguments: delts, the last step; hl, hp, hr, ht, the current model;
    m, the number of shells; pdelt (set). Returns pdelt, the proposed step.

    """

    o = oldmod_mod
    atime = nmlpar_mod.atime

    # The changes at the first shell (Fortran lines 80–90).
    hmax[1] = abs(o.hpo[1] - hp[1])
    hmax[2] = abs(o.hto[1] - ht[1])
    hmax[3] = abs(o.hro[1] - hr[1])
    if hl[1] + hl[2] > 0.0:
        hmax[4] = abs((o.hlo[1] - hl[1]) * 2.0 / (hl[2] + hl[1]))
    else:
        hmax[4] = 0.0

    # The largest changes over all shells; NHMAX, where each occurs, is
    # read nowhere (Fortran lines 91–116).
    for i in range(2, m + 1):  # loop 40
        testp = abs(o.hpo[i] - hp[i])
        if hmax[1] <= testp:
            hmax[1] = testp
            nhmax[1] = i
        testt = abs(o.hto[i] - ht[i])
        if hmax[2] <= testt:
            hmax[2] = testt
            nhmax[2] = i
        testr = abs(o.hro[i] - hr[i])
        if hmax[3] <= testr:
            hmax[3] = testr
            nhmax[3] = i
        if hl[i] + hl[i - 1] > 0.0:
            testl = abs((o.hlo[i] - hl[i]) * 2.0 / (hl[i] + hl[i - 1]))
        else:
            testl = 0.0
        if hmax[4] <= testl:
            hmax[4] = testl
            nhmax[4] = i

    # The factor that brings every change within its limit (Fortran lines
    # 119–127).
    fac = hmax[1] / atime[9]
    if hmax[2] / atime[8] > fac:
        fac = hmax[2] / atime[8]
    if hmax[3] / atime[10] > fac:
        fac = hmax[3] / atime[10]
    # Question (see the page): ATIME(11) times ln 10, a limit on d log L.
    if hmax[4] / (atime[11] * 2.3026) > fac:
        fac = hmax[4] / (atime[11] * 2.3026)
    if fac == 0.0:
        fac = 1.0

    # The factor kept within [1/atime[13], atime[13]] (Fortran lines
    # 132–135).
    temp = atime[13]
    if fac > temp:
        fac = temp
    if fac < 1.0 / temp:
        fac = 1.0 / temp
    pdelt = delts / fac

    return pdelt
