"""Python version of ptime.f (YREC 09c3b9b): the time step limited by the
structure's change since the previous model."""

import numpy as np

from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# Locals kept between calls (SAVE).
hmax = np.zeros(5)
nhmax = np.zeros(5, dtype=int)


def ptime(delts, hl, hp, hr, ht, m, pdelt):
    """Propose a time step from the largest changes, over all shells, of
    log P, log T, log r and the relative luminosity since the previous
    model (OLDMOD_MOD), against the limits ATIME(8)-ATIME(11), with the
    factor bounded by ATIME(13).

    Translates: ptime.f#PTIME (09c3b9b)

    Arguments: delts, the last step; hl, hp, hr, ht, the current model;
    m, the number of shells; pdelt (set). Returns pdelt, the proposed step.
    """
    o = oldmod_mod
    atime = nmlpar_mod.atime
    hmax[1] = abs(o.hpo[1] - hp[1])
    hmax[2] = abs(o.hto[1] - ht[1])
    hmax[3] = abs(o.hro[1] - hr[1])
    if hl[1] + hl[2] > 0.0:
        hmax[4] = abs((o.hlo[1] - hl[1]) * 2.0 / (hl[2] + hl[1]))
    else:
        hmax[4] = 0.0
    # Question (see the page): NHMAX is filled but never used.
    for i in range(2, m + 1):
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
    temp = atime[13]
    if fac > temp:
        fac = temp
    if fac < 1.0 / temp:
        fac = 1.0 / temp
    pdelt = delts / fac
    return pdelt
