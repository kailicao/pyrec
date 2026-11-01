"""Python version of mixcz.f (YREC 09c3b9b): every convection zone
homogenized, then light clean-ups of the composition."""

import numpy as np

from pyrec import comp_mod
from pyrec import nmlpar_mod

# The Fortran's local arrays (SAVE): up to 12 zones as pairs (first, last
# shell) in mxzone[1..24], and the sums.
mxzone = np.zeros(25, dtype=int)
xsum = np.zeros(16)


def mixcz(hcomp, hs2, lcz, m):
    """Mix each convection zone to its mass-weighted mean composition
    (species 1..11, or 1..15 with LEXCOM), then at every shell floor X at
    zero, cap Z at 1 - X and floor O16 (row 9) at 0.99e-3 (Z - ZENVM).

    Translates: mixcz.f#MIXCZ (09c3b9b)

    Arguments: the model hcomp[k, j], the shell masses hs2 and the
    convective flags lcz, with m shells (arrays with index 0 unused; lcz
    needs room for the sentinel lcz[m + 1]).

    Sets in place: hcomp, and lcz[m + 1] = False. Returns nothing.
    """
    # The standard convection zones.
    j = 1
    lcc = False
    mp1 = m + 1
    # Question (see the page): the sentinel is written past the model.
    lcz[mp1] = False
    i1 = 0
    full = False
    for i in range(1, mp1 + 1):                          # loop 11
        if lcz[i]:
            if not lcc:
                # Start of a convection zone.
                lcc = True
                i1 = i
        elif lcc:
            # End of a convection zone (label 10).
            lcc = False
            mxzone[j] = i1
            mxzone[j + 1] = i - 1
            j = j + 2
            # Question (see the page): more than 12 zones are ignored.
            if j >= 24:
                full = True
                break                                    # GOTO 12
    if not full:
        mxzone[j] = 0
    # label 12
    iend = 11
    if nmlpar_mod.lexcom:
        iend = 15
    # Mix all convective zones.
    for j in range(1, 25, 2):                            # loop 100
        if mxzone[j] <= 0:
            break                                        # GOTO 110
        i1 = mxzone[j]
        i2 = min(m, mxzone[j + 1])
        # One shell is enough only at the centre.
        if i1 != 1 and i1 >= i2:
            continue
        wtsum = 0.0
        for i in range(1, iend + 1):
            xsum[i] = 0.0
        for jj in range(i1, i2 + 1):
            wtsum = wtsum + hs2[jj]
            for i in range(1, iend + 1):
                xsum[i] = xsum[i] + hcomp[i, jj] * hs2[jj]
        for i in range(1, iend + 1):
            xsum[i] = xsum[i] / wtsum
        for jj in range(i1, i2 + 1):
            for i in range(1, iend + 1):
                hcomp[i, jj] = xsum[i]
    # label 110: renormalize the composition if necessary.
    # Question (see the page): Y is not recomputed after these.
    for i in range(1, m + 1):                            # loop 140
        hcomp[1, i] = max(hcomp[1, i], 0.0)
        hcomp[3, i] = min(hcomp[3, i], 1.0 - hcomp[1, i])
        # Question (see the page): the oxygen floor.
        hcomp[9, i] = max(hcomp[9, i],
                          0.99e-3 * (hcomp[3, i] - comp_mod.zenvm))
    # The convective turnover time is not computed here: it comes from
    # OVRTRN_MOD (G. Somers 3/17).
