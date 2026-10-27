"""Python version of op95xtab.f (YREC development, 09c3b9b): the OPAL 1995
opacity table at a fixed X, from the table at the model's Z."""

import numpy as np

from pyrec import llot95_mod
from pyrec import llot95a_mod
from pyrec import llot95e_mod
from pyrec.intrp2 import intrp2

NUMT = 70
NUMD = 19
NUMX = 10

# OP95XTAB has SAVE: its locals are kept between calls.
a = np.zeros(5)
b = np.zeros(5)
idx = 0


def op95xtab(x):
    """Build CAPPAE, the OPAL95 table at the hydrogen abundance x (the
    surface X) and the model's Z, by four-point Lagrange interpolation in X
    between the tables of CAPPA (MHP 7/98).

    Translates: op95xtab.f#OP95XTAB (09c3b9b)

    Arguments: x, the desired hydrogen abundance.

    Returns nothing; sets llot95e_mod.cappae and llot95e_mod.xop95.
    """
    global idx
    xxg = llot95a_mod.xxg
    ztab = llot95_mod.ztab
    llot95e_mod.xop95 = x
    # Find the 4 nearest tables in X to the desired value.
    if x <= 0.8:
        # No missing X tables at high Z to worry about.
        for i in range(3, NUMX):                     # I = 3, NUMX-1
            if xxg[i] >= x:
                idx = i - 2
                break
        else:
            idx = NUMX - 3
        # label 10
        for i in range(1, 5):
            a[i] = xxg[idx + i - 1]
    elif ztab <= 0.04:
        # Question (see the page): table 9 at Z = 0.04.
        # High X tables present at low Z.
        for i in range(3, NUMX):
            if xxg[i] >= x:
                idx = i - 2
                break
        else:
            idx = NUMX - 3
        # label 20
        for i in range(1, 5):
            a[i] = xxg[idx + i - 1]
    elif ztab >= 0.1:
        # Question (see the page): tables 6-9, though table 9 is unfilled.
        idx = 6
        for i in range(1, 5):
            a[i] = xxg[idx + i - 1]
    else:
        # Z between 0.04 and 0.1: tables 1-8 and 10 exist; X > 0.8, so use
        # 6-8 and 10. Question (see the page): the sentinel IDX = 11.
        idx = 11
        a[1] = xxg[6]
        a[2] = xxg[7]
        a[3] = xxg[8]
        a[4] = xxg[10]
    # Interpolation factors in X (the Fortran comment says Z: Question).
    intrp2(a, b, x)
    # Indices of the 4 compositions.
    if idx < 10:
        i1 = idx
        i2 = idx + 1
        i3 = idx + 2
        i4 = idx + 3
    else:
        i1 = 6
        i2 = 7
        i3 = 8
        i4 = 10
    cappa = llot95_mod.cappa
    cappae = llot95e_mod.cappae
    for j in range(1, NUMT + 1):
        for k in range(1, NUMD + 1):
            cappae[j, k] = (b[1] * cappa[i1, j, k] + b[2] * cappa[i2, j, k]
                            + b[3] * cappa[i3, j, k]
                            + b[4] * cappa[i4, j, k])
