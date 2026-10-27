"""Python version of op95ztab.f (YREC 09c3b9b): the OPAL95 opacity tables
at the model's fixed Z, from the full set."""

import numpy as np

from pyrec import llot95_mod
from pyrec import llot95a_mod
from pyrec.intrp2 import intrp2

NUMT = 70
NUMD = 19
NUMZ = 13

# Locals kept between calls (SAVE): IDZ2 and BB matter if label 30 is
# reached without setting them (question 2).
a = np.zeros(5)
b = np.zeros(5)
bb = np.zeros(5)
idz = 0
idz2 = 0


def _fill_table(itab, idzz, w):
    """CAPPA(itab, j, k) = sum_n w(n) * CAPPA2(NZ(idzz+n-1)+itab, j, k),
    summed in the Fortran's order."""
    m = llot95a_mod
    i1 = m.nz[idzz] + itab
    i2 = m.nz[idzz + 1] + itab
    i3 = m.nz[idzz + 2] + itab
    i4 = m.nz[idzz + 3] + itab
    for j in range(1, NUMT + 1):
        for k in range(1, NUMD + 1):
            llot95_mod.cappa[itab, j, k] = (
                w[1] * m.cappa2[i1, j, k] + w[2] * m.cappa2[i2, j, k]
                + w[3] * m.cappa2[i3, j, k] + w[4] * m.cappa2[i4, j, k])


def op95ztab(z):
    """Interpolate the full OPAL95 set (LLOT95A_MOD) in Z, with four-point
    Lagrange weights, to the ten X tables at the desired Z (CAPPA and ZTAB
    in LLOT95_MOD). Tables 1-8 always; table 9 (X = 0.95) only for
    Z < 0.04; table 10 (X = 1 - Z) only for Z < 0.1. Stops the run if Z
    lies above the twelfth grid value, 0.08.

    Translates: op95ztab.f#OP95ZTAB (09c3b9b)

    Arguments: z, the desired metallicity. Returns nothing.
    """
    global idz, idz2
    zzg = llot95a_mod.zzg
    llot95_mod.ztab = z
    # Find the four Z values around z; the loop stops at NUMZ-1 (question 1).
    for i in range(3, NUMZ):
        if zzg[i] >= z:
            idz = i - 2
            break
    else:
        print(f" DESIRED Z{z:10.6f}OUTSIDE OP95 TABLE RANGE\n RUN STOPPED")
        raise SystemExit()
    # label 10
    for i in range(1, 5):
        a[i] = zzg[idz + i - 1]
    intrp2(a, b, z)
    for i in range(1, 9):
        _fill_table(i, idz, b)
    # Table 9 (X = 0.95) only for Z < 0.04.
    if not z >= 0.04:
        if zzg[idz + 3] > 0.04:
            idz2 = NUMZ - 6     # Question (see the page): fixed grid position
            for i in range(1, 5):
                a[i] = zzg[idz2 + i - 1]
            intrp2(a, bb, z)
        else:
            idz2 = idz
            for i in range(1, 5):
                bb[i] = b[i]
        _fill_table(9, idz2, bb)
    # label 20: table 10 (X = 1 - Z), meant to be skipped for Z >= 0.1.
    if not z >= 0.1:
        if zzg[idz + 3] >= 0.1:
            idz2 = NUMZ - 4
            for i in range(1, 5):
                a[i] = zzg[idz2 + i - 1]
            intrp2(a, bb, z)
        else:
            idz2 = idz
            for i in range(1, 5):
                bb[i] = b[i]
    # label 30
    # Question (see the page): label 30 stands before the loop, so table 10
    # is filled even for Z >= 0.1, with IDZ2 and BB from earlier.
    _fill_table(10, idz2, bb)
