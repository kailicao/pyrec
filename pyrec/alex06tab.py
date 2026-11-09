"""Python version of alex06tab.f (YREC 7be4133): the Alexander (2006)
low-temperature opacity table interpolated to one composition (X, Z)."""

import numpy as np

from pyrec import alot06_mod
from pyrec import alot06all_mod
from pyrec import galot06_mod
from pyrec import iounits_mod
from pyrec.intrp2 import intrp2

NUMX = 9
NUMZ = 16
NUMT = 85
NUMD = 19

# SAVE: the work arrays keep their values between calls.
qr = np.zeros(5)
fx = np.zeros(5)
captemp = np.zeros((5, NUMT + 1, NUMD + 1))


def _stop_composition():
    """FORMAT 5 to the terminal and the .short file, then STOP 1."""

    a = alot06_mod
    text = (f"ILLEGAL COMPOSITION (X,Z) = {a.xe:6.2f}{a.ze:6.2f}"
            f" IN ALEX06.RUN STOPPED")
    print(text)
    iounits_mod.write(iounits_mod.ishort, text)
    raise SystemExit(1)    # STOP 1, the run-time error code


def _combine(kk, kk2, kk3, kk4):
    """FX(1)*T(KK) + FX(2)*T(KK2) + FX(3)*T(KK3) + FX(4)*T(KK4), over all
    (T, R), in the Fortran's order of the sum."""

    t = alot06all_mod.captalex06
    return (fx[1] * t[kk, 1:, 1:] + fx[2] * t[kk2, 1:, 1:]
            + fx[3] * t[kk3, 1:, 1:] + fx[4] * t[kk4, 1:, 1:])


def alex06tab() -> None:
    """Interpolate the stored Alexander tables CAPTALEX06 to the composition
    XE, ZE of ALOT06_MOD, first in Z, then in X, each with 4-point Lagrange
    weights, giving CAPALEX06 (and setting IX).

    Translates: alex06tab.f#ALEX06TAB (7be4133)

    No arguments; nothing is returned. The tables are stored as
    CAPTALEX06(JJ, T, R) with JJ = NUMZ*(IX-1) + IZ.

    """

    a = alot06_mod
    g = galot06_mod

    # XE = desired X; ZE = desired Z. Check that the requested composition
    # is inside the table bounds (Fortran lines 28–37).
    xmax = 1.0 - a.ze
    if a.xe < 0.0 or a.xe > xmax:
        _stop_composition()

    # Permit extrapolation in Z by up to 1 table element (Fortran lines
    # 38–45).
    # Question (see the page): is the range 0.10 < Z <= 0.12 meant for runs?
    zmax = g.zzg06[NUMZ] + (g.zzg06[NUMZ] - g.zzg06[NUMZ - 1])
    if a.ze < 0.0 or a.ze > zmax:
        _stop_composition()

    # Find the 4 nearest tables in Z (Fortran lines 46–54).
    for i in range(3, NUMZ - 2 + 1):
        if a.ze <= g.zzg06[i]:
            iz = i - 2
            break
    else:
        iz = NUMZ - 3

    # Label 7: find the 4 nearest tables in X (Fortran lines 55–68).
    for i in range(3, NUMX - 2 + 1):
        if a.xe <= g.xxg06[i]:
            a.ix = i - 2
            break
    else:
        # No table for X > 0.9 if Z = 0.10 or more.
        if a.ze >= 0.1:
            a.ix = NUMX - 4
        else:
            a.ix = NUMX - 3
    ix = a.ix

    # Label 10: interpolation factors for Z (Fortran lines 69–75).
    for i in range(1, 5):
        qr[i] = g.zzg06[iz + i - 1]
    intrp2(qr, fx, a.ze)

    # The first 3 X cases can be treated normally (Fortran lines 76–90).
    for k in range(1, 4):
        kk = NUMZ * (ix + k - 2) + iz
        captemp[k, 1:, 1:] = _combine(kk, kk + 1, kk + 2, kk + 3)

    # In the high-Z, high-X domain the top table is X = 1-Z (entry NUMX)
    # except for Z = 0.10, where the X = 0.9 table doubles for it
    # (Fortran lines 91–111).
    kk = NUMZ * (ix + 2) + iz
    if ix == NUMX - 3 and iz == NUMZ - 3:
        kk4 = NUMZ * (ix + 2)
    else:
        kk4 = kk + 3
    captemp[4, 1:, 1:] = _combine(kk, kk + 1, kk + 2, kk4)

    # Interpolation factors for X (Fortran lines 112–123).
    for i in range(1, 4):
        qr[i] = g.xxg06[ix + i - 1]
    if iz == NUMZ - 3:
        # Suspected bug (see the page): 1-Z is used whatever IX is.
        qr[4] = 1.0 - a.ze
    else:
        # Suspected bug (see the page): XXG06(9) = 1.0 stands for X = 1-Z.
        qr[4] = g.xxg06[ix + 3]
    intrp2(qr, fx, a.xe)

    # The interpolation in X (Fortran lines 124–131).
    a.capalex06[1:, 1:] = (fx[1] * captemp[1, 1:, 1:] + fx[2] * captemp[2, 1:, 1:]
                           + fx[3] * captemp[3, 1:, 1:] + fx[4] * captemp[4, 1:, 1:])
