"""Python version of getalex06.f (YREC 7be4133): the Ferguson-Alexander
2006 low-temperature opacity."""

import math

import numpy as np

from pyrec import alot06_mod
from pyrec import constants_mod
from pyrec import galot06_mod
from pyrec.alex06tab import alex06tab
from pyrec.interp import interp

numt = 85
numd = 19


def getalex06(dl: float, tl: float, x: float, z: float, o: float, ol: float, qod: float, qot: float) -> tuple[float,
              float, float, float]:
    """Interpolate the Ferguson-Alexander 2006 opacity in log T and
    log R = log rho - 3 (log T - 6), with four-point Lagrange weights,
    from the table for the current X, Z (rebuilt by ALEX06TAB when they
    change); extrapolate linearly above the largest log R, and convert
    the temperature derivative to fixed rho.

    Translates: getalex06.f#GETALEX06 (7be4133)

    Arguments: dl, tl (log rho, log T), x, z; o, ol, qod, qot are only
    set.

    Returns (o, ol, qod, qot): kappa, log kappa, dlog kappa/dlog rho at
    fixed T, dlog kappa/dlog T at fixed rho.

    """

    a = alot06_mod
    g = galot06_mod
    qr = np.zeros(4 + 1)
    fxt = np.zeros(4 + 1)
    fxtd = np.zeros(4 + 1)
    fxr = np.zeros(4 + 1)
    fxrd = np.zeros(4 + 1)
    op = np.zeros(4 + 1)
    qopr = np.zeros(4 + 1)

    # Ensure that the table has the point's X, Z (Fortran lines 26–33).
    # Question (see the page): the table is rebuilt at every change of
    # composition.
    difz = abs(z - a.ze)
    difx = abs(x - a.xe)
    if difz > 1.0e-8 or difx > 1.0e-8:
        a.ze = z
        a.xe = x
        alex06tab()

    # R for the grid point (Fortran line 35).
    r = dl - 3.0 * (tl - 6.0)

    # The 4-point stencil in log T, searched from the last one (Fortran
    # lines 36–59).
    # Question (see the page): no check below or above the table.
    if a.it < 1:
        a.it = 1
    if a.it + 2 > numt:
        a.it = numt - 2
    if tl < g.tgr06[a.it + 2]:
        for i in range(a.it + 1, 1, -1):
            if tl > g.tgr06[i]:
                a.it = i - 1
                break  # GOTO 10
        else:
            a.it = 1
    else:
        for i in range(a.it + 3, numt + 1):
            if tl < g.tgr06[i]:
                a.it = i - 2
                a.it = min(numt - 3, a.it)
                break  # GOTO 20
        else:
            a.it = numt - 3

    # Interpolation factors in log T (Fortran lines 60–64).
    for i in range(1, 4 + 1):
        qr[i] = g.tgr06[a.it + i - 1]
    interp(qr, fxt, fxtd, tl)

    # The same in log R (Fortran lines 65–88).
    if a.id < 1:
        a.id = 1
    if a.id + 2 > numd:
        a.id = numd - 2
    if r < g.r0gr06[a.id + 2]:
        for i in range(a.id + 1, 1, -1):
            if r > g.r0gr06[i]:
                a.id = i - 1
                break  # GOTO 30
        else:
            a.id = 1
    else:
        for i in range(a.id + 3, numd + 1):
            if r < g.r0gr06[i]:
                a.id = i - 2
                a.id = min(numd - 3, a.id)
                break  # GOTO 40
        else:
            a.id = numd - 3

    # Above the table's largest log R: interpolate at the edge, then
    # extrapolate linearly (Fortran lines 89–100).
    if r > g.r0gr06[numd]:
        linear = True
        rr = r
        r = g.r0gr06[numd]
    else:
        linear = False
    for i in range(1, 4 + 1):
        qr[i] = g.r0gr06[a.id + i - 1]
    interp(qr, fxr, fxrd, r)

    # Interpolate in log R at fixed T (Fortran lines 101–112).
    cap = a.capalex06
    idd = a.id
    for i in range(1, 4 + 1):
        ii = a.it + i - 1
        op[i] = (fxr[1] * cap[ii, idd] + fxr[2] * cap[ii, idd + 1]
                 + fxr[3] * cap[ii, idd + 2] + fxr[4] * cap[ii, idd + 3])
        qopr[i] = (fxrd[1] * cap[ii, idd] + fxrd[2] * cap[ii, idd + 1]
                   + fxrd[3] * cap[ii, idd + 2] + fxrd[4] * cap[ii, idd + 3])

    # The linear extrapolation in log R (Fortran lines 113–120).
    if linear:
        # Suspected bug (see the page): qopr, and so qod, is not
        # extrapolated.
        for i in range(1, 4 + 1):
            ii = a.it + i - 1
            op[i] = op[i] + (rr - r) * (cap[ii, numd] - cap[ii, numd - 1]) / (
                g.r0gr06[numd] - g.r0gr06[numd - 1])
        r = rr

    # Interpolate in T (Fortran lines 121–122).
    ol = fxt[1] * op[1] + fxt[2] * op[2] + fxt[3] * op[3] + fxt[4] * op[4]

    # dlog kappa/dlog T at fixed R, then at fixed rho (Fortran lines
    # 123–130).
    qot = (fxtd[1] * op[1] + fxtd[2] * op[2] + fxtd[3] * op[3]
           + fxtd[4] * op[4])
    qod = (fxt[1] * qopr[1] + fxt[2] * qopr[2] + fxt[3] * qopr[3]
           + fxt[4] * qopr[4])
    qot = qot - 3.0 * qod
    o = math.exp(constants_mod.cln * ol)
    return o, ol, qod, qot
