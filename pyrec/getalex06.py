"""Python version of getalex06.f (YREC 09c3b9b): the Ferguson-Alexander
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


def getalex06(dl, tl, x, z, o, ol, qod, qot):
    """Interpolate the Ferguson-Alexander 2006 opacity in log T and
    log R = log rho - 3 (log T - 6), with four-point Lagrange weights,
    from the table for the current X, Z (rebuilt by ALEX06TAB when they
    change); extrapolate linearly above the largest log R, and convert
    the temperature derivative to fixed rho.

    Translates: getalex06.f#GETALEX06 (09c3b9b)

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

    # Question (see the page): the table is rebuilt at every change of
    # composition.
    difz = abs(z - a.ze)
    difx = abs(x - a.xe)
    if difz > 1.0e-8 or difx > 1.0e-8:
        a.ze = z
        a.xe = x
        alex06tab()

    r = dl - 3.0 * (tl - 6.0)

    # The 4-point stencil in log T, searched from the last one.
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
    for i in range(1, 4 + 1):
        qr[i] = g.tgr06[a.it + i - 1]
    interp(qr, fxt, fxtd, tl)

    # The same in log R.
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
    # extrapolate linearly.
    if r > g.r0gr06[numd]:
        linear = True
        rr = r
        r = g.r0gr06[numd]
    else:
        linear = False
    for i in range(1, 4 + 1):
        qr[i] = g.r0gr06[a.id + i - 1]
    interp(qr, fxr, fxrd, r)

    cap = a.capalex06
    idd = a.id
    for i in range(1, 4 + 1):
        ii = a.it + i - 1
        op[i] = (fxr[1] * cap[ii, idd] + fxr[2] * cap[ii, idd + 1]
                 + fxr[3] * cap[ii, idd + 2] + fxr[4] * cap[ii, idd + 3])
        qopr[i] = (fxrd[1] * cap[ii, idd] + fxrd[2] * cap[ii, idd + 1]
                   + fxrd[3] * cap[ii, idd + 2] + fxrd[4] * cap[ii, idd + 3])
    if linear:
        # Question (see the page): qod is not extrapolated.
        for i in range(1, 4 + 1):
            ii = a.it + i - 1
            op[i] = op[i] + (rr - r) * (cap[ii, numd] - cap[ii, numd - 1]) / (
                g.r0gr06[numd] - g.r0gr06[numd - 1])
        r = rr

    ol = fxt[1] * op[1] + fxt[2] * op[2] + fxt[3] * op[3] + fxt[4] * op[4]
    # dlog kappa/dlog T at fixed R, then at fixed rho.
    qot = (fxtd[1] * op[1] + fxtd[2] * op[2] + fxtd[3] * op[3]
           + fxtd[4] * op[4])
    qod = (fxt[1] * qopr[1] + fxt[2] * qopr[2] + fxt[3] * qopr[3]
           + fxt[4] * qopr[4])
    qot = qot - 3.0 * qod
    o = math.exp(constants_mod.cln * ol)
    return o, ol, qod, qot
