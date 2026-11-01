"""Python version of getopal95.f (YREC development, 09c3b9b): the driver for
the OPAL 1995 opacity tables."""

import numpy as np

from pyrec import acdpth_mod
from pyrec import llot95_mod
from pyrec import llot95a_mod
from pyrec import llot95e_mod
from pyrec import op95interp_mod
from pyrec.interp import interp
from pyrec.intrp2 import intrp2
from pyrec.op952d import op952d
from pyrec.op953d import op953d
from pyrec.op954d import op954d

NUMT = 70
NUMD = 19
NUMX = 10
NUMZ = 13

# GETOPAL95 has SAVE: its locals keep their values between calls. Only
# LSHIFT is read before it is set (Question 2); the work arrays A, B, C are
# kept here as well, as in the Fortran.
a = np.zeros(5)
b = np.zeros(5)
c = np.zeros(5)
lshift = False
lchk = False


def _x_stencil(k, x, zz):
    """The X stencil JX(K,1..4) and its weights BX(K,1..4) at the Z value zz.

    The same block appears twice in the Fortran: for the 3D path (K = 1,
    with ZZ = ZTAB, labels 70 and 80) and for each Z node of the 4D path
    (labels 75 and 85). It is one helper here; the statements are those of
    the Fortran.
    """
    jx = op95interp_mod.jx
    xxg = llot95a_mod.xxg
    if x > xxg[jx[k, 1] + 2]:
        for i in range(jx[k, 1] + 3, NUMX):         # I = JX(K,1)+3, NUMX-1
            if x < xxg[i]:
                jx[k, 1] = i - 2
                break
        else:
            jx[k, 1] = NUMX - 3
    else:
        for i in range(jx[k, 1] + 1, 1, -1):        # I = JX(K,1)+1, 2, -1
            if x >= xxg[i]:
                jx[k, 1] = i - 1
                break
        else:
            jx[k, 1] = 1
    # Check for absent tables at high X.
    if zz >= 0.04 and x >= 0.8:
        if zz < 0.1:
            # No table 9 (X = 0.95).
            jx[k, 4] = NUMX
            jx[k, 3] = NUMX - 2
            jx[k, 2] = NUMX - 3
            jx[k, 1] = NUMX - 4
        else:
            # No table 9 or table 10.
            jx[k, 4] = NUMX - 2
            jx[k, 3] = NUMX - 3
            jx[k, 2] = NUMX - 4
            jx[k, 1] = NUMX - 5
    else:
        jxx = jx[k, 1] - 1
        for i in range(1, 5):
            jx[k, i] = jxx + i
    # X interpolation factors.
    for i in range(1, 5):
        a[i] = xxg[jx[k, i]]
    if jx[k, 4] == NUMX:
        a[4] = 1.0 - zz
    intrp2(a, b, x)
    for j in range(1, 5):
        op95interp_mod.bx[k, j] = b[j]


def getopal95(dl, tl, x, z, o, ol, qod, qot):
    """Opacity from the OPAL 1995 tables, given rho, T, X and Z.

    Finds the interpolation stencils and weights in T, R = rho/T6**3, X and
    Z, handles the tables' empty high-density corner, and calls the 2D, 3D
    or 4D interpolator.

    Translates: getopal95.f#GETOPAL95 (09c3b9b)

    Parameters
    ----------
    dl, tl : float
        log10 density and log10 temperature.
    x, z : float
        Hydrogen and metal mass fractions.
    o, ol, qod, qot : float
        Set here (passed in with any value).

    Returns
    -------
    o, ol, qod, qot
        The opacity, its log10, and dlog(kappa)/dlog(rho) and
        dlog(kappa)/dlog(T), as from the interpolators.
    """
    global lshift, lchk
    im = op95interp_mod
    tgr = llot95a_mod.tgr
    r0gr = llot95a_mod.r0gr
    zzg = llot95a_mod.zzg
    cappa = llot95_mod.cappa
    cappa2 = llot95a_mod.cappa2
    jd = im.jd
    jx = im.jx

    # log R = log rho - 3 (log T - 6).
    im.rho_t63l = dl - 3.0 * (tl - 6.0)
    if tl < tgr[1] or tl > tgr[NUMT]:
        print(f" LOG T OF{tl:11.6f}OUT OF OPAL 95 TABLE RANGE")   # format 5
        raise SystemExit
    # Extrapolation below the first density of the table?
    im.lextl = im.rho_t63l < r0gr[1]

    # Stencil in T, starting from the last one.
    if tl >= tgr[im.jt + 2]:
        for i in range(im.jt + 3, NUMT):
            if tl < tgr[i]:
                im.jt = i - 2
                break                                   # label 10
        else:
            im.jt = NUMT - 3
    else:
        for i in range(im.jt + 1, 1, -1):
            if tl > tgr[i]:
                im.jt = i - 1
                break                                   # label 20
        else:
            im.jt = 1
    jt = im.jt
    for i in range(1, 5):
        a[i] = tgr[jt + i - 1]
    interp(a, b, c, tl)
    for i in range(1, 5):
        im.bt[i] = b[i]
        im.btt[i] = c[i]

    # Stencil in R (for the first temperature).
    if im.rho_t63l >= r0gr[jd[1] + 2]:
        for i in range(jd[1] + 3, NUMD):
            if im.rho_t63l < r0gr[i]:
                jd[1] = i - 2
                break                                   # label 30
        else:
            jd[1] = NUMD - 3
    else:
        for i in range(jd[1] + 1, 1, -1):
            if im.rho_t63l > r0gr[i]:
                jd[1] = i - 1
                break                                   # label 40
        else:
            jd[1] = 1

    # The empty corner. At low X and low T the check is postponed (and
    # never made): GOTO 60 skips the reset of LSHIFT and LEXTH.
    if x < 0.2 and jt <= 5:
        # Question (see the page): LSHIFT and LEXTH keep their values from
        # the last call on this path; LCHK is never used.
        lchk = True
    else:
        lchk = False
        lshift = False
        im.lexth = False
        # Question (see the page): JX(1,1) is the X stencil of the last call.
        if cappa[jx[1, 1], jt + 3, jd[1] + 3] >= 9.9:
            lshift = True
            for k in range(1, 5):
                if cappa[jx[1, 1], jt + k - 1, jd[1] + 1] >= 9.9:
                    im.lexth = True
                    im.lexh[k] = True
                else:
                    im.lexh[k] = False
            # The four nearest densities with data at each temperature.
            ixx = jx[1, 1]
            idd = jd[1]
            for i in range(1, 5):
                itt = jt + i - 1
                if cappa[ixx, itt, idd + 3] > 9.9:
                    for j in range(idd + 2, 0, -1):
                        # Question (see the page): CAPPA2 here, CAPPA above.
                        if cappa2[ixx, itt, j] <= 9.9:
                            # Question (see the page): can be zero or less.
                            jd[i] = j - 3
                            idd = j - 3
                            break                       # label 50
                else:
                    jd[i] = idd
    # label 60

    # Weights in R for each of the four temperatures.
    if not lshift:
        for i in range(1, 5):
            a[i] = r0gr[jd[1] + i - 1]
        if not im.lextl:
            interp(a, b, c, im.rho_t63l)
        else:
            im.rho_t63ll = r0gr[1]
            interp(a, b, c, im.rho_t63ll)
        for j in range(1, 5):
            jd[j] = jd[1]
            for i in range(1, 5):
                im.bd[j, i] = b[i]
                im.bdd[j, i] = c[i]
    else:
        for j in range(1, 5):
            for i in range(1, 5):
                a[i] = r0gr[jd[j] + i - 1]
            # Outside the data: weights at the edge, linear extrapolation
            # later.
            if not im.lexh[j]:
                interp(a, b, c, im.rho_t63l)
            else:
                rr = r0gr[jd[j] + 3]
                interp(a, b, c, rr)
                im.rlh[j] = rr
            for i in range(1, 5):
                im.bd[j, i] = b[i]
                im.bdd[j, i] = c[i]

    # 2D, 3D or 4D? The 4D path is forced for acoustic-depth output (JVS 04/11).
    if not (acdpth_mod.ladon and acdpth_mod.lacout):
        ztab = llot95_mod.ztab
        if abs(z - ztab) / max(ztab, 1.0e-6) <= 1.0e-4:
            if abs(x - llot95e_mod.xop95) <= 1.0e-4:
                # 2D in the surface table.
                o, ol, qod, qot = op952d(o, ol, qod, qot)
                return o, ol, qod, qot
            # 3D in the table at fixed Z (X, T, rho).
            _x_stencil(1, x, ztab)
            o, ol, qod, qot = op953d(o, ol, qod, qot)
            return o, ol, qod, qot

    # label 153: 4D in Z, X, T, rho. Nearest tables in Z.
    if z > zzg[im.jz + 2]:
        for i in range(im.jz + 3, NUMZ):  # up to NUMZ, a typo fixed (mhp 7/12)
            if z < zzg[i]:
                im.jz = i - 2
                break                                   # label 110
        else:
            im.jz = NUMZ - 3
    else:
        for i in range(im.jz + 1, 1, -1):
            if z >= zzg[i]:
                im.jz = i - 1
                break                                   # label 120
        else:
            im.jz = 1
    im.jz = min(im.jz, NUMZ - 3)  # failsafe: four Z points (mhp 7/12)
    for i in range(1, 5):
        a[i] = zzg[im.jz + i - 1]
    intrp2(a, b, z)
    for i in range(1, 5):
        im.bz[i] = b[i]
    # Nearest tables in X at each value of Z.
    for k in range(1, 5):
        _x_stencil(k, x, zzg[im.jz + k - 1])
    o, ol, qod, qot = op954d(o, ol, qod, qot)
    return o, ol, qod, qot
