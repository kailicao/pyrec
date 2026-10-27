"""Python version of alsurfp.f (YREC 09c3b9b): the outer boundary
condition from the Allard model atmospheres, by 4x4 Lagrange interpolation
in log g and log Teff."""

import numpy as np

from pyrec import alatm01_mod
from pyrec import alatm02_mod
from pyrec import atmprt_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.locate import locate
from pyrec.polint import polint

# SAVE: the work arrays keep their values between calls.
igl = np.zeros(5, dtype=int)
pl1 = np.zeros(5)
p100l1 = np.zeros(5)
t100l1 = np.zeros(5)
pl2 = np.zeros(5)
p100l2 = np.zeros(5)
t100l2 = np.zeros(5)


def _both(text):
    """The same list-directed line to ISHORT and to the screen."""
    iounits_mod.write(iounits_mod.ishort, text)
    print(text)


def alsurfp(teffl, gl, lprt, lalfail):
    """Interpolate the Allard table at (log Teff, log g): log P at
    T = Teff (AT = log Teff), or, with LALTPTau100, log P and log T at
    tau = 100, into AP and AT of ATMPRT_MOD. Sets LAlFail above the table
    in Teff or outside it in log g; stops the run below it in Teff.

    Translates: alsurfp.f#ALSURFP (09c3b9b)

    Arguments: teffl (log Teff), gl (log g), lprt (print to ISTOR too) and
    lalfail (set here). Returns lalfail.
    """
    a1 = alatm01_mod
    a2 = alatm02_mod
    at = atmprt_mod
    ltau100 = nmlpar_mod.laltptau100
    ishort = iounits_mod.ishort
    istor = iounits_mod.istor

    lbadpoint = False  # presume we have started at a good point
    if teffl >= a2.tefflmax:
        lbadpoint = True
        for show in (lambda t: iounits_mod.write(ishort, t), print):
            show("")
            show(" ALSURFP: TEFFL greater than TeffLmax: ")
            show(f"     TEFFLmin,TEFFL,TEFFLmax;  {a2.tefflmin} {teffl}"
                 f" {a2.tefflmax}")
    if gl >= a2.glxmax or gl <= a2.glxmin:
        lbadpoint = True
        iounits_mod.write(ishort, "")
        print("")
        for show in (lambda t: iounits_mod.write(ishort, t), print):
            show(" ALSURFP: GL out of max range: ")
            show(f"     GLXmin,GL,GLXmax;  {a2.glxmin} {gl} {a2.glxmax}")
    if lbadpoint:
        # Above the table in Teff, or outside it in log g: the caller
        # takes action.
        return True
    lalfail = False

    # Question (see the page): below the table in Teff the run stops.
    if teffl <= a2.tefflmin:
        for show in (lambda t: iounits_mod.write(ishort, t), print):
            show("")
            show(" ALSURFP: TEFFL less than TEFFLmin: ")
            show(f"     TEFFLmin,TEFFL,TEFFLmax;  {a2.tefflmin} {teffl}"
                 f" {a2.tefflmax}")
        _error_exit()

    # The first of the 4 rows in log Teff, kept inside the table.
    iteffl = locate(a1.teffls, a1.nteff, teffl, 0)
    iteffl = iteffl - 1
    if iteffl < 1:
        iteffl = 1
    if iteffl > a1.nteff - 3:
        iteffl = a1.nteff - 3

    # The first of the 4 columns in log g, in each row's valid range.
    for i in range(1, 5):
        i1 = iteffl + i - 1
        j1 = a2.iglmin[i1]
        j2 = a2.iglmax[i1]
        n = j2 - j1 + 1
        # GLs(j1) passed as an array: the section whose index 1 is GLs(j1).
        j3 = locate(a1.gls[j1 - 1:], n, gl, 0)
        igl[i] = j1 + j3 - 1
        igl[i] = igl[i] - 1
        if igl[i] < j1:
            igl[i] = j1
        if igl[i] > j2 - 3:
            igl[i] = j2 - 3
        if gl < a2.glmin[i1] or gl > a2.glmax[i1]:
            for show in (lambda t: iounits_mod.write(ishort, t), print):
                show("")
                show(" ALSURFP: GL out of extended range: ")
                show(f"     GLmin,GL,GLmax;  {a2.glmin[i1]} {gl}"
                     f" {a2.glmax[i1]}")
            return True

    # Interpolate each row in log g.
    for i in range(1, 5):
        i1 = iteffl + i - 1
        for j in range(1, 5):
            j1 = igl[i] + j - 1
            pl1[j] = a1.pls[i1, j1]
            p100l1[j] = a1.p100ls[i1, j1]
            t100l1[j] = a1.t100ls[i1, j1]
        j1 = igl[i]
        xa = a1.gls[j1 - 1:]
        pl2[i], dy = polint(xa, pl1, 4, gl, pl2[i], 0.0)
        p100l2[i], dy = polint(xa, p100l1, 4, gl, p100l2[i], dy)
        t100l2[i], dy = polint(xa, t100l1, 4, gl, t100l2[i], dy)

    # Then the 4 row values in log Teff.
    xa = a1.teffls[iteffl - 1:]
    if not ltau100:
        at.ap, dy = polint(xa, pl2, 4, teffl, at.ap, dy)
        at.at = teffl
    else:
        at.ap, dy = polint(xa, p100l2, 4, teffl, at.ap, dy)
        at.at, dy = polint(xa, t100l2, 4, teffl, at.at, dy)

    # Question (see the page): two lines to ISHORT at every call.
    if lprt:
        if ltau100:
            text = ("********PL,TL at Tau=100 INTERPOLATED"
                    " FROM ALLARD TABULATED VALUES********")  # 70
        else:
            text = ("********PL,TL at T=TEFF,GL INTERPOLATED "
                    "FROM ALLARD TABULATED VALUES********")  # 71
        iounits_mod.write(ishort, text)
        iounits_mod.write(istor, text)
        line = (f" {'':20s}LOG(Teffl) ={teffl:10.5f} ,LOG(G) ={gl:10.5f}"
                f", LOG(T) ={at.at:10.5f}, Log(P) ={at.ap:10.5f}")  # 72
        iounits_mod.write(ishort, line)
        iounits_mod.write(istor, line)
    if not lprt:
        if ltau100:
            text = ("ALSURFP:  PL,TL at Tau=100 INTERPOLATED"
                    " FROM ALLARD TABULATED VALUES:")  # 73
        else:
            text = ("ALSURFP: PL,TL at T=TEFF,GL INTERPOLATED "
                    "FROM ALLARD TABULATED VALUES:")  # 74
        iounits_mod.write(ishort, text)
        iounits_mod.write(ishort,
                          f" {'':10s}LOG(Teffl) ={teffl:10.5f}, LOG(G) ="
                          f"{gl:10.5f}, LOG(T) ={at.at:10.5f}, Log(P) ="
                          f"{at.ap:10.5f}")  # 75
    return lalfail


def _error_exit():
    """label 9999: the error exit (STOP)."""
    print("")
    print(" ******** ALSURFP: Program Terminated ********")
    print("")
    iounits_mod.write(iounits_mod.ishort, "")
    iounits_mod.write(iounits_mod.ishort,
                      " ******** ALSURF: Program Terminated ********")
    iounits_mod.write(iounits_mod.ishort, "")
    raise SystemExit
