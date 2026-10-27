"""Python version of setcal.f (YREC 09c3b9b): the 48 runs of an automatic
solar calibration, from the first three."""

from pyrec import nmlpar_mod
from pyrec import runloop_mod


def setcal(facage, lsol):
    """Fill the run list for a solar calibration: 16 repetitions of a
    three-run pattern (a new starting model, a run to 1e8 yr, a run to
    the solar age), copied from runs 1-3, with ages and steps scaled by
    facage.

    Translates: setcal.f#SETCAL (09c3b9b)

    Arguments: facage, the age factor (1 except in Monte Carlo runs);
    lsol, the calibration's convergence flag (only set).

    Returns lsol, set to False.
    """
    n = nmlpar_mod
    r = runloop_mod
    # Question (see the page): the comment says two runs and 24.
    # mhp 5/96: solar models in 3 runs rather than 2.
    n.numrun = 48
    for i in range(2, 49):
        n.xenv0a[i] = n.xenv0a[1]
        n.zenv0a[i] = n.zenv0a[1]
        n.cmixla[i] = n.cmixla[1]
        n.lsenv0a[i] = n.lsenv0a[1]
        n.senv0a[i] = n.senv0a[1]
    # Runs 4, 7, ..., 46 repeat run 1 (a new starting model).
    for i in range(4, 47, 3):
        r.iresca[i] = r.iresca[1]
        n.lfirst[i] = True
        n.nmodls[i] = n.nmodls[1]
        n.rsclzc[i] = n.rsclzc[1]
        n.rsclzm1[i] = n.rsclzm1[1]
        n.rsclzm2[i] = n.rsclzm2[1]
        for j in range(1, 5):
            r.rescal[j, i] = r.rescal[j, 1]
    # MHP 06/13: run 2 hardwired to 1e8 yr, run 3 to CALSOLAGE.
    # Question (see the page): the hard-wired 1e8 yr.
    n.endage[2] = 1.0e8 * facage
    n.setdt[2] = n.setdt[2] * facage
    n.endage[3] = n.calsolage * facage
    n.setdt[3] = n.setdt[3] * facage
    # Question (see the page): END_DCEN, END_XCEN, END_YCEN not copied.
    for i in range(5, 48, 3):
        r.iresca[i] = 1
        n.lfirst[i] = False
        n.nmodls[i] = n.nmodls[2]
        n.endage[i] = n.endage[2]
        r.lendag[i] = r.lendag[2]
        n.setdt[i] = n.setdt[2]
        r.lsetdt[i] = r.lsetdt[2]
    for i in range(6, 49, 3):
        r.iresca[i] = 1
        n.lfirst[i] = False
        n.nmodls[i] = n.nmodls[3]
        n.endage[i] = n.endage[3]
        r.lendag[i] = r.lendag[3]
        n.setdt[i] = n.setdt[3]
        r.lsetdt[i] = r.lsetdt[3]
    lsol = False
    return lsol
