"""Python version of altabinit.f (YREC 09c3b9b): checks the Allard
atmosphere table just read by ALFILEIN and derives the limits ALSURFP
needs."""

from pyrec import alatm01_mod
from pyrec import alatm02_mod
from pyrec import alatm05_mod
from pyrec import iounits_mod
from pyrec.alprint import alprint


def altabinit():
    """Set the table's range in log Teff (ALATM05_MOD: its own ends;
    ALATM02_MOD: one row's width beyond them), each row's first and last
    valid column in log g and its allowed range of log g, and the extremes
    GLXmin, GLXmax; then check the table, and stop the run (after ALPRINT)
    if it is bad.

    Translates: altabinit.f#ALTABINIT (09c3b9b)

    No arguments; nothing is returned.
    """
    a1 = alatm01_mod
    a2 = alatm02_mod
    a5 = alatm05_mod
    ishort = iounits_mod.ishort
    w = iounits_mod.write
    nteff = a1.nteff
    ngl = a1.ngl

    a5.al_tefflmin = a1.teffls[1]
    a5.al_tefflmax = a1.teffls[nteff]
    a2.tefflmin = a1.teffls[1] - (a1.teffls[2] - a1.teffls[1])
    a2.tefflmax = a1.teffls[nteff] + (a1.teffls[nteff] - a1.teffls[nteff - 1])

    # The first and last valid entries (above -998) of each row. Possible
    # bug (update-findings.md, batch 14): a row with none keeps its earlier
    # values (0 at first), and the loop below then reads GLs out of bounds
    # (here gls[0] and, through index -1, the array's last element).
    for i in range(1, nteff + 1):
        for j in range(1, ngl + 1):
            if a1.pls[i, j] > -998.0:
                a2.iglmin[i] = j
                break
        # label 100
        for j in range(ngl, 0, -1):
            if a1.pls[i, j] > -998.0:
                a2.iglmax[i] = j
                break
        # label 110

    # The allowed range of log g in each row, and over all rows.
    a2.glxmin = 999.0
    a2.glxmax = -999.0
    for i in range(1, nteff + 1):
        j1 = a2.iglmin[i]
        j2 = a2.iglmax[i]
        # Question (see the page): four column widths below, one above.
        a2.glmin[i] = a1.gls[j1] - 4.0 * (a1.gls[j1 + 1] - a1.gls[j1])
        if a2.glmin[i] < a2.glxmin:
            a2.glxmin = a2.glmin[i]
        a2.glmax[i] = a1.gls[j2] + (a1.gls[j2] - a1.gls[j2 - 1])
        if a2.glmax[i] > a2.glxmax:
            a2.glxmax = a2.glmax[i]

    # Question (see the page): only PLs is checked, not P100Ls or T100Ls.
    lbadtable = False
    for i in range(1, nteff + 1):
        j1 = a2.iglmin[i]
        j2 = a2.iglmax[i]
        for j in range(j1, j2 + 1):
            if a1.pls[i, j] < -998.0:
                lbadtable = True
                w(ishort, "ALTABINIT: Bad input Allard Table: TEFF, GL: "
                          f"{10.0 ** a1.teffls[i]:5.0f}{a1.gls[j]:7.2f}")  # 900

    if nteff < 4:
        lbadtable = True
        # Question (see the page): FORMAT 910, (A,I3), has one A for two
        # strings; the Fortran would stop here with a run-time format error.
        w(ishort, "ALTABINIT: Bad input Allard Table: ")
        raise SystemExit("Fortran run-time error: FORMAT 910 expects an "
                         "integer, got a string")
    for i in range(1, nteff + 1):
        j1 = a2.iglmin[i]
        j2 = a2.iglmax[i]
        if (j2 - j1 + 1) < 4:
            lbadtable = True
            w(ishort, "ALTABINIT: Bad input Allard Table: "
                      "Row with less that 4 elements: i,#,Teff,GLMin,GLmax: "
                      f"{i:4d}{j2 - j1 + 1:4d}  "
                      f"{10.0 ** a1.teffls[i]:7.2f}{a2.glmin[i]:7.2f}"
                      f"{a2.glmax[i]:7.2f}")  # 920

    if not lbadtable:
        return

    # label 9999: error exit.
    for unit_write in (print, lambda text: w(ishort, text)):
        unit_write("")
        unit_write(" ******** ALTABINIT: Program Terminated ********")
        unit_write("")
    alprint()
    raise SystemExit
