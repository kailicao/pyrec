"""Python version of alprint.f (YREC 09c3b9b): a diagnostic print of the
Allard atmosphere table, called by ALTABINIT before it stops the run."""

from pyrec import alatm01_mod
from pyrec import alatm02_mod
from pyrec import iounits_mod


def _row(i, values):
    """FORMAT 910, (I5,F6.0,2I4,2F6.2,1P20E16.8): one row per Teff."""
    a1 = alatm01_mod
    a2 = alatm02_mod
    text = (f"{i:5d}{10.0 ** a1.teffls[i]:6.0f}{a2.iglmin[i]:4d}"
            f"{a2.iglmax[i]:4d}{a2.glmin[i]:6.2f}{a2.glmax[i]:6.2f}")
    text += ''.join(f"{values[i, j]:16.8E}" for j in range(1, a1.ngl + 1))
    return text


def alprint():
    """Print the Allard table to ISHORT: its sizes, its limits, and, one row
    per Teff with its valid range in log g, the tables of log P at Teff,
    log P at tau = 100 and log T at tau = 100. Also recomputes TEFFLmin and
    TEFFLmax in ALATM02_MOD, with ALTABINIT's formulas.

    Translates: alprint.f#ALPRINT (09c3b9b)

    No arguments; nothing is returned. List-directed WRITEs are rendered
    with Python's default formatting.
    """
    a1 = alatm01_mod
    a2 = alatm02_mod
    unit = iounits_mod.ishort
    w = iounits_mod.write

    # Question (see the page): a side effect on ALATM02_MOD.
    a2.tefflmin = a1.teffls[1] - (a1.teffls[2] - a1.teffls[1])
    a2.tefflmax = (a1.teffls[a1.nteff]
                   + (a1.teffls[a1.nteff] - a1.teffls[a1.nteff - 1]))

    w(unit, "")
    w(unit, " **************** Allard Table  ***************")
    w(unit, "")
    w(unit, f" nTeff= {a1.nteff} ,nGL= {a1.ngl}")
    w(unit, f" nFeH= {a1.nfeh} ,nAlpha= {a1.nalpha}")
    w(unit, "")
    w(unit, f" TEFFLmin= {a2.tefflmin} ,TEFFLmax= {a2.tefflmax}")
    w(unit, f" GLXmin= {a2.glxmin} GLXmax= {a2.glxmax}")
    w(unit, "")
    # FORMAT 900, (A10,16X,A,20(5X,F6.2,5X)).
    w(unit, f"{' PL @ Teff':10s}{'':16s}GL ->"
            + ''.join(f"{'':5s}{a1.gls[j]:6.2f}{'':5s}"
                      for j in range(1, a1.ngl + 1)))
    for title, table in (("PL at Teff,GL", a1.pls),
                         ("PL at Tau=100", a1.p100ls),
                         ("TL at Tau=100", a1.t100ls)):
        w(unit, "")
        w(unit, f" {title}")
        w(unit, "")
        for i in range(1, a1.nteff + 1):
            w(unit, _row(i, table))
    w(unit, "")
