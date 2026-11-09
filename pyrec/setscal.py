"""Python version of setscal.f (YREC 7be4133): the 50 runs of a
calibration of a star to a given luminosity and Teff (or radius)."""

import math

from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod


def setscal(alri: float, lstar: bool, lpassr: bool) -> tuple[float, bool, bool]:
    """Compute the target's missing radius or Teff from L = 4 pi R^2
    sigma Teff^4, reset the search's state, and fill the run list: 25
    repetitions of runs 1 (a new starting model) and 2 (evolution).

    Translates: setscal.f#SETSCAL (7be4133)

    Arguments: alri, log R/Rsun of the previous model; lstar, lpassr, the
    search's flags (all only set).

    Returns (alri, lstar, lpassr) = (0, False, False).

    """

    c = constants_mod
    n = nmlpar_mod
    r = runloop_mod

    # The search's flags (Fortran lines 100–101).
    lstar = False
    lpassr = False

    # The target's missing radius or Teff (Fortran lines 102–107).
    if n.lteff:
        n.sr = math.sqrt(n.xls * c.clsun / (c.c4pi * c.csig)) / (
            n.steff * n.steff * c.crsun)
    else:
        n.steff = ((n.xls * c.clsun)
                   / (c.c4pi * c.csig * n.sr * n.sr * c.crsun * c.crsun)
                   )**0.25
    alri = 0.0

    # 50 runs; run 1's composition and mixing length copied to all
    # (Fortran lines 108–119).
    # Question (see the page): the comment says 24 runs.
    n.numrun = 50
    for i in range(2, 51):
        n.xenv0a[i] = n.xenv0a[1]
        n.zenv0a[i] = n.zenv0a[1]
        n.cmixla[i] = n.cmixla[1]
        n.lsenv0a[i] = n.lsenv0a[1]
        n.senv0a[i] = n.senv0a[1]

    # Runs 3, 5, ..., 49 repeat run 1 (a new starting model) (Fortran lines 120–130).
    for i in range(3, 50, 2):
        r.iresca[i] = r.iresca[1]
        n.lfirst[i] = True
        n.nmodls[i] = n.nmodls[1]
        n.rsclzc[i] = n.rsclzc[1]
        n.rsclzm1[i] = n.rsclzm1[1]
        n.rsclzm2[i] = n.rsclzm2[1]
        for j in range(1, 5):
            r.rescal[j, i] = r.rescal[j, 1]

    # Runs 4, 6, ..., 50 repeat run 2 (evolution) (Fortran lines 131–139).
    # Question (see the page): END_DCEN, END_XCEN, END_YCEN not copied.
    for i in range(4, 51, 2):
        r.iresca[i] = 1
        n.lfirst[i] = False
        n.nmodls[i] = n.nmodls[2]
        n.endage[i] = n.endage[2]
        r.lendag[i] = r.lendag[2]
        n.setdt[i] = n.setdt[2]
        r.lsetdt[i] = r.lsetdt[2]

    # Fortran lines 140–141.
    print(" Evolve to R*, L* = ", n.sr, n.xls)
    iounits_mod.write(iounits_mod.itrack,
                      f" #Evolve to R*, L* = {n.sr} {n.xls}")
    return alri, lstar, lpassr
