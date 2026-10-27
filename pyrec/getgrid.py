"""Python version of getgrid.f (YREC 09c3b9b): the equally spaced grid
for the rotational mixing and angular momentum diffusion."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import egrid_mod
from pyrec import nmlpar_mod
from pyrec.osplin import osplin

# SAVE: the work array of log m on the region's model points.
ytab = np.zeros(dimensions_mod.json + 1)


def getgrid(hl, hp, hs, ibeg, iend, m):
    """Build the equally spaced grid in
    chi = log m / D_m + L / D_L - log P / D_P over the region ibeg..iend,
    with as many points as the region has, and the mass at each grid
    point by OSPLIN's spline in chi.

    Translates: getgrid.f#GETGRID (09c3b9b)

    Arguments: hl (L, solar units), hp (log P), hs (log m), ibeg, iend
    (the region), m (number of shells).
    Changed in place: egrid_mod.chi, echi, es1; egrid_mod.dchi and
    egrid_mod.ntot are set.

    Returns nothing.
    """
    e = egrid_mod
    hpttol = nmlpar_mod.hpttol

    dm = hpttol[2]
    dl = hpttol[9] * hl[m]
    dp = hpttol[11]
    ntab = iend - ibeg + 1
    for ii in range(1, ntab + 1):
        i = ibeg + ii - 1
        e.chi[ii] = hs[i] / dm + hl[i] / dl - hp[i] / dp
    # Question (see the page): as many grid points as model points.
    e.ntot = ntab
    e.ntot = max(e.ntot, 3)
    e.dchi = (e.chi[ntab] - e.chi[1]) / float(e.ntot - 1)
    e.echi[1] = e.chi[1]
    for i in range(2, e.ntot + 1):
        e.echi[i] = e.echi[i - 1] + e.dchi
    # log m at the grid points, then unlogged.
    for i in range(1, ntab + 1):
        ytab[i] = hs[ibeg + i - 1]
    osplin(e.echi, e.es1, e.chi, ytab, ntab, e.ntot)
    for i in range(1, e.ntot + 1):
        e.es1[i] = math.exp(constants_mod.cln * e.es1[i])
