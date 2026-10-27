"""Python version of bursmix.f (YREC 09c3b9b): burning and mixing over one
time step in progressively smaller substeps, extrapolated to zero substep."""

import numpy as np

from pyrec import dimensions_mod
from pyrec import liov_mod
from pyrec import newrat_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldrat_mod
from pyrec.bsrotmix import bsrotmix
from pyrec.liburn2 import liburn2
from pyrec.ndifcom import ndifcom
from pyrec.rotmix import rotmix

JSON = dimensions_mod.json

# Locals kept between calls (SAVE). ECOD2, EM, IBEG and IEND are filled
# only by NDIFCOM (through MIXGRID) when LOK is true.
# Question (see the page): with LOK false they are read as left here.
ecod2 = np.zeros(JSON + 1)
em = np.zeros(JSON + 1)
ibeg = 0
iend = 0
hcsav = np.zeros((16, JSON + 1))
hkeep = np.zeros((4, JSON + 1))
# DATA NSEQ: the substep numbers of the 11 sequences (index 0 unused).
nseq = np.array([0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])


def bursmix(cod2, delts, hcomp, hd, hl, hp, hr, hs, hs1, hs2, hstot, ht,
            hv, imax, imin, jenv0, jenv1, lok, lcz, m, mrzone, mxzone,
            nrzone, nzone):
    """Alternate nuclear burning (with settling) and rotational mixing in
    sequences of 2, 4, ..., 20 substeps, and extrapolate the abundances
    from species 4 on to zero substep with bsrotmix; then restore X and Z
    to their values on entry and set Y = 1 - X - Z - X(He3).

    Translates: bursmix.f#BURSMIX (09c3b9b)

    Arguments: the mixing coefficients cod2; the time step delts (s); the
    composition hcomp[species, shell] after the caller's full step
    (changed in place); the model hd, hl, hp, hr, hs, hs1, hs2, hstot,
    ht, hv; imax, imin; jenv0, jenv1 (for the light elements); lok (from
    SECULR); lcz; m shells; the zones mrzone, mxzone, nrzone, nzone.
    oldmod_mod.hcompp is overwritten and left at the end-of-step
    composition.

    Returns jenv0: set to jenv1 when nmlpar_mod.lexcom, else as passed in.
    """
    global ibeg, iend
    n = nmlpar_mod
    hcompp = oldmod_mod.hcompp
    if n.lexcom:
        nspec = 15
    else:
        nspec = 11
    jbeg = 4
    jend = nspec
    nburs = 11
    # Store all elements: HCOMPP is the start-of-step composition in
    # KEMCOM, so it is overwritten during the substeps.
    for i in range(1, m + 1):
        for j in range(1, nspec + 1):
            hcsav[j, i] = hcompp[j, i]
        for j in range(1, 4):
            hkeep[j, i] = hcomp[j, i]
    # The caller's full step is the first entry of the extrapolation.
    iii = 1
    lconv = False
    lconv = bsrotmix(delts, hcomp, iii, m, jbeg, nspec, nseq, lconv)
    converged = False
    for iii in range(2, nburs + 1):
        # Reset the composition to its initial value.
        for i in range(1, m + 1):
            for j in range(1, nspec + 1):
                hcomp[j, i] = hcsav[j, i]
                hcompp[j, i] = hcsav[j, i]
        nstep = int(nseq[iii])
        dt = delts / float(nstep)
        for jjj in range(1, nstep + 1):
            # Nuclear burning (and settling).
            rotmix(dt, hcomp, hs2, ht, m, mrzone, mxzone, nrzone, nzone,
                   hstot, hd, hs, hr, hp, lcz, hs1)
            # Mixing.
            ibeg, iend = ndifcom(dt, cod2, ecod2, em, hd, hl, hp, hr, hs,
                                 hs1, hs2, hstot, hv, ibeg, iend, imax,
                                 imin, lcz, lok, m, hcomp, jbeg, jend)
            if n.lexcom:
                liburn2(dt, hcomp, hr, hs1, hs2, ht, jenv1, jenv0, m)
            for i in range(1, m + 1):
                for j in range(1, nspec + 1):
                    hcompp[j, i] = hcomp[j, i]
        lconv = bsrotmix(delts, hcomp, iii, m, jbeg, nspec, nseq, lconv)
        if lconv:
            for i in range(1, m + 1):
                for j in range(1, nspec + 1):
                    hcompp[j, i] = hcomp[j, i]
            converged = True
            break  # GOTO 10
    if not converged:
        # Failed to converge: the unextrapolated final composition is used.
        print(' WARNING - EXTRAPOLATION IN BSBURN DID NOT CONVERGE')
    # label 10
    # Question (see the page): X and Z restored to their values on entry.
    for i in range(1, m + 1):
        hcomp[1, i] = hkeep[1, i]
        hcomp[3, i] = hkeep[3, i]
        hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
    if n.lexcom:
        for i in range(1, m + 1):
            oldrat_mod.rli60[i] = newrat_mod.rli6[i]
            oldrat_mod.rli70[i] = newrat_mod.rli7[i]
            oldrat_mod.rbe90[i] = newrat_mod.rbe9[i]
        jenv0 = jenv1
        liov_mod.pscal0 = liov_mod.pscal1
    return jenv0
