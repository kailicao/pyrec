"""Python version of mixcom.f (YREC 7be4133): composition diffused across
one unstable region on the equally spaced grid in chi."""

import numpy as np

from pyrec import dimensions_mod
from pyrec import egrid_mod
from pyrec import splin_mod
from pyrec import tridi_mod
from pyrec.ccoeft import ccoeft
from pyrec.ctridi import ctridi
from pyrec.osplin import osplin

JSON = dimensions_mod.json

# The Fortran's locals (SAVE).
ecomp = np.zeros(JSON + 1)
id_ = np.zeros(16, dtype=int)
i0 = 0
i1 = 0


def mixcom(dt: float, ecod2: np.ndarray, em: np.ndarray, hs2: np.ndarray, ibeg: int, iend: int, lcz: np.ndarray,
           lok: bool, m: int, hcomp: np.ndarray, jbeg: int, jend: int) -> None:
    """Diffuse each species in rows jbeg..jend (helium, row 2, by closure)
    that varies over the region ibeg..iend: interpolate it to the equal
    grid in chi, solve the implicit diffusion equation there, interpolate
    the change back and add it; adjoining convection zones take the edge
    abundance if the edge zone is convective, or the end cell's change if
    it is radiative. Before the last iteration (lok False),
    Y = 1 - X - Z - He3.

    Translates: mixcom.f#MIXCOM (7be4133)

    Arguments: the step dt (s); ecod2 and em on the equal grid (from
    MIXGRID); the shell masses hs2; the region ibeg..iend; the convective
    flags lcz; lok; m shells; the model hcomp[k, j]; the rows jbeg..jend.
    The grid is in egrid_mod (chi, echi, dchi, ntot). Arrays have index 0
    unused.

    Sets in place: hcomp. Returns nothing.

    """

    global i0, i1
    chi = egrid_mod.chi
    echi = egrid_mod.echi
    xval = splin_mod.xval
    yval = splin_mod.yval
    xtab = splin_mod.xtab
    ytab = splin_mod.ytab
    u = tridi_mod.u

    # Which species vary over the unstable region (Fortran lines 59–74).
    nspec = 0
    for j in range(jbeg, jend + 1):
        test = hcomp[j, ibeg]
        # Row 2 (He4) is 1 - X - Z - He3, not solved for.
        if j != 2:
            for i in range(ibeg + 1, iend + 1):
                if abs(hcomp[j, i] - test) > 1.0e-14:
                    nspec += 1
                    id_[nspec] = j
                    break                                # GOTO 20
    if nspec == 0:
        return                                           # GOTO 200

    ntab = iend - ibeg + 1
    ntot = egrid_mod.ntot
    for n in range(1, nspec + 1):
        # The composition at the equal grid points at the start of the step (Fortran lines 77–87).
        for i in range(1, ntab + 1):
            ii = i + ibeg - 1
            xtab[i] = chi[i]
            ytab[i] = hcomp[id_[n], ii]
        osplin(echi, ecomp, xtab, ytab, ntab, ntot)

        # The implicit diffusion equation, solved for u in tridi_mod (Fortran lines 89–91).
        ccoeft(ecod2, egrid_mod.dchi, dt, ecomp, em, ntot)
        ctridi(ntot)

        # The change, not the new abundance, is interpolated back (Fortran lines 99–110).
        dcomp = 0.0
        for i in range(1, ntot + 1):
            xtab[i] = echi[i]
            ytab[i] = u[i] - ecomp[i]
            dcomp = dcomp + ytab[i] * em[i]
        for i in range(1, ntab + 1):
            xval[i] = chi[i]
        osplin(xval, yval, xtab, ytab, ntot, ntab)

        # A lower convection zone (Fortran lines 119–130).
        if ibeg > 1:
            for i in range(ibeg - 1, 0, -1):
                if not lcz[i]:
                    i0 = i + 1
                    break                                # GOTO 90
            else:
                i0 = 1
        else:
            i0 = ibeg

        # An upper convection zone (Fortran lines 132–143).
        if iend < m:
            for i in range(iend + 1, m + 1):
                if not lcz[i]:
                    i1 = i - 1
                    break                                # GOTO 97
            else:
                i1 = m
        else:
            i1 = iend

        # The species' mass before the update; the sums dcomp2, sumspec and
        # sumspec2 feed only a disabled check (unused?) (Fortran lines 144–149).
        dcomp2 = 0.0
        sumspec = 0.0
        for i in range(i0, i1 + 1):
            sumspec = sumspec + hs2[i] * hcomp[id_[n], i]

        # Update the composition (Fortran lines 151–157).
        sumspec2 = 0.0
        for i in range(1, ntab + 1):
            j = ibeg + i - 1
            hcomp[id_[n], j] = hcomp[id_[n], j] + yval[i]
            dcomp2 = dcomp2 + yval[i] * hs2[j]
            sumspec2 = sumspec2 + hs2[j] * hcomp[id_[n], j]

        # The inner convection zone, if any (Fortran lines 161–171).
        if i0 < ibeg:
            for i in range(ibeg - 1, i0 - 1, -1):
                if lcz[ibeg]:
                    hcomp[id_[n], i] = hcomp[id_[n], ibeg]
                else:
                    hcomp[id_[n], i] = hcomp[id_[n], i] + yval[1]
                sumspec2 = sumspec2 + hs2[i] * hcomp[id_[n], i]
                dcomp2 = dcomp2 + yval[1] * hs2[i]

        # The outer convection zone, if any (Fortran lines 173–183).
        if iend < i1:
            for i in range(iend + 1, i1 + 1):
                if lcz[iend]:
                    hcomp[id_[n], i] = hcomp[id_[n], iend]
                else:
                    hcomp[id_[n], i] = hcomp[id_[n], i] + yval[ntab]
                sumspec2 = sumspec2 + hs2[i] * hcomp[id_[n], i]
                dcomp2 = dcomp2 + yval[ntab] * hs2[i]
        # The conservation check (rescaling to sumspec) is commented out (Fortran lines 184–193).

    # He4 adjusted for the changes in X, Z and He3 (Fortran lines 196–200).
    if not lok:
        for i in range(i0, i1 + 1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
