"""Python version of rotgrid.f (YREC 7be4133): one radiative region mapped
onto an equally spaced grid for the transport of angular momentum."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import difad2_mod
from pyrec import difad_mod
from pyrec import dimensions_mod
from pyrec import egrid_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import splin_mod
from pyrec.getgrid import getgrid
from pyrec.osplin import osplin

# The local array EI0 (SAVE): I/(m r^2) at the new points.
ei0 = np.zeros(dimensions_mod.json + 1)


def _edge_table(ytab: np.ndarray, values: np.ndarray, ibeg: int, iend: int, ntab: int,
                fac: float = 1.0) -> None:
    """The edge values for OSPLIN: fac*values at IBEG+1, IBEG+1..IEND
    and IEND, beside XTAB's ends and half-way points."""

    ytab[1] = fac * values[ibeg + 1]
    for i in range(2, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = fac * values[ii]
    ytab[ntab + 1] = fac * values[iend]


def rotgrid(cod: np.ndarray, cod2: np.ndarray, hd: np.ndarray, hi: np.ndarray, hjmsav: np.ndarray, hl: np.ndarray,
            hp: np.ndarray, hr: np.ndarray, hru: np.ndarray, hs: np.ndarray, hs1: np.ndarray, hs2: np.ndarray,
            hstot: float, ibeg: int, iend: int, lcz: np.ndarray, m: int, wsav: np.ndarray, dr: float, ecod: np.ndarray,
            ecod2: np.ndarray, ei: np.ndarray, ej: np.ndarray, em: np.ndarray, ew: np.ndarray,
            ldum2: bool) -> tuple[float, bool]:
    """Map the radiative region ibeg..iend onto points equally spaced in
    the mesh coordinate chi (GETGRID), and compute there the zones'
    masses, moments of inertia, angular momenta and angular velocities
    (at the centers) and the diffusion coefficients times the Jacobian
    rho r^2 dchi/dr (at the edges). End zones take in any neighboring
    convection zone whole.

    Translates: rotgrid.f#ROTGRID (7be4133)

    Arguments: the model's cod, cod2, hd, hi, hjmsav, hl, hp, hr, hru, hs,
    hs1, hs2 (index 0 unused); hstot; ibeg, iend, the region; lcz; m;
    wsav; dr (set); ecod, ecod2, ei, ej, em, ew, set in place on the new
    grid; ldum2 (set). Returns (dr, ldum2): dr = DCHI, and ldum2 True
    (nothing done) when the region has at most one interface.

    """

    c = constants_mod
    nml = nmlpar_mod
    xval = splin_mod.xval
    yval = splin_mod.yval
    xtab = splin_mod.xtab
    ytab = splin_mod.ytab

    # A single unstable interface: nothing to do (Fortran lines 66–73).
    if iend - ibeg <= 1:
        ldum2 = True
        return dr, ldum2        # goto 9999
    ldum2 = False

    # The grid of equally spaced points (Fortran lines 74–75).
    getgrid(hl, hp, hs, ibeg, iend, m)
    chi = egrid_mod.chi
    echi = egrid_mod.echi
    es1 = egrid_mod.es1
    ntot = egrid_mod.ntot
    ntab = iend - ibeg + 1

    # 1. Zone centers: log r, I/(m r^2), j and Omega on the new points
    # (Fortran lines 76–113).
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        xtab[i] = chi[i]
        ytab[i] = hr[ii]
    for i in range(1, ntot + 1):
        xval[i] = echi[i]
    osplin(xval, yval, xtab, ytab, ntab, ntot)
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = hi[ii] / (hs2[ii] * hru[ii]**2)
    osplin(xval, ei0, xtab, ytab, ntab, ntot)
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = hjmsav[ii]
    osplin(xval, ej, xtab, ytab, ntab, ntot)
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = wsav[ii]
    osplin(xval, ew, xtab, ytab, ntab, ntot)

    # Total I and J of the interior zones (Fortran lines 114–122).
    for i in range(2, ntot):
        em[i] = 0.5 * (es1[i + 1] - es1[i - 1])
        ei[i] = ei0[i] * em[i] * math.exp(c.cln * 2.0 * yval[i])
        ej[i] = ej[i] * em[i]
    i = max(2, ntot)    # the Fortran's I after the loop

    # 2. The inner end zone, with any convection zone below it (Fortran
    # lines 123–154).
    emtop = 0.5 * (es1[2] + es1[1])
    if ibeg > 1:
        embot = 0.5 * (hs1[ibeg] + hs1[ibeg - 1])
    else:
        embot = 0.0
    em[1] = emtop - embot
    # r^2 at the shell IBEG, the same as the interpolated value.
    ei[1] = ei0[1] * em[1] * hru[ibeg]**2
    ej[1] = ej[1] * em[1]
    if ibeg > 1:
        for ii in range(ibeg - 1, 0, -1):
            if not lcz[ii]:
                # I, not II; i0 is unused?
                i0 = i + 1  # noqa: F841  (unused in the Fortran too)
                break           # goto 10
            em[1] = em[1] + hs2[ii]
            ei[1] = ei[1] + hi[ii]
            ej[1] = ej[1] + hjmsav[ii] * hs2[ii]
        else:
            i0 = 1
    else:
        i0 = 1  # noqa: F841  (unused in the Fortran too)

    # The outer end zone, with any convection zone above it (Fortran
    # lines 155–181).
    embot = 0.5 * (es1[ntot] + es1[ntot - 1])
    if iend < m:
        emtop = 0.5 * (hs1[iend] + hs1[iend + 1])
    else:
        emtop = math.exp(c.cln * hstot)
    em[ntot] = emtop - embot
    ei[ntot] = ei0[ntot] * em[ntot] * hru[iend]**2
    ej[ntot] = ej[ntot] * em[ntot]
    if iend < m:
        for ii in range(iend + 1, m + 1):
            if not lcz[ii]:
                i1 = i - 1  # noqa: F841  (unused in the Fortran too; I, not II)
                break           # goto 20
            em[ntot] = em[ntot] + hs2[ii]
            ei[ntot] = ei[ntot] + hi[ii]
            ej[ntot] = ej[ntot] + hjmsav[ii] * hs2[ii]
        else:
            i1 = m
    else:
        i1 = m  # noqa: F841  (unused in the Fortran too)

    # 3. Zone edges: the diffusion coefficients at chi - dchi/2 (Fortran
    # lines 182–211).
    xtab[1] = chi[1]
    for i in range(2, ntab + 1):
        xtab[i] = 0.5 * (chi[i] + chi[i - 1])
    ntabb = ntab + 1
    xtab[ntabb] = chi[ntab]
    _edge_table(ytab, cod, ibeg, iend, ntab)
    xval[1] = chi[1]
    for i in range(2, ntot + 1):
        xval[i] = echi[i] - 0.5 * egrid_mod.dchi
    osplin(xval, ecod, xtab, ytab, ntabb, ntot)
    _edge_table(ytab, cod2, ibeg, iend, ntab)
    osplin(xval, ecod2, xtab, ytab, ntabb, ntot)

    # With LDIFAD, the advective and diffusive terms; unused? (DADCOEFT
    # sets ecod3, ecod4 afresh) (Fortran lines 212–238).
    if nml.ldifad:
        fac = 0.2 * c.c4pi * nml.fw
        _edge_table(ytab, difad2_mod.vesa, ibeg, iend, ntab, fac)
        osplin(xval, difad_mod.ecod3, xtab, ytab, ntabb, ntot)
        _edge_table(ytab, difad2_mod.vesd, ibeg, iend, ntab, fac)
        osplin(xval, difad_mod.ecod4, xtab, ytab, ntabb, ntot)
        for i in range(1, ntot + 1):
            difad_mod.ecod4[i] = difad_mod.ecod4[i] + ecod2[i]

    # 4. The Jacobian rho r^2 dchi/dr, interpolated in log10; question
    # (see the page): the comment omits P (Fortran lines 239–255).
    dm = nml.hpttol[2]
    dl = nml.hpttol[9] * hl[m] * c.clsun
    dp = nml.hpttol[11]
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        xtab[i] = chi[i]
        fac = c.c4pi * math.exp(c.cln * (hd[ii] + 2.0 * hr[ii]))
        qchir = (fac / (c.cln * dm * hs1[ii]) + fac * mdphy_mod.epsm[ii] / dl
                 + math.exp(c.cln * (c.cgl + hd[ii] + hs[ii] - hp[ii]
                                     - 2.0 * hr[ii])) / (c.cln * dp))
        ytab[i] = hd[ii] + math.log10(qchir) + 2.0 * hr[ii]
    osplin(xval, yval, xtab, ytab, ntab, ntot)

    # The factor for the composition coefficient (Fortran lines 256–261).
    for i in range(1, ntot + 1):
        ecod2[i] = ecod2[i] * math.exp(c.cln * yval[i])

    # rho r^2 (I/m) dchi/dr for the angular-momentum coefficients
    # (Fortran lines 262–280).
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = ytab[i] + math.log10(hi[ii] / hs2[ii])
    osplin(xval, yval, xtab, ytab, ntab, ntot)
    for i in range(1, ntot + 1):
        ecod[i] = ecod[i] * math.exp(c.cln * yval[i])
    if nml.ldifad:
        for i in range(1, ntot + 1):
            difad_mod.ecod3[i] = difad_mod.ecod3[i] * math.exp(c.cln * yval[i])
            difad_mod.ecod4[i] = difad_mod.ecod4[i] * math.exp(c.cln * yval[i])

    # DR is DCHI (Fortran lines 281–283; label 9999).
    dr = egrid_mod.dchi
    return dr, ldum2
