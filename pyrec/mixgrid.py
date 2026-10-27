"""Python version of mixgrid.f (YREC 09c3b9b): the equally spaced grid in
chi for mixing one unstable region, with its zone masses and coefficients,
for MIXCOM."""

import math

from pyrec import constants_mod
from pyrec import egrid_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import splin_mod
from pyrec.getgrid import getgrid
from pyrec.osplin import osplin

# The Fortran's locals (SAVE).
i0 = 0
i1 = 0


def mixgrid(cod2, hd, hl, hp, hr, hs, hs1, hs2, hstot, ibeg, iend, lcz, m,
            ecod2, em, ldum2):
    """Set up the equal grid in chi over the region ibeg..iend (GETGRID),
    the mass em of each grid zone (the end zones with any adjoining
    convection zone), and the coefficient ecod2: cod2 interpolated to the
    zone edges times rho r^2 dchi/dr. Does nothing but set ldum2 = True
    when the region spans a single unstable interface.

    Translates: mixgrid.f#MIXGRID (09c3b9b)

    Arguments: the mixing coefficients cod2 at the model interfaces; the
    model hd, hl, hp, hr, hs (logs), hs1 (mass, g), hs2 (shell masses),
    hstot (log total mass); the region ibeg..iend; the flags lcz; m shells.
    ldum2 is passed in with any value. Arrays have index 0 unused.

    Sets in place: ecod2 and em on the egrid_mod.ntot grid points.

    Returns ldum2.
    """
    global i0, i1
    c = constants_mod
    chi = egrid_mod.chi
    echi = egrid_mod.echi
    es1 = egrid_mod.es1
    xval = splin_mod.xval
    yval = splin_mod.yval
    xtab = splin_mod.xtab
    ytab = splin_mod.ytab
    # The special case of a single unstable interface.
    if iend - ibeg <= 1:
        ldum2 = True
        return ldum2                                     # GOTO 9999
    ldum2 = False
    # The grid of equally spaced points (chi, echi, es1, dchi, ntot).
    getgrid(hl, hp, hs, ibeg, iend, m)
    ntot = egrid_mod.ntot
    dchi = egrid_mod.dchi
    # Zone masses, intermediate points.
    for i in range(2, ntot):
        em[i] = 0.5 * (es1[i + 1] - es1[i - 1])
    i = max(ntot, 2)   # the Fortran's I after DO I = 2, NTOT-1
    # The centre: halfway to the next grid point, plus any convection zone
    # below.
    emtop = 0.5 * (es1[2] + es1[1])
    if ibeg > 1:
        embot = 0.5 * (hs1[ibeg] + hs1[ibeg - 1])
    else:
        embot = 0.0
    em[1] = emtop - embot
    if ibeg > 1:
        for ii in range(ibeg - 1, 0, -1):
            if not lcz[ii]:
                # Question (see the page): I, not II (I0 is unused).
                i0 = i + 1
                break                                    # GOTO 10
            em[1] = em[1] + hs2[ii]
        else:
            i0 = 1
    else:
        i0 = 1
    # The surface.
    embot = 0.5 * (es1[ntot] + es1[ntot - 1])
    if iend < m:
        emtop = 0.5 * (hs1[iend] + hs1[iend + 1])
    else:
        emtop = math.exp(c.cln * hstot)
    em[ntot] = emtop - embot
    if iend < m:
        for ii in range(iend + 1, m + 1):
            if not lcz[ii]:
                # Question (see the page): I, not II (I1 is unused).
                i1 = i - 1
                break                                    # GOTO 20
            em[ntot] = em[ntot] + hs2[ii]
        else:
            i1 = m
    else:
        i1 = m
    # The coefficient at the zone edges, constant beyond the end interfaces.
    ntab = iend - ibeg + 1
    xtab[1] = chi[1]
    for i in range(2, ntab + 1):
        xtab[i] = 0.5 * (chi[i] + chi[i - 1])
    ntabb = ntab + 1
    xtab[ntabb] = chi[ntab]
    ytab[1] = cod2[ibeg + 1]
    for i in range(2, ntab + 1):
        ii = ibeg + i - 1
        ytab[i] = cod2[ii]
    ytab[ntabb] = cod2[iend]
    xval[1] = chi[1]
    for i in range(2, ntot + 1):
        xval[i] = echi[i] - 0.5 * dchi
    osplin(xval, ecod2, xtab, ytab, ntabb, ntot)
    # log10 of rho r^2 dchi/dr, with FAC = 4 pi rho r^2.
    hpttol = nmlpar_mod.hpttol
    dm = hpttol[2]
    dl = hpttol[9] * hl[m] * c.clsun
    dp = hpttol[11]
    for i in range(1, ntab + 1):
        ii = ibeg + i - 1
        xtab[i] = chi[i]
        fac = c.c4pi * math.exp(c.cln * (hd[ii] + 2.0 * hr[ii]))
        # Question (see the page): the comment omits P in the last term.
        qchir = (fac / (c.cln * dm * hs1[ii]) + fac * mdphy_mod.epsm[ii] / dl
                 + math.exp(c.cln * (c.cgl + hd[ii] + hs[ii] - hp[ii]
                                     - 2.0 * hr[ii])) / (c.cln * dp))
        ytab[i] = hd[ii] + math.log10(qchir) + 2.0 * hr[ii]
    osplin(xval, yval, xtab, ytab, ntab, ntot)
    # The multiplicative factors.
    for i in range(1, ntot + 1):
        ecod2[i] = ecod2[i] * math.exp(c.cln * yval[i])
    return ldum2
