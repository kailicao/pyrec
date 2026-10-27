"""Python version of eq2mod.f (YREC 09c3b9b): a change of specific angular
momentum on the even grid, mapped back onto the model's points."""

from pyrec import egrid_mod
from pyrec import splin_mod
from pyrec.osplin import osplin


def eq2mod(dj, ej, hs2, ibeg, iend, lcz, m, sumdj, hjm):
    """Interpolate Delta j / j from the even grid back to the model points
    of the region ibeg..iend (extended to the edges of convection zones it
    touches), apply it to hjm, and rescale the region so that its angular
    momentum changes by exactly sumdj.

    Translates: eq2mod.f#EQ2MOD (09c3b9b)

    Arguments: dj, ej, the change and the angular momentum on the even
    grid; hs2, the shell masses; ibeg, iend, the region; lcz, the
    convective flags; m, the number of shells; sumdj, the total change;
    hjm, the specific angular momentum, changed in place.

    Returns nothing.
    """
    chi = egrid_mod.chi
    echi = egrid_mod.echi
    ntot = egrid_mod.ntot
    xval = splin_mod.xval
    yval = splin_mod.yval
    xtab = splin_mod.xtab
    ytab = splin_mod.ytab

    # Extend the region to the edges of the convection zones at its ends.
    if not lcz[ibeg] or ibeg == 1:
        i0 = ibeg
    else:
        i0 = 1
        for i in range(ibeg - 1, 0, -1):
            if not lcz[i]:
                i0 = i + 1
                break  # label 10
    if not lcz[iend] or iend == m:
        i1 = iend
    else:
        i1 = m
        for i in range(iend + 1, m + 1):
            if not lcz[i]:
                i1 = i - 1
                break  # label 20

    # Interpolate y = dj/ej from the even grid to the model points.
    nmod = iend - ibeg + 1
    for i in range(1, nmod + 1):
        xval[i] = chi[i]
    for i in range(1, ntot + 1):
        xtab[i] = echi[i]
        ytab[i] = dj[i] / ej[i]
    osplin(xval, yval, xtab, ytab, ntot, nmod)

    # Apply it; the extended convective shells take the edge values.
    # Question (see the page): this assumes those zones rotate as solid
    # bodies already.
    sumjmod = 0.0
    sumdjmod = 0.0
    if i0 < ibeg:
        for i in range(i0, ibeg):
            sumjmod = sumjmod + hjm[i] * hs2[i]
            sumdjmod = sumdjmod + yval[1] * hjm[i] * hs2[i]
            hjm[i] = hjm[i] * (1.0 + yval[1])
    for i in range(ibeg, iend + 1):
        ii = i - ibeg + 1
        sumjmod = sumjmod + hjm[i] * hs2[i]
        sumdjmod = sumdjmod + yval[ii] * hjm[i] * hs2[i]
        hjm[i] = hjm[i] * (1.0 + yval[ii])
    if i1 > iend:
        for i in range(iend + 1, i1 + 1):
            sumjmod = sumjmod + hjm[i] * hs2[i]
            sumdjmod = sumdjmod + yval[nmod] * hjm[i] * hs2[i]
            hjm[i] = hjm[i] * (1.0 + yval[nmod])

    # Rescale to conserve angular momentum.
    # Question (see the page): with sumdj = 0, any mismatch triggers it.
    test = abs(sumdjmod - sumdj)
    val = 1.0e-5 * abs(sumdj)
    if test > val:
        ratio = (sumjmod + sumdj) / (sumjmod + sumdjmod)
        for i in range(i0, i1 + 1):
            hjm[i] = hjm[i] * ratio
