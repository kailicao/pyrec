"""Python version of bsrotmix.f (YREC 09c3b9b): extrapolation of burn-mix
sequences to zero substep."""

import math

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Local arrays and scalars kept between calls by the Fortran's blanket SAVE
# (the tableau HCOMPI, HCOMPA, HCOMPE was COMMON/BSBURN/ until 09c3b9b).
# Index 0 of every dimension is unused.
_hcompi = np.zeros((12, 16, _JSON + 1))
_hcompa = np.zeros((16, _JSON + 1))
_hcompe = np.zeros((16, _JSON + 1))
_d = np.zeros((16, _JSON + 1))
_x = np.zeros(12)
_jj = np.zeros(16, dtype=int)
_ldo = np.zeros(16, dtype=bool)
_lcno = np.zeros(_JSON + 1, dtype=bool)
_lhe3 = np.zeros(_JSON + 1, dtype=bool)
_saved = {'nmax': 0, 'lcncheck': False}


def _transform(target, hcomp, i):
    """Fill the CNO variables 5, 6, 7, 9 of target[:, i]: C13/C12, N14/C12,
    O16/N14 and C12/12 + C13/13 + N14/14 + O16/16 (written out three times
    in the Fortran)."""
    x12 = hcomp[5, i]
    x13 = hcomp[6, i]
    x14 = hcomp[7, i]
    x16 = hcomp[9, i]
    target[5, i] = x13 / x12
    target[6, i] = x14 / x12
    target[7, i] = x16 / x14
    target[9, i] = x12 / 12.0 + x13 / 13.0 + x14 / 14.0 + x16 / 16.0


def bsrotmix(delts, hcomp, iest, m, jbeg, nspec, nseq, lconv):
    """Polynomial (Richardson) extrapolation to zero substep of the
    abundances from successive burn-mix sequences, with He3 in logarithm
    and CNO as ratios plus a total; sets lconv when the largest relative
    correction is below 1e-3.

    Translates: bsrotmix.f#BSROTMIX (09c3b9b)

    Arguments: the time step delts; the abundances hcomp[species, shell]
    after sequence iest (changed in place to the extrapolated abundances
    once converged); the number of shells m; the species range jbeg..nspec;
    the substep numbers nseq.

    Returns lconv: True once converged; otherwise the value passed in
    (False after the first sequence).
    """
    hcompi, hcompa, hcompe, d, x = _hcompi, _hcompa, _hcompe, _d, _x
    jj, ldo, lcno, lhe3 = _jj, _ldo, _lcno, _lhe3

    # Determine which species require calculation.
    if iest == 1:
        for j in range(jbeg, nspec + 1):
            ldo[j] = False
            for i in range(m, 0, -1):
                if hcomp[j, i] > 1.0e-14:
                    ldo[j] = True
                    break  # to label 5
        if ldo[4]:
            for i in range(1, m + 1):
                lhe3[i] = hcomp[4, i] > 1.0e-24
        # Alternate method for CNO: solve for C13/C12, N14/C12, O16/N14 and
        # the sum of CNO nuclei rather than the individual abundances.
        if ldo[5] and ldo[6] and ldo[7] and ldo[9]:
            _saved['lcncheck'] = True
            for i in range(1, m + 1):
                lcno[i] = (hcomp[5, i] > 1.0e-24 and hcomp[6, i] > 1.0e-24
                           and hcomp[7, i] > 1.0e-24
                           and hcomp[9, i] > 1.0e-24)
        else:
            _saved['lcncheck'] = False
            for i in range(1, m + 1):
                lcno[i] = False
        icount = 1
        for i in range(jbeg, nspec + 1):
            if ldo[i]:
                jj[icount] = i
                icount = icount + 1
        _saved['nmax'] = icount - 1
        # Question (see the page): the run stops without a message.
        if _saved['nmax'] < 1:
            raise SystemExit(999)
    nmax = _saved['nmax']
    lcncheck = _saved['lcncheck']

    # Store current results in the composition (hcompa) and the error
    # (hcompe; the Fortran comment says DCOMPA).
    for j1 in range(1, nmax + 1):
        j = jj[j1]
        for i in range(1, m + 1):
            hcompa[j, i] = hcomp[j, i]
            hcompe[j, i] = hcomp[j, i]
    # Extrapolate in the log of He3 rather than the abundance.
    if ldo[4]:
        for i in range(1, m + 1):
            if lhe3[i]:
                hcompa[4, i] = math.log(hcomp[4, i])
                hcompe[4, i] = hcompa[4, i]
    if lcncheck:
        for i in range(1, m + 1):
            if lcno[i]:
                _transform(hcompa, hcomp, i)
                for k in (5, 6, 7, 9):
                    hcompe[k, i] = hcompa[k, i]

    x[iest] = (delts / float(nseq[iest]))**2
    if iest == 1:
        # Current results are the initial extrapolation; exit.
        for j1 in range(1, nmax + 1):
            j = jj[j1]
            for i in range(1, m + 1):
                hcompi[iest, j, i] = hcomp[j, i]
        if ldo[4]:
            for i in range(1, m + 1):
                if lhe3[i]:
                    hcompi[iest, 4, i] = math.log(hcomp[4, i])
        if lcncheck:
            for i in range(1, m + 1):
                if lcno[i]:
                    for k in (5, 6, 7, 9):
                        hcompi[iest, k, i] = hcompa[k, i]
        lconv = False
    else:
        # Polynomial extrapolation (Neville's scheme, as in PZEXTR).
        xest = x[iest]
        for j1 in range(1, nmax + 1):
            j = jj[j1]
            for i in range(1, m + 1):
                d[j, i] = hcomp[j, i]
        if ldo[4]:
            for i in range(1, m + 1):
                if lhe3[i]:
                    d[4, i] = math.log(hcomp[4, i])
        if lcncheck:
            for i in range(1, m + 1):
                if lcno[i]:
                    _transform(d, hcomp, i)
        for k1 in range(1, iest):
            delta = 1.0 / (x[iest - k1] - xest)
            f1 = xest * delta
            f2 = x[iest - k1] * delta
            for j1 in range(1, nmax + 1):
                j = jj[j1]
                for i in range(1, m + 1):
                    q = hcompi[k1, j, i]
                    hcompi[k1, j, i] = hcompe[j, i]
                    delta = d[j, i] - q
                    hcompe[j, i] = f1 * delta
                    d[j, i] = f2 * delta
                    hcompa[j, i] = hcompa[j, i] + hcompe[j, i]
        for j1 in range(1, nmax + 1):
            j = jj[j1]
            for i in range(1, m + 1):
                hcompi[iest, j, i] = hcompe[j, i]

        # Now check if the error is within tolerance.
        errmax = 1.0e-30
        imax = 1
        jmax = 1
        for j1 in range(1, nmax + 1):
            j = jj[j1]
            for i in range(1, m + 1):
                # Question (see the page): ln X3 < 0 and the ratios are
                # tested on their own size.
                if hcompa[j, i] > 1.0e-12:
                    q = hcompe[j, i] / hcompa[j, i]
                    if abs(q) > errmax:
                        imax = i
                        jmax = j
                        errmax = abs(q)
        print(f"     ITER{iest:3d} MAX ERR {errmax:10.2E} SPECIES {jmax:3d}"
              f" SHELL {imax:5d} ABUND{hcomp[jmax, imax]:12.4E}"
              f"{hcompa[jmax, imax]:12.4E}")  # format 10
        # Question (see the page): the tolerance 1e-3 is hard-coded.
        if errmax < 1.0e-3:
            lconv = True
        if lconv:
            # Convert back to absolute abundances.
            if ldo[4]:
                for i in range(1, m + 1):
                    if lhe3[i]:
                        hcompa[4, i] = math.exp(hcompa[4, i])
            if lcncheck:
                for i in range(1, m + 1):
                    if lcno[i]:
                        total = (1.0 / 12.0 + hcompa[5, i] / 13.0
                                 + hcompa[6, i] / 14.0
                                 + hcompa[6, i] * hcompa[7, i] / 16.0)
                        x12 = hcompa[9, i] / total
                        x13 = hcompa[5, i] * x12
                        x14 = hcompa[6, i] * x12
                        x16 = hcompa[7, i] * x14
                        hcompa[5, i] = x12
                        hcompa[6, i] = x13
                        hcompa[7, i] = x14
                        hcompa[9, i] = x16
            for j1 in range(1, nmax + 1):
                j = jj[j1]
                for i in range(1, m + 1):
                    hcomp[j, i] = hcompa[j, i]
                    if hcomp[j, i] < 1.0e-24:
                        hcomp[j, i] = 0.0
    return lconv
