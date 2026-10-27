"""Python version of reduce.f (YREC 09c3b9b): one step of the Henyey
elimination, for the pair of shells (i-1, i)."""

import numpy as np

from pyrec import constants_mod

# The work matrix Q(4,4) (SAVE; only columns 1-2 are used).
q = np.zeros((5, 5))


def reduce(i, hha, hhb, hl, hmax, hp, hr, hs, ht, qp0, qp, qpr0, qpr,
           qpp0, qpp, qt0, qt, qtr0, qtr, qtl0, qtl, qtp0, qtp, qtt0, qtt,
           qr0, qr, qrr0, qrr, qrp0, qrp, qrt0, qrt, ql0, ql, qlp0, qlp,
           qlt0, qlt):
    """Form the residuals of the four trapezoidal difference equations
    (P, T, R, L) between shells i-1 and i, and eliminate the unknowns of
    shell i-1, leaving x_k = HHB(k,i) - HHA(k,1,i) dP_i - HHA(k,2,i) dT_i
    for x = (dP_{i-1}, dT_{i-1}, dR_i, dL_i).

    Translates: reduce.f#REDUCE (09c3b9b)

    Arguments: i, the shell; hha[4+1, 2+1, JSON+1], hhb[4+1, JSON+1], set
    in place at shell i (and read at i-1); hl, hp, hr, hs, ht, the model;
    hmax[1..4], the largest residuals, updated in place; then each
    equation's right-hand side at i-1 (…0) and i, and their derivatives
    (qpr is d f_P / d R, and so on). Returns nothing.
    """
    # Residuals of the trapezoidal rule.
    sd05 = 0.5 * (hs[i] - hs[i - 1])
    hhb[1, i] = sd05 * (qp0 + qp) - (hp[i] - hp[i - 1])
    hmax[1] = max(hmax[1], abs(hhb[1, i]))
    hhb[2, i] = sd05 * (qt0 + qt) - (ht[i] - ht[i - 1])
    hmax[2] = max(hmax[2], abs(hhb[2, i]))
    hhb[3, i] = sd05 * (qr0 + qr) - (hr[i] - hr[i - 1])
    hmax[3] = max(hmax[3], abs(hhb[3, i]))
    hhb[4, i] = sd05 * (ql0 + ql) - (hl[i] - hl[i - 1])
    hmax[4] = max(hmax[4], abs(hhb[4, i]))
    sd05 = sd05 * constants_mod.cln

    # Columns 1 and 2: substitute dR_{i-1}, dL_{i-1} from the last call.
    pivot = -sd05 * qpr0
    q[1, 1] = (-1.0 - sd05 * qpp0) - pivot * hha[3, 1, i - 1]
    q[1, 2] = - pivot * hha[3, 2, i - 1]
    hhb[1, i] = hhb[1, i] - pivot * hhb[3, i - 1]
    pivot = -sd05 * qtr0
    pivotb = -sd05 * qtl0
    q[2, 1] = (-sd05 * qtp0 - pivot * hha[3, 1, i - 1]
               - pivotb * hha[4, 1, i - 1])
    q[2, 2] = ((-1.0 - sd05 * qtt0) - pivot * hha[3, 2, i - 1]
               - pivotb * hha[4, 2, i - 1])
    hhb[2, i] = hhb[2, i] - pivot * hhb[3, i - 1] - pivotb * hhb[4, i - 1]
    pivot = -1.0 - sd05 * qrr0
    q[3, 1] = -sd05 * qrp0 - pivot * hha[3, 1, i - 1]
    q[3, 2] = -sd05 * qrt0 - pivot * hha[3, 2, i - 1]
    hhb[3, i] = hhb[3, i] - pivot * hhb[3, i - 1]
    q[4, 1] = -sd05 * qlp0 + hha[4, 1, i - 1]
    q[4, 2] = -sd05 * qlt0 + hha[4, 2, i - 1]
    hhb[4, i] = hhb[4, i] + hhb[4, i - 1]
    hha[4, 1, i] = -sd05 * qlp
    hha[4, 2, i] = -sd05 * qlt

    # Pivot on row 4 and column 6 (dL_i).
    pivot = -sd05 * qtl
    q[2, 1] = q[2, 1] - pivot * q[4, 1]
    q[2, 2] = q[2, 2] - pivot * q[4, 2]
    hhb[2, i] = hhb[2, i] - pivot * hhb[4, i]
    hha[2, 1, i] = -sd05 * qtp - pivot * hha[4, 1, i]
    hha[2, 2, i] = (1.0 - sd05 * qtt) - pivot * hha[4, 2, i]

    # Pivot on row 3 and column 5 (dR_i).
    div = 1.0 / (1.0 - sd05 * qrr)
    q[3, 1] = q[3, 1] * div
    q[3, 2] = q[3, 2] * div
    hhb[3, i] = hhb[3, i] * div
    div = -sd05 * div
    hha[3, 1, i] = div * qrp
    hha[3, 2, i] = div * qrt
    pivot = -sd05 * qpr
    q[1, 1] = q[1, 1] - pivot * q[3, 1]
    q[1, 2] = q[1, 2] - pivot * q[3, 2]
    hhb[1, i] = hhb[1, i] - pivot * hhb[3, i]
    hha[1, 1, i] = (1.0 - sd05 * qpp) - pivot * hha[3, 1, i]
    hha[1, 2, i] = - pivot * hha[3, 2, i]
    pivot = -sd05 * qtr
    q[2, 1] = q[2, 1] - pivot * q[3, 1]
    q[2, 2] = q[2, 2] - pivot * q[3, 2]
    hhb[2, i] = hhb[2, i] - pivot * hhb[3, i]
    hha[2, 1, i] = hha[2, 1, i] - pivot * hha[3, 1, i]
    hha[2, 2, i] = hha[2, 2, i] - pivot * hha[3, 2, i]

    # Pivot on row 2 and column 4 (dT_{i-1}).
    # Question (see the page): no check for a zero pivot.
    div = 1.0 / q[2, 2]
    q[2, 1] = q[2, 1] * div
    hhb[2, i] = hhb[2, i] * div
    hha[2, 1, i] = hha[2, 1, i] * div
    hha[2, 2, i] = hha[2, 2, i] * div
    for j in range(1, 5):
        if j == 2:
            continue
        q[j, 1] = q[j, 1] - q[j, 2] * q[2, 1]
        hhb[j, i] = hhb[j, i] - q[j, 2] * hhb[2, i]
        hha[j, 1, i] = hha[j, 1, i] - q[j, 2] * hha[2, 1, i]
        hha[j, 2, i] = hha[j, 2, i] - q[j, 2] * hha[2, 2, i]

    # Pivot on row 1 and column 3 (dP_{i-1}).
    div = 1.0 / q[1, 1]
    hhb[1, i] = hhb[1, i] * div
    hha[1, 1, i] = hha[1, 1, i] * div
    hha[1, 2, i] = hha[1, 2, i] * div
    for j in range(2, 5):
        hhb[j, i] = hhb[j, i] - q[j, 1] * hhb[1, i]
        hha[j, 1, i] = hha[j, 1, i] - q[j, 1] * hha[1, 1, i]
        hha[j, 2, i] = hha[j, 2, i] - q[j, 1] * hha[1, 2, i]
