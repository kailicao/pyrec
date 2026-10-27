"""Python version of eqscvg.f (YREC development, 09c3b9b): the SCV
equation of state for a general composition, from the pure-H, pure-He
and metal tables."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import scveos_mod
from pyrec.interp import interp

NTS = 63                  # temperatures in the SCV tables

# Locals kept between calls (SAVE).
qr = np.zeros(5)
ft = np.zeros(5)
ftd = np.zeros(5)
fp = np.zeros(5)
fpd = np.zeros(5)
temp = np.zeros((5, 4))


def _rows(table, col, ncol=3):
    """The block the Fortran writes out for each quantity: Lagrange
    interpolation in log Pgas (weights fp) of column col of a table, at
    the 4 temperatures of the stencil. Returns the values at rows 1-4.
    ncol=2 for a table with no column index (SMIX)."""
    scv = scveos_mod
    out = np.zeros(5)
    for i in range(1, 5):
        ii = scv.idtscv + i - 1
        idp = scv.idp
        if ncol == 2:
            out[i] = sum(fp[k] * table[ii, idp + k - 1] for k in range(1, 5))
        else:
            out[i] = sum(fp[k] * table[ii, idp + k - 1, col]
                         for k in range(1, 5))
    return out


def _in_t(vals, w=None):
    """Interpolate in temperature: sum of the weights ft (or w) times the
    values at the 4 temperatures."""
    w = ft if w is None else w
    return w[1] * vals[1] + w[2] * vals[2] + w[3] * vals[3] + w[4] * vals[4]


def eqscvg(tl, t, p, dl, d, x, z, beta, fxion, qdt, qdp, qcp, dela, lcalc):
    """The SCV equation of state for any X, Z: density by additive volumes
    of the H, He and metal tables, derivatives as volume-weighted means,
    c_P from the tabulated entropies, then radiation and the H and He
    ionization fractions (equations on the page). Called by eqscve.

    Translates: eqscvg.f#EQSCVG (09c3b9b)

    Arguments: tl, t (log T, T), p (total pressure), x, z, beta (Pgas/P);
    dl, d, qdt, qdp, qcp, dela and lcalc are only set. Array changed in
    place: fxion (H+, He+, He++ fractions, index 0 unused).

    Returns (dl, d, qdt, qdp, qcp, dela, lcalc). Outside the table, lcalc
    is False and the others are returned as passed in.
    """
    c = constants_mod
    scv = scveos_mod
    # The gas pressure is the tabulated quantity.
    betm1 = 1.0 - beta
    pp = math.log10(beta * p)
    # Is the point within the table?
    if tl < scv.tlogx[1] or tl > scv.tlogx[NTS] or pp < 4.0:
        lcalc = False
        return dl, d, qdt, qdp, qcp, dela, lcalc

    # Nearest 4 points in temperature, searching from the last stencil
    # (left by this routine or eqscve).
    if tl < scv.tlogx[scv.idtscv + 1]:
        for i in range(scv.idtscv, 0, -1):
            if tl > scv.tlogx[i]:
                ii = i - 1
                break                              # go to 10
        else:
            ii = 1
    else:
        for i in range(scv.idtscv + 2, NTS + 1):
            if tl < scv.tlogx[i]:
                ii = i - 2
                break                              # go to 20
        else:
            ii = NTS - 3
    scv.idtscv = max(1, ii)
    scv.idtscv = min(NTS - 3, scv.idtscv)
    idt = scv.idtscv

    # Nearest 4 points in pressure, searched in the row idt (the lowest
    # temperature of the stencil, with the smallest range in P).
    nptsx = scv.nptsx[idt]
    jjj = min(nptsx - 3, scv.idp)
    if pp < scv.tablenv[idt, jjj + 1, 1]:
        for j in range(jjj, 0, -1):
            if pp > scv.tablenv[idt, j, 1]:
                jj = j - 1
                break                              # go to 30
        else:
            jj = 1
        scv.idp = max(1, jj)
        scv.idp = min(nptsx - 3, scv.idp)
    else:
        for j in range(jjj + 2, nptsx + 1):
            if pp < scv.tablenv[idt, j, 1]:
                jj = j - 2
                break                              # go to 40
        else:
            # The point is above the table: return.
            lcalc = False
            return dl, d, qdt, qdp, qcp, dela, lcalc
        scv.idp = min(nptsx - 3, jj)               # label 40
    idp = scv.idp
    lcalc = True

    # Lagrange weights in log T and log Pgas.
    for k in range(1, 5):
        qr[k] = scv.tlogx[k + idt - 1]
    interp(qr, ft, ftd, tl)
    for k in range(1, 5):
        # Question (see the page): row idt's pressures for all 4 rows.
        qr[k] = scv.tablenv[idt, idp + k - 1, 1]
    interp(qr, fp, fpd, pp)
    y = 1.0 - x - z
    # Radiation pressure.
    prad = betm1 * p
    pgas = beta * p

    # Densities for X=1, Y=1, Z=1 (column 4; log10 in the tables).
    denxl = _in_t(_rows(scv.tablex, 4))
    denx = math.exp(c.cln * denxl)
    denyl = _in_t(_rows(scv.tabley, 4))
    deny = math.exp(c.cln * denyl)
    denzl = _in_t(_rows(scv.tablez, 4))
    denz = math.exp(c.cln * denzl)
    # Density by additive volumes, and the volume weights.
    d = 1.0 / (x / denx + y / deny + z / denz)
    dl = math.log10(d)
    ddx = x * d / denx
    ddy = y * d / deny
    ddz = z * d / denz

    # d ln rho / d ln P (gas): columns 8, 8 and 13.
    qdpx = _in_t(_rows(scv.tablex, 8))
    qdpy = _in_t(_rows(scv.tabley, 8))
    qdpz = _in_t(_rows(scv.tablez, 13))
    qdp0 = qdpx * ddx + qdpy * ddy + qdpz * ddz
    qdp = qdp0 / beta

    # d ln rho / d ln T (gas): columns 7, 7 and 10; then radiation.
    qdtx = _in_t(_rows(scv.tablex, 7))
    qdty = _in_t(_rows(scv.tabley, 7))
    qdtz = _in_t(_rows(scv.tablez, 10))
    qdt0 = qdtx * ddx + qdty * ddy + qdtz * ddz
    qpt0 = -qdt0 / qdp0
    qpt = qpt0 * beta + 4.0 * betm1
    qdt = -qpt * qdp

    # Entropies (column 5, log10) and the entropy of mixing.
    temp[:, 1] = _rows(scv.tablex, 5)
    temp[:, 2] = _rows(scv.tabley, 5)
    temp[:, 3] = _rows(scv.smix, None, ncol=2)
    sxl = _in_t(temp[:, 1])
    sx = math.exp(c.cln * sxl)
    syl = _in_t(temp[:, 2])
    sy = math.exp(c.cln * syl)
    sm = _in_t(temp[:, 3])
    # Question (see the page): log10 of SMIX fails if it is <= 0.
    qsmt = (ftd[1] * math.log10(temp[1, 3]) + ftd[2] * math.log10(temp[2, 3])
            + ftd[3] * math.log10(temp[3, 3])
            + ftd[4] * math.log10(temp[4, 3]))
    s0 = x * sx + y * sy + sm                      # unused, as in the Fortran

    # d ln S / d ln T for X and Y (column 9), dU/dT for Z (column 7, log10).
    qstx = _in_t(_rows(scv.tablex, 9))
    qsty = _in_t(_rows(scv.tabley, 9))
    qutzl = _in_t(_rows(scv.tablez, 7))
    qutz = math.exp(c.cln * qutzl)
    qst0 = x * sx * qstx + y * sy * qsty + sm * qsmt
    qcp0 = qst0 + z * (qutz + pgas * qdtz ** 2 / qdpz / denz / t)
    # dU/dT for the gas, then corrected for radiation.
    qut = qcp0 - pgas * qdt0 ** 2 / qdp0 / d / t
    ur = 3.0 * prad / d
    qut = qut + 4.0 * ur / t
    # c_P with radiation pressure, and the adiabatic gradient.
    qcp = qut + p * qdt ** 2 / qdp / d / t
    dela = -p * qdt / d / t / qcp

    # Fractions of all particles (electrons included): neutral H2 and
    # neutral He (column 2), neutral H atoms and He+ (column 3).
    xtf_h2 = _in_t(_rows(scv.tablex, 2))
    xtf_he = _in_t(_rows(scv.tabley, 2))
    xtf_h1 = _in_t(_rows(scv.tablex, 3))
    xtf_hep = _in_t(_rows(scv.tabley, 3))
    # The rest from conservation of particles and of charge.
    xtf_h_e = 0.5 * (1.0 - xtf_h2 - xtf_h1)
    xtf_hp = xtf_h_e
    xtf_he_e = (1.0 / 3.0) * (2.0 - 2.0 * xtf_he - xtf_hep)
    xtf_hepp = 1.0 - xtf_he - xtf_hep - xtf_he_e
    # Fractions of H nuclei singly ionized, of He nuclei singly and doubly.
    # Question (see the page): H2 counted once among the H nuclei.
    xhp = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
    xhep = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
    xhepp = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)
    fxion[1] = xhp
    fxion[2] = xhep
    fxion[3] = xhepp
    return dl, d, qdt, qdp, qcp, dela, lcalc
