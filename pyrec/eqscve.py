"""Python version of eqscve.f (YREC development, 09c3b9b): the SCV
equation of state at the envelope composition."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import scveos_mod
from pyrec.eqscvg import eqscvg
from pyrec.interp import interp
from pyrec.kspline import kspline
from pyrec.ksplint import ksplint

NTS = 63                  # temperatures in the SCV tables
TOLSP = 0.08              # smoothing zone in log P (DATA)
TOLST = 0.032             # smoothing zone in log T (DATA)


def _stencil_spline(ismt, ismp, pp, tl):
    """The interpolation block the Fortran writes out four times: for the
    4x4 stencil shifted by ismt rows in T and ismp columns in P, a natural
    spline in log Pgas at each of the 4 temperatures, then one in log T,
    for columns 2-6 of TABLENV. The pressure abscissae always come from
    the row idtscv (unshifted), as in the Fortran. Returns tempt[1..5]."""
    scv = scveos_mod
    idt = scv.idtscv
    idp = scv.idp
    qr = np.zeros(5)
    ytab = np.zeros(5)
    y2 = np.zeros(5)
    temp = np.zeros((5, 6))
    tempt = np.zeros(6)
    # Interpolate in pressure at 4 different temperature points.
    for i in range(1, 5):
        ii = idt + i - 1 + ismt
        for k in range(1, 5):
            # Question (see the page): row idt's pressures for all 4 rows.
            qr[k] = scv.tablenv[idt, idp + k - 1 + ismp, 1]
        for j in range(1, 6):
            for k in range(1, 5):
                ytab[k] = scv.tablenv[ii, idp + k - 1 + ismp, j + 1]
            kspline(qr, ytab, y2)
            temp[i, j] = ksplint(qr, ytab, y2, pp, 0.0)
    # Interpolate in temperature.
    for k in range(1, 5):
        qr[k] = scv.tlogx[k + idt - 1 + ismt]
    for j in range(1, 6):
        for k in range(1, 5):
            ytab[k] = temp[k, j]
        kspline(qr, ytab, y2)
        tempt[j] = ksplint(qr, ytab, y2, tl, 0.0)
    return tempt


def eqscve(tl, t, p, dl, d, x, z, beta, fxion, qdt, qdp, qcp, dela, lcalc):
    """Interpolate the Saumon-Chabrier-Van Horn equation of state in its
    table for the envelope composition (2D splines, smoothed where the
    stencil changes), add radiation, and get the H and He ionization
    fractions; any other composition goes to eqscvg.

    Translates: eqscve.f#EQSCVE (09c3b9b)

    Arguments: tl, t (log T, T), p (total pressure), x, z, beta (Pgas/P);
    dl ... dela and lcalc are only set. Array changed in place: fxion (H+,
    He+, He++ fractions, index 0 unused).

    Returns (dl, d, qdt, qdp, qcp, dela, lcalc). Outside the table, lcalc
    is False and the others are returned as passed in.
    """
    c = constants_mod
    scv = scveos_mod
    if (abs(x - comp_mod.xenv) > 1.0e-5
            or abs(z - comp_mod.zenv) > 1.0e-5):
        dl, d, qdt, qdp, qcp, dela, lcalc = eqscvg(
            tl, t, p, dl, d, x, z, beta, fxion, qdt, qdp, qcp, dela, lcalc)
        return dl, d, qdt, qdp, qcp, dela, lcalc
    # The gas pressure is the tabulated quantity.
    betm1 = 1.0 - beta
    pp = math.log10(beta * p)
    # Is the point within the table?
    if tl < scv.tlogx[1] or tl > scv.tlogx[NTS] or pp < 4.0:
        lcalc = False
        return dl, d, qdt, qdp, qcp, dela, lcalc

    # Nearest 4 points in temperature, searching from the last stencil.
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
    # Smoothing in T near the edges of the middle interval.
    lsmt = False
    ismt = 0
    fst = 0.0
    if idt == ii:
        dtabove = tl - scv.tlogx[idt + 1]
        dtbelow = scv.tlogx[idt + 2] - tl
        if dtabove < TOLST:
            if idt > 1:
                lsmt = True
                ismt = -1
                fst = (dtabove + TOLST) / (2.0 * TOLST)
        elif dtbelow < TOLST:
            if idt < NTS - 3:
                lsmt = True
                ismt = 1
                fst = 0.5 * dtbelow / TOLST

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
    # Smoothing in P near the edges of the middle interval.
    lsmp = False
    ismp = 0
    fsp = 0.0
    if idp == jj:
        dpabove = pp - scv.tablenv[idt, idp + 1, 1]
        dpbelow = scv.tablenv[idt, idp + 2, 1] - pp
        if dpabove < TOLSP:
            if idp > 1:
                lsmp = True
                ismp = -1
                fsp = (dpabove + TOLSP) / (2.0 * TOLSP)
        elif dpbelow < TOLSP:
            if idp < nptsx - 3:
                lsmp = True
                ismp = 1
                fsp = 0.5 * dpbelow / TOLSP

    y = 1.0 - x - z                                # set, never used
    # Include radiation pressure in the equation of state.
    prad = betm1 * p
    pgas = beta * p                                # set, never used
    # Columns 2-6 of TABLENV: log rho, Q_rhoT, Q_rhoP, cP, dU/dT.
    tempt = _stencil_spline(0, 0, pp, tl)
    tempt1t = np.zeros(6)
    tempt2p = np.zeros(6)
    if lsmt:
        # Interpolate at the stencil shifted in T.
        tempt1 = _stencil_spline(ismt, 0, pp, tl)
        if ismt == -1:
            tempt1t[1:] = tempt1[1:] + fst * (tempt[1:] - tempt1[1:])
        else:
            tempt1t[1:] = tempt[1:] + fst * (tempt1[1:] - tempt[1:])
    if lsmp:
        # Interpolate at the stencil shifted in P.
        tempt2 = _stencil_spline(0, ismp, pp, tl)
        if ismp == -1:
            tempt2p[1:] = tempt2[1:] + fsp * (tempt[1:] - tempt2[1:])
        else:
            tempt2p[1:] = tempt[1:] + fsp * (tempt2[1:] - tempt[1:])
    if lsmp and lsmt:
        # Interpolate at the stencil shifted in both.
        tempt3 = _stencil_spline(ismt, ismp, pp, tl)
    if lsmt:
        if lsmp:
            # Blend in T at both pressures, then in P.
            for j in range(1, 6):
                if ismt == -1:
                    chgt1 = tempt1[j] + fst * (tempt[j] - tempt1[j])
                    chgt2 = tempt3[j] + fst * (tempt2[j] - tempt3[j])
                else:
                    chgt1 = tempt[j] + fst * (tempt1[j] - tempt[j])
                    chgt2 = tempt2[j] + fst * (tempt3[j] - tempt2[j])
                if ismp == -1:
                    tempt[j] = chgt2 + fsp * (chgt1 - chgt2)
                else:
                    tempt[j] = chgt1 + fsp * (chgt2 - chgt1)
        else:
            tempt[1:] = tempt1t[1:]
    elif lsmp:
        tempt[1:] = tempt2p[1:]

    # Density.
    dl = tempt[1]
    d = math.exp(c.cln * dl)
    # Q_rhoP: radiation adds pressure but no density, so divide by beta.
    qdp0 = tempt[3]
    qdp = qdp0 / beta
    # Q_rhoT = -chi_T Q_rhoP, with chi_T (qpt) corrected for radiation.
    qdt0 = tempt[2]
    qpt0 = -qdt0 / qdp0
    qpt = qpt0 * beta + 4.0 * betm1
    qdt = -qpt * qdp
    # Question (see the page): the tabulated cP is not used.
    qcp0 = tempt[4]
    # dU/dT from the table, corrected for radiation.
    qut = tempt[5]
    ur = 3.0 * prad / d
    qut = qut + 4.0 * ur / t
    # cP corrected for radiation pressure, and the adiabatic gradient.
    qcp = qut - p * qdt * qpt / d / t
    dela = -p * qdt / d / t / qcp

    # Ionization: particle fractions from the pure-H (TABLEX) and pure-He
    # (TABLEY) tables, by 4-point Lagrange interpolation on the same stencil.
    qr = np.zeros(5)
    ft = np.zeros(5)
    ftd = np.zeros(5)
    fp = np.zeros(5)
    fpd = np.zeros(5)
    temp = np.zeros((5, 3))
    for k in range(1, 5):
        qr[k] = scv.tlogx[k + idt - 1]
    interp(qr, ft, ftd, tl)
    for k in range(1, 5):
        # Question (see the page): TABLENV's pressures for TABLEX, TABLEY.
        qr[k] = scv.tablenv[idt, idp + k - 1, 1]
    interp(qr, fp, fpd, pp)

    def lagrange(col):
        """Columns col of TABLEX and TABLEY at (tl, pp)."""
        for i in range(1, 5):
            ii = idt + i - 1
            temp[i, 1] = sum(fp[k] * scv.tablex[ii, idp + k - 1, col]
                             for k in range(1, 5))
            temp[i, 2] = sum(fp[k] * scv.tabley[ii, idp + k - 1, col]
                             for k in range(1, 5))
        return (sum(ft[k] * temp[k, 1] for k in range(1, 5)),
                sum(ft[k] * temp[k, 2] for k in range(1, 5)))

    # Column 2: fractions of all particles that are H2, and neutral He.
    xtf_h2, xtf_he = lagrange(2)
    # Column 3: fractions that are neutral H atoms, and He+.
    xtf_h1, xtf_hep = lagrange(3)
    # Particle and charge conservation give the rest.
    xtf_h_e = 0.5 * (1.0 - xtf_h2 - xtf_h1)
    xtf_hp = xtf_h_e
    xtf_he_e = (1.0 / 3.0) * (2.0 - 2.0 * xtf_he - xtf_hep)
    xtf_hepp = 1.0 - xtf_he - xtf_hep - xtf_he_e
    # Fractions of H and He nuclei in each stage.
    # Question (see the page): H2 counted once in the denominator.
    xhp = xtf_hp / (xtf_h2 + xtf_h1 + xtf_hp)
    xhep = xtf_hep / (xtf_he + xtf_hep + xtf_hepp)
    xhepp = xtf_hepp / (xtf_he + xtf_hep + xtf_hepp)
    fxion[1] = xhp
    fxion[2] = xhep
    fxion[3] = xhepp
    return dl, d, qdt, qdp, qcp, dela, lcalc
