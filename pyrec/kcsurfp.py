"""Python version of kcsurfp.f (YREC development, 09c3b9b)."""

import numpy as np

from pyrec import atmos2c_mod
from pyrec import atmprt_mod
from pyrec import fac_mod
from pyrec import iounits_mod
from pyrec.inter3 import inter3
from pyrec.kspline import kspline
from pyrec.ksplint import ksplint

# Parameters NTC and NGC for the tabulated surface pressures.
NTC = 76
NGC = 11


def _spline4(qg, ptab, x):
    """KSPLINE then KSPLINT on four points: the spline's value at x."""
    yg = np.zeros(5)
    kspline(qg, ptab, yg)
    return ksplint(qg, ptab, yg, x, 0.0)


def kcsurfp(teffl, gl, lprt):
    """Gas pressure at T = Teff from the Kurucz/Castelli atmosphere table.

    The version of SURFP for the newer Kurucz/Castelli tables (76 x 11):
    log P is interpolated in log g at four table temperatures around
    log Teff (linearly, by 3-point Lagrange or by a 4-point natural spline,
    depending on the temperature), then in log Teff by a 4-point natural
    spline. Called only by ENVINT.

    Translates: kcsurfp.f#KCSURFP (09c3b9b)

    Parameters
    ----------
    teffl : float
        log Teff.
    gl : float
        log g.
    lprt : bool
        Write the result to the short output and the model-storage unit.

    Sets ``atmprt_mod.ap`` (log P) and ``atmprt_mod.at`` (log Teff), and
    ``fac_mod.jp`` and ``fac_mod.kp[1..4]`` (the table corner used).
    Nothing is returned.
    """
    qt = np.zeros(5)
    pp = np.zeros(5)
    qs = np.zeros(4)
    qgg = np.zeros(4)
    dum = np.zeros(4)
    qg = np.zeros(5)
    ptab = np.zeros(5)
    atmplc = atmos2c_mod.atmplc
    atmtlc = atmos2c_mod.atmtlc
    atmglc = atmos2c_mod.atmglc

    # Check that the point is within the table.
    # Question (see the page): hard-coded bounds; no upper bounds checked.
    if teffl < 3.5 or gl < -0.5:
        for unit in (iounits_mod.iowr, iounits_mod.ishort):
            iounits_mod.write(unit, " DESIRED ATMOSPHERE OUTSIDE TABLE RANGE")
            iounits_mod.write(unit, f" LOG TEFF{teffl:10.6f} LOG G{gl:10.6f}")
            iounits_mod.write(unit, " RUN STOPPED")  # format 911
        raise SystemExit

    # Temperature interpolation: the four table points around teffl.
    for j in range(1, NTC + 1):
        if teffl <= atmtlc[j]:
            break  # to label 10
    else:
        j = NTC
    jj = max(1, j - 2)  # label 10
    jj = min(NTC - 3, jj)
    for k in range(1, 5):
        qt[k] = atmtlc[jj + k - 1]
    fac_mod.jp = jj

    # Gravity interpolation at each of the four temperatures (loop 20).
    for j in range(jj, jj + 4):
        n = j - jj + 1
        if atmtlc[j] > 4.5:
            if atmtlc[j] > 4.55:
                # Linear interpolation (extrapolation) in the two highest
                # gravities. Question (see the page): unbounded in log g.
                fx = (gl - atmglc[NGC - 1]) / (atmglc[NGC] - atmglc[NGC - 1])
                pp[n] = (atmplc[j, NGC - 1]
                         + fx * (atmplc[j, NGC] - atmplc[j, NGC - 1]))
                fac_mod.kp[n] = NGC - 1
            else:
                # 3-point Lagrange interpolation in the three highest.
                for k in range(1, 4):
                    qs[k] = atmglc[NGC - 3 + k]
                inter3(qs, qgg, dum, gl)
                pp[n] = (atmplc[j, NGC - 2] * qgg[1]
                         + atmplc[j, NGC - 1] * qgg[2]
                         + atmplc[j, NGC] * qgg[3])
                fac_mod.kp[n] = NGC - 2
            continue

        # Question (see the page): this branch and the search below use
        # NGC, not IMINMAX2(J), so missing top entries can enter.
        if gl >= atmglc[NGC - 1]:
            # Above the second-highest gravity: the top four.
            for kk in range(1, 5):
                qg[kk] = atmglc[NGC - 4 + kk]
                ptab[kk] = atmplc[j, NGC - 4 + kk]
            pp[n] = _spline4(qg, ptab, gl)
            fac_mod.kp[n] = NGC - 3
            continue

        # General case: the four nearest valid gravities.
        found = False
        for k in range(fac_mod.iminmax2[j] - 3, fac_mod.imin2[j] - 1, -1):
            if atmglc[k + 2] > gl >= atmglc[k + 1]:
                kk = max(fac_mod.imin2[j], k)
                kk = min(NGC - 3, kk)
                for kkk in range(1, 5):
                    qg[kkk] = atmglc[kk + kkk - 1]
                    ptab[kkk] = atmplc[j, kk + kkk - 1]
                pp[n] = _spline4(qg, ptab, gl)
                fac_mod.kp[n] = kk
                found = True
                break
        if found:
            continue

        # Below the second valid gravity: the lowest four valid ones.
        for k in range(1, 5):
            qg[k] = atmglc[k + fac_mod.imin2[j] - 1]
            ptab[k] = atmplc[j, k + fac_mod.imin2[j] - 1]
        pp[n] = _spline4(qg, ptab, gl)
        fac_mod.kp[n] = fac_mod.imin2[j]

    # Interpolate in temperature to find log P.
    atmprt_mod.ap = _spline4(qt, pp, teffl)
    atmprt_mod.at = teffl

    # Write the result to the model file.
    if lprt:
        for unit in (iounits_mod.ishort, iounits_mod.istor):
            iounits_mod.write(unit, "********PRESSURE AT T=TEFF INTERPOLATED "
                              "FROM TABULATED VALUES********")  # format 70
        for unit in (iounits_mod.ishort, iounits_mod.istor):
            iounits_mod.write(unit, f" {'':20s}LOG (Teff) ={teffl:10.5f}"
                              f" LOG P ={atmprt_mod.ap:10.5f}")  # format 71
