"""Python version of surfp.f (YREC 09c3b9b): the pressure at T = Teff from
the Kurucz model-atmosphere table (KTTAU = 3), the outer boundary of the
envelope integration."""

import numpy as np

from pyrec import atmos2_mod
from pyrec import atmprt_mod
from pyrec import fac_mod
from pyrec import iounits_mod
from pyrec.inter3 import inter3
from pyrec.kspline import kspline
from pyrec.ksplint import ksplint

# PARAMETER (NT=57, NG=11): the table's temperatures and gravities.
NT = 57
NG = 11


def _spline4(qg, ptab, x):
    """KSPLINE then KSPLINT on four points: the spline's value at x."""
    yg = np.zeros(5)
    kspline(qg, ptab, yg)
    return ksplint(qg, ptab, yg, x, 0.0)


def surfp(teffl, gl, lprt):
    """log P at T = Teff, interpolated in log g at four table temperatures
    around log Teff (linearly, by 3-point Lagrange or by a 4-point natural
    spline, by temperature), then in log Teff by a 4-point natural spline.

    Translates: surfp.f#SURFP (09c3b9b)

    Arguments: teffl (log Teff), gl (log g); lprt, to print the result on
    ISHORT and ISTOR. Sets atmprt_mod.ap (log P) and atmprt_mod.at
    (log Teff), and fac_mod.jp and fac_mod.kp[1..4] (the table cell used).
    None is returned.
    """
    qt = np.zeros(5)
    pp = np.zeros(5)
    qs = np.zeros(4)
    qgg = np.zeros(4)
    dum = np.zeros(4)
    qg = np.zeros(5)
    ptab = np.zeros(5)
    atmpl = atmos2_mod.atmpl
    atmtl = atmos2_mod.atmtl
    atmgl = atmos2_mod.atmgl

    # Question (see the page): only the lower limits are checked, and they
    # are hard-coded.
    if teffl < 3.5 or gl < -0.5:
        for unit in (iounits_mod.iowr, iounits_mod.ishort):
            iounits_mod.write(unit, " DESIRED ATMOSPHERE OUTSIDE TABLE RANGE")
            iounits_mod.write(unit, f" LOG TEFF{teffl:10.6f} LOG G{gl:10.6f}")
            iounits_mod.write(unit, " RUN STOPPED")  # format 911
        raise SystemExit()

    # The four table temperatures around teffl.
    for j in range(1, NT + 1):
        if teffl <= atmtl[j]:
            break  # to label 10
    else:
        j = NT
    jj = max(1, j - 2)  # label 10
    jj = min(NT - 3, jj)
    for k in range(1, 5):
        qt[k] = atmtl[jj + k - 1]
    fac_mod.jp = jj

    # log P in log g at each of the four temperatures (loop 20).
    for j in range(jj, jj + 4):
        n = j - jj + 1
        if atmtl[j] > 4.5:
            # Question (see the page): hot temperatures use only the
            # highest gravities, without checking IMINMAX(J).
            if atmtl[j] > 4.55:
                fx = (gl - atmgl[NG - 1]) / (atmgl[NG] - atmgl[NG - 1])
                pp[n] = atmpl[j, NG - 1] + fx * (atmpl[j, NG]
                                                 - atmpl[j, NG - 1])
                fac_mod.kp[n] = NG - 1
            else:
                for k in range(1, 4):
                    qs[k] = atmgl[NG - 3 + k]
                inter3(qs, qgg, dum, gl)
                pp[n] = (atmpl[j, NG - 2] * qgg[1]
                         + atmpl[j, NG - 1] * qgg[2]
                         + atmpl[j, NG] * qgg[3])
                fac_mod.kp[n] = NG - 2
            continue

        if gl >= atmgl[NG - 1]:
            # The top four gravities.
            # Question (see the page): IMINMAX(J) is not checked here.
            for kk in range(1, 5):
                qg[kk] = atmgl[NG - 4 + kk]
                ptab[kk] = atmpl[j, NG - 4 + kk]
            pp[n] = _spline4(qg, ptab, gl)
            fac_mod.kp[n] = NG - 3
            continue

        # General case: the four valid gravities bracketing gl.
        found = False
        for k in range(fac_mod.iminmax[j] - 3, fac_mod.imin[j] - 1, -1):
            if gl < atmgl[k + 2] and gl >= atmgl[k + 1]:
                kk = max(fac_mod.imin[j], k)
                kk = min(NG - 3, kk)
                for kkk in range(1, 5):
                    qg[kkk] = atmgl[kk + kkk - 1]
                    ptab[kkk] = atmpl[j, kk + kkk - 1]
                pp[n] = _spline4(qg, ptab, gl)
                fac_mod.kp[n] = kk
                found = True
                break
        if found:
            continue

        # The lowest four valid gravities.
        for k in range(1, 5):
            qg[k] = atmgl[k + fac_mod.imin[j] - 1]
            ptab[k] = atmpl[j, k + fac_mod.imin[j] - 1]
        pp[n] = _spline4(qg, ptab, gl)
        fac_mod.kp[n] = fac_mod.imin[j]

    # log P in log Teff.
    atmprt_mod.ap = _spline4(qt, pp, teffl)
    atmprt_mod.at = teffl

    if lprt:
        for unit in (iounits_mod.ishort, iounits_mod.istor):
            iounits_mod.write(unit, "********PRESSURE AT T=TEFF INTERPOLATED "
                              "FROM TABULATED VALUES********")  # format 70
        for unit in (iounits_mod.ishort, iounits_mod.istor):
            iounits_mod.write(unit, f" {'':20s}LOG (Teff) ={teffl:10.5f}"
                              f" LOG P ={atmprt_mod.ap:10.5f}")  # format 71
