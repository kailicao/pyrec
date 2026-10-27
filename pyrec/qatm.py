"""Python version of qatm.f (YREC 09c3b9b): the right-hand side of the
atmosphere integration, d log P / d log tau."""

import math

import numpy as np

from pyrec import atmos_mod
from pyrec import atmprt_mod
from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import pulse1_mod
from pyrec import pulse2_mod
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.hra import hra

# Locals kept between calls (SAVE): TL keeps its value when KTTAU is not
# 0, 1 or 2.
tl = 0.0


def qatm(x0, y, dydx, b, fpl, ftl, gl, latmo, lderiv, locond, lsave,
         rl, teffl, x, z, katm, ksaha):
    """Derivative of log P with respect to log tau in the atmosphere, from
    a T-tau relation (Eddington, Krishna Swamy or the Harvard reference
    atmosphere, by KTTAU) and hydrostatic equilibrium dP/dtau = g/kappa.
    Also stores the point for output in atmprt_mod and pulse2_mod.

    Translates: qatm.f#QATM (09c3b9b)

    Arguments: x0 = log tau; y (array of 3, index 0 unused), y[1] = log P;
    the surface gravity gl, distortion factor fpl, teffl, x, z; lsave
    (store the point for output). b, ftl, locond and rl are not used: the
    argument list is that of the derivative routine called by bsstep and
    mmid. Arrays changed in place: dydx (dydx[1] set).

    Returns (katm, ksaha): the atmosphere-call counter, incremented, and
    the counter of Saha calls, from eqstat.
    """
    global tl
    fxion = np.zeros(4)   # index 0 unused

    # Statement functions: Eddington, and Krishna Swamy (1966, ApJ 145, 176).
    def ttaul0(yy):
        return teffl - 0.031235 + 0.25 * math.log10(yy + constants_mod.cc23)

    def ttaul1(yy):
        return teffl - 0.031235 + 0.25 * math.log10(
            yy + 1.39 - 0.815 * math.exp(-2.54 * yy)
            - 0.025 * math.exp(-30.0 * yy))

    g = math.exp(constants_mod.cln * gl) * fpl
    atmprt_mod.taul = x0
    tau = math.exp(constants_mod.cln * atmprt_mod.taul)
    if nmlpar_mod.kttau == 0:
        tl = ttaul0(tau)
    elif nmlpar_mod.kttau == 1:
        tl = ttaul1(tau)
    elif nmlpar_mod.kttau == 2:
        tl = teffl + hra(tau) - atmos_mod.hras
    # Question (see the page): for KTTAU 3-5 TL keeps its last value.
    pl = y[1]
    (t, p, dl, d, beta, beta1, beta14, rmu, amu, emu, eta, qdt, qdp, qcp,
     dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
        tl, 0.0, pl, 0.0, 0.0, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        lderiv, latmo, ksaha)
    o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
    dydx[1] = g * tau / (p * o)
    katm = katm + 1
    atmprt_mod.ap = pl
    atmprt_mod.at = tl
    if lsave or pulse1_mod.lpumod:
        atmprt_mod.ad = dl
        atmprt_mod.ao = o
        atmprt_mod.afxion[1] = fxion[1]
        atmprt_mod.afxion[2] = fxion[2]
        atmprt_mod.afxion[3] = fxion[3]
        q = pulse2_mod
        q.qtl = tl
        q.qt = math.exp(constants_mod.cln * tl)
        q.qpl = pl
        q.qp = math.exp(constants_mod.cln * pl)
        q.qdl = dl
        q.qd = math.exp(constants_mod.cln * dl)
        q.qo = o
        q.qol = ol
        q.qqdp = qdp
        q.qqed = 0.0
        q.qqod = qod
        q.qqot = qot
        q.qdel = 0.0
        q.qqdt = qdt
        q.qdela = dela
        q.qqcp = qcp
        q.qrmu = rmu
        q.qemu = emu
    return katm, ksaha
