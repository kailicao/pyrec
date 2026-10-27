"""Python version of qenv.f (YREC 09c3b9b): the right-hand side of the
envelope integration in log P."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import dpmix_mod
from pyrec import envprt_mod
from pyrec import pulse1_mod
from pyrec import pulse2_mod
from pyrec import qenvmix_mod
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.tpgrad import tpgrad


def qenv(x0, y, dydx, b, fpl, ftl, gl, latmo, lderiv, locond, lsave,
         rl, teffl, x, z, kenv, ksaha):
    """Derivatives of log(m/M), log T and log r with respect to log P in
    the envelope, from hydrostatic equilibrium, mass conservation and the
    temperature gradient of tpgrad. Also stores the point in envprt_mod
    (always) and pulse2_mod (when saving).

    Translates: qenv.f#QENV (09c3b9b)

    Arguments: x0 = log P; y (array of 3, index 0 unused) = log(m/M),
    log T, log r; the luminosity b, distortion factors fpl, ftl, teffl,
    x, z; lsave (store the point for output). The mixing length comes from
    qenvmix_mod.cmixlq, since the argument list is that of the derivative
    routine called by bsstep and mmid (gl and locond are not used).
    Arrays changed in place: dydx[1..3].

    Returns (rl, kenv, ksaha): rl, set to y[3]; the envelope-call
    counter, incremented; the counter of Saha calls, from eqstat.
    """
    c = constants_mod
    fxion = np.zeros(4)   # index 0 unused

    pl = x0
    sl = y[1] + comp_mod.stotal
    tl = y[2]
    rl = y[3]
    (t, p, dl, d, beta, beta1, beta14, rmu, amu, emu, eta, qdt, qdp, qcp,
     dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
        tl, 0.0, pl, 0.0, 0.0, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        lderiv, latmo, ksaha)
    o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
    dpmix_mod.iovim = -1
    (dele, delr, qact, qacp, qacr, vel, lconv) = tpgrad(
        tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, 0.0,
        0.0, dela, qdtt, qdtp, qat, qap, 0.0, 0.0, 0.0, qcpt, qcpp,
        0.0, lderiv, False, fpl, ftl, teffl, qenvmix_mod.cmixlq)
    dydx[1] = -math.exp(c.cln * (c.c4pil + 4.0 * rl + pl - c.cgl - sl - sl)) / fpl
    dydx[2] = dele
    # Question (see the page): FPL multiplies here, where COEFFT's form
    # would divide.
    dydx[3] = -math.exp(c.cln * (pl + rl - c.cgl - sl - dl)) * fpl
    kenv = kenv + 1
    # 07/02: always store the basic structure variables.
    e = envprt_mod
    e.ep = pl
    e.et = tl
    e.es = sl - comp_mod.stotal
    e.er = rl
    e.ed = dl
    e.evel = vel
    # JVS 08/13: always store the gradients (for tracking the CZ).
    e.edel[1] = delr
    e.edel[2] = dela
    e.edel[3] = dele
    e.ebeta = beta
    # JVS 08/25: also always save these, for the profile.
    e.efxion[1] = fxion[1]
    e.efxion[2] = fxion[2]
    e.efxion[3] = fxion[3]
    q = pulse2_mod
    q.qqdp = qdp
    q.qqdt = qdt
    q.qqcp = qcp

    if lsave or pulse1_mod.lpumod:
        e.eo = o
        # Question (see the page): these repeat the assignments above.
        e.efxion[1] = fxion[1]
        e.efxion[2] = fxion[2]
        e.efxion[3] = fxion[3]
        q.qtl = tl
        q.qt = math.exp(c.cln * tl)
        q.qpl = pl
        q.qp = math.exp(c.cln * pl)
        q.qdl = dl
        q.qd = math.exp(c.cln * dl)
        q.qo = o
        q.qol = ol
        q.qfs = math.exp(c.cln * (sl - comp_mod.stotal))
        q.qqdp = qdp
        q.qqed = 0.0
        q.qqod = qod
        q.qqot = qot
        q.qdel = dele
        q.qqdt = qdt
        q.qdela = dela
        q.qqcp = qcp
        q.qrmu = rmu
        q.qemu = emu
    return rl, kenv, ksaha
