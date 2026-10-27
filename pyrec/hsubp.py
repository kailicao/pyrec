"""Python version of hsubp.f (YREC 09c3b9b): the pressure scale height at
the edge of a convective region, for the overshoot."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import nmlpar_mod
from pyrec.eqstat import eqstat

# SAVE: the ionization fractions from the equation of state.
fxion = np.zeros(3 + 1)


def hsubp(hcomp, hd, hp, hr, hs, ht, jedge, pscahe):
    """Return the hydrostatic pressure scale height P r**2 / (rho G m) at
    point jedge, with rho from the equation of state at the point's P,
    T, X, Z (setting the Debye-Hueckel composition first, with
    nmlpar_mod.ldh).

    Translates: hsubp.f#HSUBP (09c3b9b)

    Arguments: hcomp, hd, hp, hr, hs, ht (the model: composition, log
    rho, log P, log r, log m, log T), jedge (the point); pscahe is only
    set.

    Returns pscahe, in cm.
    """
    c = constants_mod
    lderiv = False
    latmo = True
    ksaha = 0
    x = hcomp[1, jedge]
    z = hcomp[3, jedge]
    pl = hp[jedge]
    tl = ht[jedge]
    dl = hd[jedge]
    if nmlpar_mod.ldh:
        # Question (see the page): overwrites the Debye-Hueckel
        # composition other EOS callers may rely on.
        debhu_mod.xxdh = hcomp[1, jedge]
        debhu_mod.yydh = hcomp[2, jedge] + hcomp[4, jedge]
        debhu_mod.zzdh = hcomp[3, jedge]
        debhu_mod.zdh[1] = hcomp[5, jedge] + hcomp[6, jedge]
        debhu_mod.zdh[2] = hcomp[7, jedge] + hcomp[8, jedge]
        debhu_mod.zdh[3] = hcomp[9, jedge] + hcomp[10, jedge] + hcomp[11, jedge]
    (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp, qcp,
     dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
        tl, 0.0, pl, 0.0, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        lderiv, latmo, ksaha)
    sl = hs[jedge]
    rl = hr[jedge]
    # Question (see the page): the ideal-gas estimate pscap is not used.
    pscap = rmu * math.exp(c.cln * (tl - c.cgl - sl + rl + rl))
    pscahe = math.exp(c.cln * (pl + 2.0 * rl - dl - c.cgl - sl))
    return pscahe
