"""Python version of physic.f (YREC development, 09c3b9b): the auxiliary
physics at every shell, and the shear-stability limit on dOmega/dln r."""

import math

import numpy as np

from pyrec import comp2_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import dpmix_mod
from pyrec import dwmax_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod
from pyrec import temp_mod
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.intrp2 import intrp2
from pyrec.tpgrad import tpgrad
from pyrec.viscos import viscos

# DATA ATOMWT: atomic weights of H, He4, the metals (as C12) and He3.
atomwt = np.array([0.0, 1.007825, 4.002603, 12.0, 3.01603])

# PHYSIC has SAVE: its locals are kept between calls. KSAHA, which EQSTAT
# updates, is never initialized here (Question, see the page).
ksaha = 0
fxion = np.zeros(4)
hf = np.zeros(5)
hold = np.zeros(5)


def physic(fp, ft, hcomp, hd, hg, hl, hp, hr, hs, ht, lc, m, teffl, cmixl):
    """Recompute, at every shell, the equation of state, opacity,
    temperature gradients and convective velocity, and the mean molecular
    weight; then the thermal diffusivity and viscosity (VISCOS), and the
    largest dOmega/dln r stable against the shear instability.

    Translates: physic.f#PHYSIC (09c3b9b)

    Arguments: fp, ft, the rotational distortion factors; the model hcomp
    [species, shell], hd, hg, hl, hp, hr, hs, ht (arrays, index 0 unused);
    lc, the convection flags (changed in place); m, the number of shells;
    teffl, log Teff; cmixl, the mixing length.

    Returns nothing; sets lc and, in the modules, SDEL, SO, SVEL
    (scrtch_mod), CP, HAMU, SQDT (temp_mod; THDIF and VISC through VISCOS),
    QWRMAX (dwmax_mod) and IOVIM (dpmix_mod).
    """
    global ksaha
    c = constants_mod
    sdel = scrtch_mod.sdel
    # Actual and adiabatic temperature gradients, opacity and mean
    # molecular weight for all shells (the Fortran comment says "radiative
    # shells": Question, see the page).
    lderiv = False
    locond = False  # Question (see the page): unused.
    latmo = False
    for im in range(1, m + 1):
        sl = hs[im]
        tl = ht[im]
        pl = hp[im]
        rl = hr[im]
        b = hl[im]
        x = hcomp[1, im]
        z = hcomp[3, im]
        dl = hd[im]
        fpl = fp[im]
        ftl = ft[im]
        # The Debye-Hueckel correction, if chosen; the MHD equation of state is
        # removed (9/22/26 MHP via Claude).
        if nmlpar_mod.ldh:
            debhu_mod.xxdh = hcomp[1, im]
            debhu_mod.yydh = hcomp[2, im] + hcomp[4, im]
            debhu_mod.zzdh = hcomp[3, im]
            debhu_mod.zdh[1] = hcomp[5, im] + hcomp[6, im]
            debhu_mod.zdh[2] = hcomp[7, im] + hcomp[8, im]
            debhu_mod.zdh[3] = hcomp[9, im] + hcomp[10, im] + hcomp[11, im]
        (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
         qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
            tl, 0.0, pl, 0.0, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            lderiv, latmo, ksaha)
        o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
        dpmix_mod.iovim = im
        (dele, delr, qact, qacp, qacr, vel, lconv) = tpgrad(
            tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, 0.0,
            0.0, dela, qdtt, qdtp, qat, qap, 0.0, 0.0, 0.0, qcpt, qcpp,
            0.0, lderiv, False, fpl, ftl, teffl, cmixl)
        # Question (see the page): the caller's LC is overwritten.
        lc[im] = lconv
        sdel[1, im] = delr
        sdel[2, im] = dele
        sdel[3, im] = dela
        # Mean molecular weight of a fully ionized gas. AMUENV is 1/(mean
        # molecular weight per ion) of the surface mixture.
        dfx1 = hcomp[1, im] - comp_mod.xenv
        dfx2 = hcomp[2, im] - comp2_mod.yenv
        dfx3 = hcomp[3, im] - comp_mod.zenv
        dfx4 = hcomp[4, im] - comp2_mod.y3env
        temp = (comp_mod.amuenv + dfx1 / atomwt[1] + dfx2 / atomwt[2]
                + dfx3 / atomwt[3] + dfx4 / atomwt[4])
        amu2 = 1.0 / temp
        temp = (hcomp[1, im] / atomwt[1]
                + 2.0 * (hcomp[4, im] / atomwt[4] + hcomp[2, im] / atomwt[2])
                + 0.5 * hcomp[3, im])
        emu2 = 1.0 / temp
        temp_mod.hamu[im] = amu2 * emu2 / (amu2 + emu2)
        scrtch_mod.so[im] = o
        temp_mod.cp[im] = qcp
        temp_mod.sqdt[im] = qdt
        # SVEL (the convective velocity) is always kept (JVS 10/13).
        scrtch_mod.svel[im] = vel
    # label 30
    # Thermometric diffusivity and kinematic viscosity.
    viscos(hcomp, hd, ht, m)
    # The largest stable dOmega/dln r at each interface (Endal & Sofia
    # 1978), from quantities interpolated to the midpoint in HS.
    qwrmax = dwmax_mod.qwrmax
    for im in range(2, m + 1):
        # Skip convective regions.
        if lc[im] and lc[im - 1]:
            qwrmax[im] = 0.0
            continue  # label 100
        if im <= 2:
            k = 1
        elif im == m:
            k = m - 3
        else:
            k = im - 2
        for i in range(1, 5):
            hold[i] = hs[i + k - 1]
        hsmid = 0.5 * (hs[im] + hs[im - 1])
        intrp2(hold, hf, hsmid)
        pmid = (hp[k] * hf[1] + hp[k + 1] * hf[2] + hp[k + 2] * hf[3]
                + hp[k + 3] * hf[4])
        dmid = (hd[k] * hf[1] + hd[k + 1] * hf[2] + hd[k + 2] * hf[3]
                + hd[k + 3] * hf[4])
        delmid = (sdel[2, k] * hf[1] + sdel[2, k + 1] * hf[2]
                  + sdel[2, k + 2] * hf[3] + sdel[2, k + 3] * hf[4])
        delamd = (sdel[3, k] * hf[1] + sdel[3, k + 1] * hf[2]
                  + sdel[3, k + 2] * hf[3] + sdel[3, k + 3] * hf[4])
        gmid = (hg[k] * hf[1] + hg[k + 1] * hf[2] + hg[k + 2] * hf[3]
                + hg[k + 3] * hf[4])
        temp = math.exp(c.cln * (dmid - pmid)) * (delamd - delmid) * gmid**2
        if temp > 0.0:
            qwrmax[im] = 2.0 * math.sqrt(temp)
        else:
            qwrmax[im] = 0.0
