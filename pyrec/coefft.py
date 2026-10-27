"""Python version of coefft.f (YREC 09c3b9b): the structure equations and
their derivatives at every point, for the Henyey solution."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import dimensions_mod
from pyrec import dpmix_mod
from pyrec import fluxes_mod
from pyrec import grab_mod
from pyrec import masschg2_mod
from pyrec import nmlpar_mod
from pyrec import pualpha_mod
from pyrec import pulse1_mod
from pyrec import pulse2_mod  # USEd by the Fortran; nothing read here
from pyrec import rotder_mod
from pyrec import roten_mod
from pyrec import scrtch_mod
from pyrec import sound_mod
from pyrec.engeb import engeb
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.reduce import reduce
from pyrec.tpgrad import tpgrad

JSON = dimensions_mod.json

# COEFFT has SAVE. The nuclear terms are set only at points above TCUT(1),
# but read at every point (SEG with HT = TCUT(1) exactly; SQET, PQED, PQET
# in the cooler outer points), where they keep the values of the last point
# that set them, in this call or an earlier one.
eg = np.zeros(7)
qed = 0.0
qet = 0.0
# ENGEB's per-point rates (local arrays in the Fortran, filled at IU).
hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12, hr13, \
    hf1, hf2 = (np.zeros(JSON + 1) for _ in range(15))


def coefft(delts, m, hd, hha, hhb, hhc, hl, hmax, hp, hpp, hr, hs, hs1,
           hs2, ht, htt, hcomp, lc, tlumx, latmo, lderiv, lmix, locond,
           qdt, qdp, ksaha, fp, ft, hkerot, hkeroto, jenv, teffl, cmixl):
    """For each point: the equation of state, opacity, temperature gradient
    and nuclear energy; the structure equations and their derivatives,
    passed to REDUCE (the central conditions at point 1); the luminosity by
    source; and per-point physics stored for output and other routines.

    Translates: coefft.f#COEFFT (09c3b9b)

    Arguments: delts, the time step (s); m, the number of points; the
    model hd, hl, hp, hr, hs, ht (log rho, L/L_sun, log P, log r, log m,
    log T), hs1, hs2 (the mass at each point and in each shell), htt, hpp
    (the changes of log T, log P over the step), hcomp[species, point];
    the matrices hha[k, l, point], hhb[k, point], hhc (the gravothermal
    energy), hmax; lc (convective or not), tlumx[1..8] (the luminosity by
    source); fp, ft, the rotational factors; hkerot, hkeroto, the
    rotational kinetic energies; jenv; teffl; cmixl, the mixing length.
    Arrays have index 0 unused and are changed in place.

    Returns (latmo, lderiv, lmix, locond, qdt, qdp, ksaha): the flags set
    here, rho_T and rho_P at the last point, and the EOS's KSAHA.
    """
    global qed, qet
    c = constants_mod
    nml = nmlpar_mod
    s = scrtch_mod
    fxion = np.zeros(4)
    # The previous point's equations, for REDUCE (set at point 1 first).
    qp0 = qpr0 = qpp0 = qt0 = qtr0 = qtl0 = qtp0 = qtt0 = 0.0
    qr0 = qrr0 = qrp0 = qrt0 = ql0 = qlp0 = qlt0 = 0.0

    if nml.lsnu:
        for j in range(1, 11):
            fluxes_mod.fluxtot[j] = 0.0
    locond = True
    lderiv = True
    latmo = False
    lmix = False
    lntrpy = delts > 0.0
    if lntrpy:
        deltsi = 1.0 / delts
        # Question (see the page): a one-year floor, hard-coded.
        time = 3.1558e7
        timei = 3.1688e-8
    for j in range(1, 9):
        tlumx[j] = 0.0

    for im in range(1, m + 1):  # label 30
        # Local variables for the physics routines.
        thl = 0.0
        sl = hs[im]
        tl = ht[im]
        pl = hp[im]
        rl = hr[im]
        b = hl[im]
        x = hcomp[1, im]
        y = hcomp[2, im]
        z = hcomp[3, im]
        xhe3 = hcomp[4, im]
        xc12 = hcomp[5, im]
        xc13 = hcomp[6, im]
        xn14 = hcomp[7, im]
        xo16 = hcomp[9, im]
        xo18 = hcomp[11, im]
        xh2 = hcomp[12, im]
        iu = im
        dl = hd[im]
        fpl = fp[im]
        ftl = ft[im]
        # Debye-Hueckel composition.
        if nml.ldh:
            debhu_mod.xxdh = hcomp[1, im]
            debhu_mod.yydh = hcomp[2, im] + hcomp[4, im]
            debhu_mod.zzdh = hcomp[3, im]
            debhu_mod.zdh[1] = hcomp[5, im] + hcomp[6, im]
            debhu_mod.zdh[2] = hcomp[7, im] + hcomp[8, im]
            debhu_mod.zdh[3] = hcomp[9, im] + hcomp[10, im] + hcomp[11, im]
        (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
         qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
            tl, 0.0, pl, 0.0, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
            0.0, 0.0, qdt, qdp, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            lderiv, latmo, ksaha)
        o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
        dpmix_mod.iovim = im
        (dele, delr, qact, qacp, qacr, vel, lconv) = tpgrad(
            tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, 0.0,
            0.0, dela, qdtt, qdtp, qat, qap, 0.0, 0.0, 0.0, qcpt, qcpp,
            0.0, lderiv, False, fpl, ftl, teffl, cmixl)
        hd[im] = dl

        # The structure equations (d/d ln m) and their derivatives.
        qtemp = c.c4pil + rl + rl + rl
        qr = math.exp(c.cln * (sl - dl - qtemp))
        qrr = -qr - qr - qr
        qrp = -qr * qdp
        qrt = -qr * qdt
        # Question (see the page): f_P only, f_T unused, and constant.
        qp = -math.exp(c.cln * (c.cgl + sl + sl - pl - qtemp - rl)) * fp[im]
        qpr = -qp - qp - qp - qp
        qpp = -qp
        lc[im] = lconv
        qt = dele * qp
        qtr = -qt - qt - qt - qt
        if not lconv:
            qtl = c.clni * qt / b
            qtp = qt * qod * qdp
            qtt = qt * (-4.0 + qot + qod * qdt)
        else:
            qtl = 0.0
            qtp = qt * (-1.0 + qacp)
            qtt = qt * qact
            qtr = qtr + qt * qacr
        ql = 0.0
        qlt = 0.0
        qlp = 0.0

        # Nuclear energy above TCUT(1).
        if tl > nml.tcut[1]:
            (egg11, egg22, egg33, egg44, egg55, qed, qet, sum1, be7mfloc,
             ealpca, enu) = engeb(
                0.0, 0.0, 0.0, 0.0, 0.0, qed, qet, 0.0, dl, tl, x, y, xhe3,
                xc12, xc13, xn14, xo16, xo18, xh2, iu, hr1, hr2, hr3, hr4,
                hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12, hr13, hf1, hf2,
                0.0, 0.0, 0.0)
            e = sum1
            eg[1] = egg11
            eg[2] = egg22
            eg[3] = egg33
            eg[4] = egg44
            eg[5] = egg55
            eg[6] = enu
            eca = ealpca
            if nml.lsnu:
                for j in range(1, 11):
                    fluxes_mod.fluxtot[j] += fluxes_mod.flux[j] * hs2[im]
            for j in range(1, 7):  # label 20
                tlumx[j] = tlumx[j] + (hs2[im] / c.clsun) * eg[j]
                thl = thl + (hs2[im] / c.clsun) * eg[j]
            grab_mod.umhe33[im] = (hs2[im] / c.clsun) * grab_mod.enghe3
            grab_mod.umheal[im] = (hs2[im] / c.clsun) * grab_mod.he3all
            tlumx[8] = tlumx[8] + (hs2[im] / c.clsun) * eca
            thl = thl + (hs2[im] / c.clsun) * eca
            ql = e
            qlt = qlt + qet + qed * qdt
            qlp = qlp + qed * qdp

        # Gravothermal energy.
        if lntrpy:
            dt = deltsi
            # Accretion: points at and above JENV add DLOGT, DLOGP.
            # Question (see the page): accretion entropy.
            if nml.lmdot and nml.dmdt0 > 0.0:
                if im >= jenv:
                    httt = htt[im] + masschg2_mod.dlogt
                    hppp = hpp[im] + masschg2_mod.dlogp
                else:
                    httt = htt[im]
                    hppp = hpp[im]
            else:
                httt = htt[im]
                hppp = hpp[im]
            if hcomp[1, im] > 0.01 and delts < time:
                dt = timei
            entr1 = p * qdt / d
            entr2 = entr1 / dela
            entr = (entr2 * httt - entr1 * hppp) * c.cln
            entr3 = entr2 * c.cln * httt
            egrav = dt * entr
            hhc[im] = egrav
            tlumx[7] = tlumx[7] + (hs2[im] / c.clsun) * egrav
            ql = ql + egrav
            qlp = qlp + dt * (entr * (1.0 - qdp + qdtp) - entr1 - entr3 * qap)
            qlt = qlt + dt * (entr * (-qdt + qdtt) + entr2 - entr3 * qat)
            if nml.lrot:
                roten_mod.derot[im] = (dt * (hkerot[im] - hkeroto[im])
                                       / hs2[im])
                ql = ql - roten_mod.derot[im]
        cccql = c.clnsun * hs1[im]
        ql = cccql * ql
        qlp = cccql * qlp
        qlt = cccql * qlt

        # REDUCE for each pair of points; the central conditions at 1.
        if im > 1:
            im1 = im
            reduce(im1, hha, hhb, hl, hmax, hp, hr, hs, ht, qp0, qp, qpr0,
                   qpr, qpp0, qpp, qt0, qt, qtr0, qtr, qtl0, qtl, qtp0, qtp,
                   qtt0, qtt, qr0, qr, qrr0, qrr, qrp0, qrp, qrt0, qrt, ql0,
                   ql, qlp0, qlp, qlt0, qlt)
        else:
            hha[3, 1, 1] = c.cc13 * qdp
            hha[3, 2, 1] = c.cc13 * qdt
            hhb[3, 1] = -c.cc13 * (c.c4pi3l + dl - sl) - rl
            hha[4, 1, 1] = -qlp
            hha[4, 2, 1] = -qlt
            hhb[4, 1] = c.clni * ql - b
        qp0 = qp
        qpr0 = qpr
        qpp0 = qpp
        qt0 = qt
        qtr0 = qtr
        qtl0 = qtl
        qtp0 = qtp
        qtt0 = qtt
        qr0 = qr
        qrr0 = qrr
        qrp0 = qrp
        qrt0 = qrt
        ql0 = ql
        qlp0 = qlp
        qlt0 = qlt

        # Per-point physics for output and other routines.
        s.svel[im] = vel
        if ht[im] < nml.tcut[1]:
            s.sesum[im] = 0.0
            s.seg[7, im] = hhc[im]
            for j in range(1, 7):
                s.seg[j, im] = 0.0
        else:
            s.sesum[im] = eg[1] + eg[2] + eg[3] + eg[4] + eg[5]
            s.seg[6, im] = eg[6]
            s.seg[7, im] = hhc[im]
            if s.sesum[im] > 1.0e-22:
                esumi = 1.0 / s.sesum[im]
            else:
                esumi = 0.0
            for j in range(1, 6):
                s.seg[j, im] = eg[j] * esumi
        s.sbeta[im] = beta
        s.seta[im] = eta
        s.locons[im] = locond
        s.so[im] = o
        s.sdel[1, im] = delr
        s.sdel[2, im] = dele
        s.sdel[3, im] = dela
        for j in range(1, 4):
            s.sfxion[j, im] = fxion[j]
        s.svel[im] = vel
        s.scp[im] = qcp
        # Gamma_1, always (JVS 01/11).
        chrh = 1.0 / qdp
        cht = -chrh * qdt
        cv = qcp - math.exp(c.cln * (hp[im] - hd[im] - ht[im])) * cht**2 / chrh
        sound_mod.gam1[im] = chrh * qcp / cv
        pulse1_mod.pqdp[im] = qdp
        pulse1_mod.pqdt[im] = qdt
        if nml.lrot:
            rotder_mod.sqot[im] = qot
            if s.sesum[im] > 0.0:
                etot = s.sesum[im]
                egneut = s.seg[6, im] + s.seg[7, im]
                rotder_mod.fracnuc[im] = (etot - egneut) / etot
            else:
                rotder_mod.fracnuc[im] = 0.0
            rotder_mod.sqet[im] = qet
        if nml.lpulse:
            p1 = pulse1_mod
            pa = pualpha_mod
            p1.pqdp[im] = qdp
            p1.pqed[im] = qed
            p1.pqet[im] = qet
            p1.pqod[im] = qod
            p1.pqot[im] = qot
            p1.pqcp[im] = qcp
            p1.prmu[im] = rmu
            p1.pemu[im] = emu
            p1.pqdt[im] = qdt
            pa.valfmlt[im] = pa.alfmlt
            pa.vphmlt[im] = pa.phmlt
            pa.vcmxmlt[im] = pa.cmxmlt
    return latmo, lderiv, lmix, locond, qdt, qdp, ksaha
