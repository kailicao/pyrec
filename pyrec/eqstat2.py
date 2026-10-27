"""Python version of eqstat2.f (YREC development, 09c3b9b): the core of
the equation of state, choosing and blending Saha, SCV, the relativistic
gas and OPAL 2006 by regime."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.eqbound06 import eqbound06
from pyrec.eqrelv import eqrelv
from pyrec.eqsaha import eqsaha
from pyrec.eqscve import eqscve
from pyrec.oeqos06 import oeqos06

# DATA, kept between calls (SAVE).
NZP1 = 12
ATOMWT = np.array([0.0, 0.9921, 0.24975, 0.08322, 0.4995])
SCALET = 0.500
SCALES = 2.000
# Atomic weights of Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He.
ATOMW2 = np.array([0.0, 23.0, 26.99, 24.32, 55.86, 28.1, 12.015,
                   1.008, 16.0, 14.01, 39.96, 20.19, 4.004])

# Locals kept between calls (SAVE).
fx = np.zeros(NZP1 + 1)
qdtt0 = qdtp0 = qcpt0 = qcpp0 = qat0 = qap0 = 0.0
qdt2 = qcp2 = dela2 = 0.0


def eqstat2(tl, t, pl, p, dl, d, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
            qcpp, lderiv, latmo, ksaha):
    """Choose and blend the equations of state by regime: composition,
    radiation, Saha (or SCV) below TSCUT, the relativistic gas above it
    with a smooth-step blend between, and OPAL 2006 on top where its
    table covers the point.

    Translates: eqstat2.f#EQSTAT2 (09c3b9b)

    Arguments: tl, pl (log T, log P), x, z, lderiv (also the second
    derivatives), latmo (passed to eqsaha), ksaha (eqsaha's iteration
    counter); dl is a starting guess. t, p, d, beta ... qcpp are only set.
    Array changed in place: fxion (H+, He+, He++ fractions, index 0
    unused).

    Returns (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt,
    qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha).
    """
    global qdtt0, qdtp0, qcpt0, qcpp0, qat0, qap0, qdt2, qcp2, dela2
    c = constants_mod
    nml = nmlpar_mod
    cmp = comp_mod

    def _out():
        return (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt,
                qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha)

    # lsaha: solve Saha (else fully ionized); lonlys: Saha alone, no
    # relativistic gas.
    lsaha = tl < nml.tscut
    lonlys = tl <= (nml.tscut - SCALET)

    # Composition: AMU (ions per unit mass) and the number fractions FX.
    z0 = z
    if nml.ldifz:
        # Metal diffusion: all metals scale together by Z/ZENV.
        rat = z / cmp.zenv
        r1a = (rat - 1.0) * cmp.amuenv
        amu = cmp.amuenv
        for i in range(1, 7):
            amu = amu + r1a * cmp.fxenv[i]
        dxx = (x - cmp.xenv) / ATOMW2[7]
        amu = amu + dxx
        for i in range(8, 12):
            amu = amu + r1a * cmp.fxenv[i]
        dyy = (cmp.xenv + cmp.zenv - x - z) / ATOMW2[12]
        amu = amu + dyy
        raa = rat * cmp.amuenv / amu
        if lsaha:
            for i in range(1, 7):
                fx[i] = raa * cmp.fxenv[i]
            fx[7] = (cmp.fxenv[7] * cmp.amuenv + dxx) / amu
            for i in range(8, 12):
                fx[i] = raa * cmp.fxenv[i]
            fx[12] = (cmp.fxenv[12] * cmp.amuenv + dyy) / amu
    else:
        dfx1 = x - cmp.xenv
        dfx12 = z - cmp.zenv
        if abs(dfx1) + abs(dfx12) < 1.0e-5:
            # The envelope abundances.
            amu = cmp.amuenv
            if lsaha:
                for i in range(1, NZP1 + 1):
                    fx[i] = cmp.fxenv[i]
        else:
            # Any excess Z is taken as carbon 12.
            dfx1 = dfx1 * ATOMWT[1]
            dfx12 = dfx12 * ATOMWT[3]
            dfx4 = (cmp.xenv + cmp.zenv - x - z) * ATOMWT[2]
            amu = cmp.amuenv + dfx1 + dfx4 + dfx12
            amui = 1.0 / amu
            if lsaha:
                fmuenv = cmp.amuenv * amui
                for i in range(1, NZP1 + 1):
                    fx[i] = fmuenv * cmp.fxenv[i]
                fx[6] = fx[6] + dfx12 * amui
                fx[7] = fx[7] + dfx1 * amui
                fx[12] = fx[12] + dfx4 * amui

    # Radiation pressure.
    t = math.exp(c.cln * tl)
    p = math.exp(c.cln * pl)
    beta14 = c.ca3 * (t**2)**2 / p
    beta = 1 - beta14
    if beta < 0.0:
        # Question (see the page): the routine goes on with beta < 0.
        iounits_mod.write(iounits_mod.ishort,
                          f" EQSTAT2: BETA is negative: TL,PL,BETA,1-BETA "
                          f"{tl} {pl} {beta} {1.0 - beta}")
    betai = 1.0 / beta
    # Question (see the page): beta14 is 4(1-beta)/beta here, 1-beta
    # after OPAL.
    beta14 = 4.0 * betai * min(1.0, beta14)
    eta = 0.0

    if not lsaha:
        # Fully ionized gas.
        emu = x * ATOMWT[1] + (1.0 - x) * ATOMWT[4]
        rmu = c.cgas * (amu + emu)
        fxion[1] = 1.0
        fxion[2] = 0.0
        fxion[3] = 1.0
        (dl, d, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
         qcpp) = eqrelv(tl, t, pl, p, dl, d, beta, amu, emu, eta, qdt, qdp,
                        qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp)
    else:
        if nml.lscv:
            # SCV replaces Saha, except for the second derivatives.
            (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
             qcpt, qcpp, ksaha) = eqsaha(
                fx, tl, t, pl, p, dl, d, beta, betai, beta14, fxion, rmu,
                amu, emu, lderiv, latmo, qdt, qdp, qcp, dela, qdtt, qdtp,
                qat, qap, qcpt, qcpp, ksaha)
            lcalc = None
            dl, d, qdt, qdp, qcp, dela, lcalc = eqscve(
                tl, t, p, dl, d, x, z, beta, fxion, qdt, qdp, qcp, dela,
                lcalc)
            # Question (see the page): this block never runs, and would
            # use qdt2, qcp2, dela2, dl1, d1 unset.
            lscvd = False
            if lscvd and lcalc:
                # Unset in the Fortran (uninitialized locals).
                dl1 = d1 = qdt1 = qdp1 = qcp1 = dela1 = lcalc1 = None
                qdp2 = lcalc2 = None
                dpl = 0.010
                dtl = 0.040
                ttl = tl + dtl
                t = 10.0**ttl
                dl1, d1, qdt1, qdp1, qcp1, dela1, lcalc1 = eqscve(
                    ttl, t, p, dl1, d1, x, z, beta, fxion, qdt1, qdp1,
                    qcp1, dela1, lcalc1)
                ttl = tl - dtl
                t = 10.0**ttl
                dl1, d1, qdt1, qdp1, qcp1, dela1, lcalc1 = eqscve(
                    ttl, t, p, dl1, d1, x, z, beta, fxion, qdt1, qdp1,
                    qcp1, dela1, lcalc1)
                dtl2 = 2.0 * dtl
                qdtt1 = (qdt1 - qdt2) / dtl2 / c.cln
                qcpt1 = (math.log10(qcp1) - math.log10(qcp2)) / dtl2
                qat1 = (math.log10(dela1) - math.log10(dela2)) / dtl2
                t = 10.0**tl
                ppl = pl + dpl
                p = 10.0**ppl
                dl1, d1, qdt1, qdp1, qcp1, dela1, lcalc1 = eqscve(
                    tl, t, p, dl1, d1, x, z, beta, fxion, qdt1, qdp1, qcp1,
                    dela1, lcalc1)
                ppl = pl - dpl
                p = 10.0**ppl
                dl1, d1, qdt2, qdp2, qcp2, dela2, lcalc2 = eqscve(
                    tl, t, p, dl1, d1, x, z, beta, fxion, qdt2, qdp2, qcp2,
                    dela2, lcalc2)
                p = 10.0**pl
                dpl2 = 2.0 * dpl
                qdtp1 = (qdt1 - qdt2) / dpl2 / c.cln
                qcpp1 = (math.log10(qcp1) - math.log10(qcp2)) / dpl2
                qap1 = (math.log10(dela1) - math.log10(dela2)) / dpl2
                qdtt = qdtt1
                qdtp = qdtp1
                qat = qat1
                qap = qap1
                qcpp = qcpp1
                qcpt = qcpt1
        else:
            # Saha alone (Prather's EOS).
            (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
             qcpt, qcpp, ksaha) = eqsaha(
                fx, tl, t, pl, p, dl, d, beta, betai, beta14, fxion, rmu,
                amu, emu, lderiv, latmo, qdt, qdp, qcp, dela, qdtt, qdtp,
                qat, qap, qcpt, qcpp, ksaha)
        # Store the Saha values for the blend with the relativistic gas.
        dl0 = dl
        qdt0 = qdt
        qdp0 = qdp
        qcp0 = qcp
        dela0 = dela
        if lderiv:
            qdtt0 = qdtt
            qdtp0 = qdtp
            qcpt0 = qcpt
            qcpp0 = qcpp
            qat0 = qat
            qap0 = qap
        if not lonlys:
            # Transition band: the fully ionized gas too, then blend.
            emu = x * ATOMWT[1] + (1.0 - x) * ATOMWT[4]
            rmu = c.cgas * (amu + emu)
            fxion[1] = 1.0
            fxion[2] = 0.0
            fxion[3] = 1.0
            (dl, d, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
             qcpp) = eqrelv(tl, t, pl, p, dl, d, beta, amu, emu, eta, qdt,
                            qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
                            qcpp)
            # Smooth step W(w) = w^3 (6w^2 - 15w + 10).
            wts = SCALES * (nml.tscut - tl)
            wts2 = wts * wts
            wt0 = wts2 * wts * (6.0 * wts2 - 15.0 * wts + 10.0)
            dl0 = dl0 - dl
            qdt0 = qdt0 - qdt
            qdp0 = qdp0 - qdp
            qcp0 = qcp0 - qcp
            dela0 = dela0 - dela
            # Question (see the page): no differences for qdtt0, qdtp0.
            qcpt0 = qcpt0 - qcpt
            qcpp0 = qcpp0 - qcpp
            qat0 = qat0 - qat
            qap0 = qap0 - qap
            dl = dl + wt0 * dl0
            d = math.exp(c.cln * dl)
            rmu = beta * p / (d * t)
            emu = rmu / c.cgas - amu
            fxion[1] = 1.0 + wt0 * (fxion[1] - 1.0)
            fxion[2] = wt0 * fxion[2]
            fxion[3] = 1.0 + wt0 * (fxion[3] - 1.0)
            qdt = qdt + wt0 * qdt0
            qdp = qdp + wt0 * qdp0
            qcp = qcp + wt0 * qcp0
            dela = dela + wt0 * dela0
            if lderiv:
                qdtt = qdtt + wt0 * qdtt0
                qdtp = qdtp + wt0 * qdtp0
                qcpt = qcpt + wt0 * qcpt0
                qcpp = qcpp + wt0 * qcpp0
                qat = qat + wt0 * qat0
                qap = qap + wt0 * qap0
    # label 200

    # OPAL 2006, replacing or blended with the Yale/SCV result.
    # Question (see the page): the rough ranges are hard-coded.
    if nml.lopale06 and t >= 1.870e3 and t <= 200e6 and dl <= 7.0:
        # oeqos06's alternate return (*998) is its last returned value,
        # alt1: True means it failed, and the Yale/SCV result stands.
        (t, p, d10o, doo, betao, betaio, beta14o, rmuo, amuo, emuo, qdto,
         qdpo, qcpo, delao, alt1) = oeqos06(
            tl, t, pl, p, 0.0, 0.0, x, z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0)
        if alt1:
            return _out()                          # go to 998
        fac, ltab, lramp = eqbound06(t, d10o, 0.0, False, False)
        if not ltab:
            return _out()                          # go to 998
        if not lramp:
            # Fully in the OPAL 2006 table.
            dl = d10o
            d = 10.0**dl
            beta = betao
            betai = 1.0 / beta
            beta14 = 1.0 - beta
            rmu = rmuo
            amu = amuo
            emu = emuo
            qdt = qdto
            qdp = qdpo
            qcp = qcpo
            dela = delao
        else:
            # On the border: ramp between OPAL 2006 and Yale/SCV.
            dl = dl + fac * (d10o - dl)
            d = 10.0**dl
            beta = beta + fac * (betao - beta)
            betai = 1.0 / beta
            beta14 = 1.0 - beta
            rmu = rmu + fac * (rmuo - rmu)
            amu = amu + fac * (amuo - amu)
            emu = emu + fac * (emuo - emu)
            qdt = qdt + fac * (qdto - qdt)
            qdp = qdp + fac * (qdpo - qdp)
            qcp = qcp + fac * (qcpo - qcp)
            dela = dela + fac * (delao - dela)
    # label 998: Z = Z0 (question 6: nothing here changes z).
    z = z0
    return _out()
