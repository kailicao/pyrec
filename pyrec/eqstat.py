"""Python version of eqstat.f (YREC 09c3b9b): the front end to the
equation of state, with numerical second derivatives."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec.eqstat2 import eqstat2

# SAVE: the counters of the derivative calls. ksahax is never set from
# ksaha, because of the misspelled kssahax (see the page's question 1).
ksahax = 0
kssahax = 0


def eqstat(tl, t, pl, p, dl, d, x, z, beta, betai, beta14, fxion, rmu,
           amu, emu, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
           qcpp, lderiv, latmo, ksaha):
    """Call EQSTAT2, the equation of state, at the point; with lderiv and
    nmlpar_mod.lnumderiv, first compute the second derivatives by central
    differences from four calls at log T +- 3e-4 and log P +- 1.5e-3.

    Translates: eqstat.f#EQSTAT (09c3b9b)

    Arguments: tl, pl (log T, log P), x, z, lderiv (second derivatives
    wanted), latmo (atmosphere), ksaha (EQSAHA's iteration counter,
    updated); dl, d are starting guesses. t, p, beta ... qcpp are set.
    Array changed in place: fxion (index 0 unused).

    Returns (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt,
    qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha). The Fortran
    also assigns TL and PL, but only back to their own values, so they are
    not returned. Without lderiv, qdtt ... qcpp come back as passed in.
    """
    global ksahax, kssahax
    cln = constants_mod.cln

    plo = pl
    tlo = tl
    latmox = latmo
    # Question (see the page): the misspelled kssahax leaves ksahax unset.
    kssahax = ksaha

    if lderiv and nmlpar_mod.lnumderiv:
        # Central differences, error of order h**2.
        dpl = 0.150      # the approximate table intervals
        dtl = 0.030
        dpl = dpl * 0.01  # empirically, scaling the intervals down by 100
        dtl = dtl * 0.01  # seemed to work best
        dpl2 = 2.0 * dpl
        dtl2 = 2.0 * dtl

        dl1 = dl  # start the table lookups from the called values
        d1 = d

        # Derivatives with respect to T.
        ttl = tl + dtl
        t = 10.0 ** ttl
        lderiv2 = False
        (t, p, dl1, d1, beta, betai, beta14, rmu, amu, emu, eta, qdt1,
         qdp1, qcp1, dela1, _, _, _, _, _, _, ksahax) = eqstat2(
            ttl, t, pl, p, dl1, d1, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, lderiv2, latmox, ksahax)

        ttl = tl - dtl
        t = 10.0 ** ttl
        lderiv2 = False
        (t, p, dl1, d1, beta, betai, beta14, rmu, amu, emu, eta, qdt2,
         qdp2, qcp2, dela2, _, _, _, _, _, _, ksahax) = eqstat2(
            ttl, t, pl, p, dl1, d1, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, lderiv2, latmox, ksahax)
        qdtt1 = (qdt1 - qdt2) / dtl2 / cln
        qcpt1 = (math.log10(qcp1) - math.log10(qcp2)) / dtl2
        qat1 = (math.log10(dela1) - math.log10(dela2)) / dtl2
        tl = tlo
        t = 10.0 ** tlo  # restore the original T

        # Derivatives with respect to P.
        ppl = pl + dpl
        p = 10.0 ** ppl
        lderiv2 = False
        (t, p, dl1, d1, beta, betai, beta14, rmu, amu, emu, eta, qdt1,
         qdp1, qcp1, dela1, _, _, _, _, _, _, ksahax) = eqstat2(
            tl, t, ppl, p, dl1, d1, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, lderiv2, latmox, ksahax)
        ppl = pl - dpl
        p = 10.0 ** ppl
        lderiv2 = False
        (t, p, dl1, d1, beta, betai, beta14, rmu, amu, emu, eta, qdt2,
         qdp2, qcp2, dela2, _, _, _, _, _, _, ksahax) = eqstat2(
            tl, t, ppl, p, dl1, d1, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, lderiv2, latmox, ksahax)
        pl = plo
        p = 10.0 ** plo  # restore the original P
        qdtp1 = (qdt1 - qdt2) / dpl2 / cln
        qcpp1 = (math.log10(qcp1) - math.log10(qcp2)) / dpl2
        qap1 = (math.log10(dela1) - math.log10(dela2)) / dpl2

        qdtt = qdtt1
        qdtp = qdtp1
        qat = qat1
        qap = qap1
        qcpp = qcpp1
        qcpt = qcpt1

    tl = tlo  # restore the original TL and PL
    pl = plo

    if lderiv and not nmlpar_mod.lnumderiv:
        # Derivatives wanted, none numerical: EQSTAT2's analytic ones.
        lderiv2 = True
        (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
         qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat2(
            tl, t, pl, p, dl, d, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
            qcpt, qcpp, lderiv2, latmo, ksaha)
    else:
        # Numerical derivatives already in hand, or none needed: the
        # second derivatives go to dummies so they cannot overwrite them.
        lderiv2 = False
        (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
         qcp, dela, _, _, _, _, _, _, ksaha) = eqstat2(
            tl, t, pl, p, dl, d, x, z, beta, betai, beta14, fxion, rmu,
            amu, emu, eta, qdt, qdp, qcp, dela, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, lderiv2, latmo, ksaha)

    return (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
            qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha)

