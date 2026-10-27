"""Python version of eqrelv.f (YREC 09c3b9b): the equation of state for
fully ionized gas with partially degenerate, possibly relativistic
electrons."""

import math

import numpy as np

from pyrec import ccr_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod


def eqrelv(tl, t, pl, p, dl, d, beta, amu, emu, eta, qdt, qdp, qcp,
           dela, qdtt, qdtp, qat, qap, qcpt, qcpp):
    """Find the density at which ions, electrons (from the Fermi table) and
    radiation, plus an optional Debye-Hueckel term, give the pressure p,
    then the thermodynamic derivatives the structure equations need.

    Translates: eqrelv.f#EQRELV (09c3b9b)

    Arguments: tl, t (log10 T and T); pl, p (log10 P and P); dl, the
    starting guess of log10 rho; d; beta (P_gas/P); amu, emu (ions and
    electrons per unit mass, in units of 1/m_H); eta, qdt, qdp, qcp, dela,
    qdtt, qdtp, qat, qap, qcpt, qcpp (set only).

    Returns (dl, d, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
    qcpp): the density and its log, the electron degeneracy, d ln rho /
    d ln T and / d ln P, c_P, nabla_ad, and their logarithmic derivatives.
    If the iteration does not converge, a message is written and the
    arguments come back as they were passed (eta possibly set).
    """
    c = constants_mod
    dval = ccr_mod.dval
    etat8 = ccr_mod.etat8
    dtable = ccr_mod.dtable
    iden = ccr_mod.iden
    ldh = nmlpar_mod.ldh
    etadh0 = nmlpar_mod.etadh0
    etadh1 = nmlpar_mod.etadh1
    df1 = np.zeros((5 + 1, 3 + 1))
    df2 = np.zeros((5 + 1, 3 + 1))
    ff = np.zeros((5 + 1, 3 + 1))
    ffe = np.zeros(5 + 1)

    pinv = 1.0 / p
    beta1 = 1.0 - beta
    tl8 = tl
    pl8 = pl
    pr = beta1 * p
    ramu = c.cgas * amu
    ramut = ramu * t
    emul = math.log10(emu)
    emul8 = emul
    # Newton interpolation in log T: three columns from jt1.
    # Question (see the page): forced into the table, extrapolating silently.
    jt1 = int(5.0 * (tl - 5.0))
    jt1 = max(1, min(jt1, 18))
    jt2 = jt1 + 1
    jt3 = jt1 + 2
    ttud = 5.0 * (tl8 - (5.0 + 0.20 * float(jt1)))
    ttcu = 0.50 * ttud * (ttud - 1.0)
    ttdcu = ttud - 0.50
    dml = dl + emul
    nden = 1
    id1p = -100
    if ldh:
        almix = 0.0
        for ii in range(1, 19):
            almix = almix + debhu_mod.dhnue[ii] * debhu_mod.zdh[ii]
        cmfdh0 = (2.0 * debhu_mod.xxdh + 1.5 * debhu_mod.yydh
                  + 0.5 * (debhu_mod.zzdh + almix))
        cmfdh1 = debhu_mod.xxdh + debhu_mod.yydh + 0.5 * almix

    # Newton iteration on log rho for the pressure balance.
    while True:  # label 10
        xx = dml - 1.50 * tl8
        pr6 = 20.0 * (xx - dval[1]) + 1.0
        kk = int(pr6)
        kk = min(261, max(1, kk))
        id1 = iden[kk]
        if id1 != id1p:  # else GO TO 30: same three table rows
            id1p = id1
            id2 = id1 + 1
            id3 = id1 + 2
            pr9 = 1.0 / (dval[id3] - dval[id2])
            pr8 = 1.0 / (dval[id2] - dval[id1])
            pr7 = 1.0 / (dval[id1] - dval[id3])
            cden1 = -pr7 * pr8
            cden2 = -pr8 * pr9
            cden3 = -pr9 * pr7
            for nn in range(1, 4):  # label 20
                id_ = id1 + nn - 1
                df1[1, nn] = dtable[1, id_, jt2] - dtable[1, id_, jt1]
                df2[1, nn] = dtable[1, id_, jt3] - dtable[1, id_, jt2] - df1[1, nn]
                ff[1, nn] = dtable[1, id_, jt1] + ttud * df1[1, nn] + ttcu * df2[1, nn]
            dpel2 = 2.0 * (ff[1, 1] * cden1 + ff[1, 2] * cden2 + ff[1, 3] * cden3)
        # label 30: Lagrange interpolation in x.
        dx1 = xx - dval[id1]
        dx2 = xx - dval[id2]
        dx3 = xx - dval[id3]
        cl1 = dx2 * dx3 * cden1
        cl2 = dx3 * dx1 * cden2
        cl3 = dx1 * dx2 * cden3
        dcl1 = (dx2 + dx3) * cden1
        dcl2 = (dx3 + dx1) * cden2
        dcl3 = (dx1 + dx2) * cden3
        if ldh:
            # Debye-Hueckel composition factor, ramped in eta.
            eta = cl1 * etat8[id1] + cl2 * etat8[id2] + cl3 * etat8[id3]
            if eta < etadh0:
                cmfdh = cmfdh0 * math.sqrt(cmfdh0)
            elif eta > etadh1:
                cmfdh = cmfdh0 * math.sqrt(cmfdh1)
            else:
                t0 = cmfdh0 * math.sqrt(cmfdh0)
                t1 = cmfdh0 * math.sqrt(cmfdh1)
                slope = (t1 - t0) / (etadh1 - etadh0)
                cmfdh = t0 + (eta - etadh0) * slope
        pel = cl1 * ff[1, 1] + cl2 * ff[1, 2] + cl3 * ff[1, 3]
        dpel1 = dcl1 * ff[1, 1] + dcl2 * ff[1, 2] + dcl3 * ff[1, 3]
        dl8 = dml - emul8
        d8 = math.exp(c.cln * dl8)
        pa = ramut * d8
        pe = math.exp(c.cln * pel)
        if ldh:
            pdh = debhu_mod.cdh * cmfdh * math.sqrt(d8 / t) * d8
        pt = pa + pe + pr
        if ldh:
            pt = pt + pdh
            pdhp = pdh / pt
        ptl = math.log10(pt)
        pr7 = (pa + pe * dpel1) / pt
        if ldh:
            pr7 = pr7 + 1.5 * pdhp
        pr8 = (ptl - pl8) / pr7
        pr9 = c.cln * ((pa + pe * (dpel1**2 + dpel2 / c.cln)) / (pt * pr7) - pr7)
        if ldh:
            pr9 = pr9 + 2.25 * pdhp / pr7
        # Second-order correction 0.5*pr9*pr8**2, formed through logs.
        factor = abs(pr8)
        if factor < 1e-15:
            factor = 1e-15
        pr8l = math.log10(factor)
        pr9l = math.log10(abs(pr9))
        factl = math.log10(0.5) + pr9l + 2.0 * pr8l
        fact = math.exp(c.cln * factl)
        fact = math.copysign(fact, pr9)
        corr = -pr8 - fact
        if corr > 4.0:
            corr = 4.0
        if corr < -4.0:
            corr = -4.0
        if abs(corr) >= 1.0e-08:
            dml = dml + corr
            nden = nden + 1
            if nden <= 20:
                continue
            # Question (see the page): returns unconverged, with dl, d and
            # the derivatives as the caller passed them.
            iounits_mod.write(
                iounits_mod.ishort,
                f"EQRELV: Did not Converge: T,P,Pcalc,Dcalc,CORR"
                f"{tl:10.6f}{pl:10.6f}{ptl:10.6f}{dml:10.6f}{corr:20.12f}")
            return (dl, d, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
                    qcpt, qcpp)
        break

    dl = dl8
    d = d8
    # The other table columns: log U_e, and the derivatives.
    for kk in range(2, 6):  # label 60
        for nn in range(1, 4):  # label 50
            id_ = id1 + nn - 1
            df1[kk, nn] = dtable[kk, id_, jt2] - dtable[kk, id_, jt1]
            df2[kk, nn] = (dtable[kk, id_, jt3] - dtable[kk, id_, jt2]) - df1[kk, nn]
            ff[kk, nn] = dtable[kk, id_, jt1] + ttud * df1[kk, nn] + ttcu * df2[kk, nn]
        ffe[kk] = cl1 * ff[kk, 1] + cl2 * ff[kk, 2] + cl3 * ff[kk, 3]

    # Pressure derivatives.
    peq = math.exp(c.cln * ffe[3])
    qped = ffe[4]
    qqpedd = (dcl1 * ff[4, 1] + dcl2 * ff[4, 2] + dcl3 * ff[4, 3]) * c.clni
    pr1 = df1[3, 1] + ttdcu * df2[3, 1]
    pr2 = df1[3, 2] + ttdcu * df2[3, 2]
    pr3 = df1[3, 3] + ttdcu * df2[3, 3]
    pr4 = dcl1 * ff[3, 1] + dcl2 * ff[3, 2] + dcl3 * ff[3, 3]
    pr5 = 5.0 * (cl1 * pr1 + cl2 * pr2 + cl3 * pr3)
    qpeqd = peq * pr4
    qpeqt = peq * (pr5 - 1.50 * pr4)
    qpt = (pa + peq) * pinv + 4.0 * beta1
    if ldh:
        qpt = qpt - 0.5 * pdhp
    qpd = (pa + pe * qped) * pinv
    if ldh:
        qpd = qpd + 1.5 * pdhp
    qqptt = (pa + qpeqt) * pinv + 16.0 * beta1 - qpt**2
    if ldh:
        qqptt = qqptt + (0.5 * qpt + 0.25) * pdhp
    qqpdt = (pa + qpeqd) * pinv - qpt * qpd
    if ldh:
        qqpdt = qqpdt + (0.5 * qpt - 0.75) * pdhp
    qqpdd = (pa + pe * (qped**2 + qqpedd)) * pinv - qpd**2
    if ldh:
        qqpdd = qqpdd + (-1.5 * qpd + 2.25) * pdhp

    # Density derivatives, by inversion.
    qdp = 1.0 / qpd
    qdt = -qpt * qdp
    pr1 = qdp**2
    qdpp = -qqpdd * qdp * pr1
    qdtp = -qpt * qdpp - qqpdt * pr1
    qdtt = -qpt * qdtp - qdp * (qqptt + qdt * qqpdt)

    # Internal energy and its derivatives.
    ue = math.exp(c.cln * (ffe[2] + emul))
    ueq = math.exp(c.cln * (ffe[5] + emul))
    pr1 = df1[5, 1] + ttdcu * df2[5, 1]
    pr2 = df1[5, 2] + ttdcu * df2[5, 2]
    pr3 = df1[5, 3] + ttdcu * df2[5, 3]
    pr4 = dcl1 * ff[5, 1] + dcl2 * ff[5, 2] + dcl3 * ff[5, 3]
    pr5 = 5.0 * (cl1 * pr1 + cl2 * pr2 + cl3 * pr3)
    queqt = ueq * (pr5 - 1.50 * pr4)
    pd = p / d
    pdt = pd / t
    ua = 1.50 * ramut
    ur = 3.0 * beta1 * pd
    if ldh:
        udh = 3.0 * pdh / d8
    u = ua + ue + ur
    if ldh:
        u = u + udh
        udhu = udh / u
    ui = 1.0 / u
    pdu = pd * ui
    qutd = (ua + ueq + 4.0 * ur) * ui
    if ldh:
        qutd = qutd - 0.5 * udhu
    qudt = pdu * (1.0 - qpt)
    if ldh:
        qudt = qudt + 0.5 * udhu
    qutp = qutd + qdt * qudt
    qqutt = (ua + queqt + 16.0 * ur) * ui - qutd**2
    if ldh:
        qqutt = qqutt + (0.5 * qutd + 0.25) * udhu
    qqudt = qudt * (qpt - qutd) - pdu * qqptt
    if ldh:
        qqudt = qqudt + (0.5 * qudt - 0.25) * udhu

    # c_P, nabla_ad and their derivatives.
    pr1 = -pdt * qdt
    pr2 = u / t
    pr3 = -(1.0 + qdt) + qdtt / qdt
    pr4 = (1.0 - qdp) + qdtp / qdt
    acp = pr1 * qpt
    bcp = pr2 * qutd
    qcp = acp + bcp
    qcpi = 1.0 / qcp
    qcpt = (acp * pr3 + pr1 * (qqptt + qdt * qqpdt)
            + bcp * (qutp - 1.0) + pr2 * (qqutt + qdt * qqudt))
    qcpp = acp * pr4 + qdp * (pr1 * qqpdt + bcp * qudt + pr2 * qqudt)
    qcpt = qcpt * qcpi
    qcpp = qcpp * qcpi
    dela = pr1 * qcpi
    qat = pr3 - qcpt
    qap = pr4 - qcpp
    eta = cl1 * etat8[id1] + cl2 * etat8[id2] + cl3 * etat8[id3]

    return (dl, d, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt,
            qcpp)
