"""Python version of setscv.f (YREC 09c3b9b): the SCV equation-of-state
table for the envelope's mixture, from the pure H, He and metal tables."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import scveos_mod
from pyrec.inter3 import inter3

NTS = 63
# DATA: masses of H and He atoms (g), Boltzmann's constant (erg/K).
CMH = 1.67357e-24
CMHE = 6.646442e-24
CBOLTZ = 1.380658e-16

# Locals kept between calls (SAVE): interpolation weights.
qr = np.zeros(4)
ft = np.zeros(4)
ftd = np.zeros(4)
fp = np.zeros(4)
fpd = np.zeros(4)


def setscv():
    """Build the envelope mixture's table tablenv in scveos_mod: dU/dT of
    the pure H and He tables (column 12), the entropy of mixing smix, the
    mixture's density, d ln rho / d ln T and d ln P, c_P and dU/dT
    (additive volumes), then their derivatives by three-point Lagrange
    differentiation (columns 7-12).

    Translates: setscv.f#SETSCV (09c3b9b)

    No arguments; works on scveos_mod (and reads comp_mod.xenv, zenv).
    """
    c = constants_mod
    s = scveos_mod
    tlogx = s.tlogx
    tablex = s.tablex
    tabley = s.tabley
    tablez = s.tablez
    tablenv = s.tablenv
    smix = s.smix
    nptsx = s.nptsx
    xenv = comp_mod.xenv
    zenv = comp_mod.zenv

    # Question (see the page): "READ IN ... TABLES" over a loop that reads
    # nothing.
    for i in range(1, NTS + 1):
        for j in range(1, nptsx[i] + 1):
            # dU/dT from c_P = S (d ln S/d ln T)_P
            #   = dU/dT + P (d ln rho/d ln T)^2 / (rho T d ln rho/d ln P).
            tablex[i, j, 12] = (
                math.exp(c.cln * tablex[i, j, 5]) * tablex[i, j, 9]
                - math.exp(c.cln * (tablex[i, j, 1] - tablex[i, j, 4]
                                    - tlogx[i]))
                * tablex[i, j, 7]**2 / tablex[i, j, 8])
            tabley[i, j, 12] = (
                math.exp(c.cln * tabley[i, j, 5]) * tabley[i, j, 9]
                - math.exp(c.cln * (tabley[i, j, 1] - tabley[i, j, 4]
                                    - tlogx[i]))
                * tabley[i, j, 7]**2 / tabley[i, j, 8])
            # Entropy of mixing: number densities of H and He.
            # Question (see the page): Y = 1 - X here, ignoring Z.
            y = 1.0 - xenv
            ym1 = 1.0 - y
            denx = math.exp(c.cln * tablex[i, j, 4])
            deny = math.exp(c.cln * tabley[i, j, 4])
            dmix = 1.0 / (ym1 / denx + y / deny)
            hnum = (2.0 * ym1 * dmix / CMH
                    / (1.0 + 3.0 * tablex[i, j, 2] + tablex[i, j, 3]))
            henum = (3.0 * y * dmix / CMHE
                     / (1.0 + 2.0 * tabley[i, j, 2] + tabley[i, j, 3]))
            # Electron number densities from H and He.
            hnume = 0.5 * hnum * (1.0 - tablex[i, j, 2] - tablex[i, j, 3])
            henume = c.cc13 * henum * (2.0 - 2.0 * tabley[i, j, 2]
                                       - tabley[i, j, 3])
            tote = hnume + henume
            hnume = max(hnume, 1.0)
            henume = max(henume, 1.0)
            tote = max(tote, 1.0)
            if hnum == 0.0 or henum == 0.0:
                smix[i, j] = 0.0
            else:
                be = (CMH / CMHE) * (y / ym1)
                ga = (1.5 * (1.0 + tablex[i, j, 3] + 3.0 * tablex[i, j, 2])
                      / (1.0 + 2.0 * tabley[i, j, 2] + tabley[i, j, 3]))
                xd = 1.0 - tablex[i, j, 2] - tablex[i, j, 3]
                yd = 2.0 - 2.0 * tabley[i, j, 2] - tabley[i, j, 3]
                if xd > 0.0:
                    de = 1.5 * yd * be * ga / xd
                else:
                    de = 0.0
                if yd > 0.0:
                    dei = xd / (1.5 * yd * be * ga)
                else:
                    dei = 0.0
                tnum = hnum + henum
                smix[i, j] = (ym1 / CMH * 2.0
                              / (1.0 + tablex[i, j, 3] + 3.0 * tablex[i, j, 2])
                              * (math.log(1.0 + be * ga)
                                 - hnume / tnum * math.log(1.0 + de)
                                 + be * ga * math.log(1.0 + 1.0 / be / ga)
                                 - henume / tnum * math.log(1.0 + dei)))
                smix[i, j] = CBOLTZ * smix[i, j]

    # The equation of state of the surface mixture.
    for i in range(1, NTS + 1):
        tl = tlogx[i]
        t = 10.0**tl
        s.idtscv = _three_point_start(i, NTS)
        idt = s.idtscv
        for k in range(1, 4):
            qr[k] = c.cln * tlogx[k + idt - 1]
        tln = c.cln * tl
        inter3(qr, ft, ftd, tln)
        for j in range(1, nptsx[i] + 1):
            pl = tablex[i, j, 1]
            dl = tablex[i, j, 4]
            x = xenv
            z = zenv
            y = 1.0 - x - z
            p = 10.0**pl
            # Radiation pressure.
            # Question (see the page): P and PL are not used afterwards.
            prad = c.ca3 * (t**2)**2
            pgas = p
            p = p + prad
            pl = math.log10(p)
            # Densities add inversely: 1/rho = X/rho_X + Y/rho_Y + Z/rho_Z.
            denx = math.exp(c.cln * tablex[i, j, 4])
            deny = math.exp(c.cln * tabley[i, j, 4])
            denz = math.exp(c.cln * tablez[i, j, 4])
            d = 1.0 / (x / denx + y / deny + z / denz)
            dl = math.log10(d)
            # d ln rho / d ln P.
            qdp0 = (tablex[i, j, 8] * x * (d / denx)
                    + tabley[i, j, 8] * y * (d / deny)
                    + tablez[i, j, 13] * z * (d / denz))
            # d ln rho / d ln T.
            qdt0 = (x * (d / denx) * tablex[i, j, 7]
                    + y * (d / deny) * tabley[i, j, 7]
                    + z * (d / denz) * tablez[i, j, 10])
            qpt0 = -qdt0 / qdp0
            sx = math.exp(c.cln * tablex[i, j, 5])
            sy = math.exp(c.cln * tabley[i, j, 5])
            # Entropy: additive, plus the entropy of mixing.
            s0 = x * sx + y * sy + smix[i, j]
            # Question (see the page): row J of the neighbor temperatures
            # may lie beyond their NPTSX.
            qsmixt = (ftd[1] * smix[idt, j] + ftd[2] * smix[idt + 1, j]
                      + ftd[3] * smix[idt + 2, j])
            # d ln S / d ln T.
            qst0 = (x * sx * tablex[i, j, 9] + y * sy * tabley[i, j, 9]
                    + qsmixt) / s0
            # c_P (gas pressure only), with the metals' term.
            qcp0 = s0 * qst0 + z * (
                math.exp(c.cln * tablez[i, j, 7])
                + pgas * tablez[i, j, 10]**2 / tablez[i, j, 13] / denz / t)
            # dU/dT.
            qut = qcp0 - pgas * qdt0**2 / qdp0 / d / t
            tablenv[i, j, 1] = tablex[i, j, 1]
            tablenv[i, j, 2] = dl
            tablenv[i, j, 3] = qdt0
            tablenv[i, j, 4] = qdp0
            tablenv[i, j, 5] = qcp0
            tablenv[i, j, 6] = qut

    # The derivatives, numerically.
    # Question (see the page): the comment names DEL(AD), the code dU/dT.
    for i in range(1, NTS + 1):
        s.idtscv = _three_point_start(i, NTS)
        idt = s.idtscv
        for k in range(1, 4):
            qr[k] = c.cln * tlogx[k + idt - 1]
        tl = c.cln * tlogx[i]
        inter3(qr, ft, ftd, tl)
        t = math.exp(tl)
        for j in range(1, nptsx[i] + 1):
            s.idp = _three_point_start(j, nptsx[i])
            idp = s.idp
            for k in range(1, 4):
                qr[k] = c.cln * tablenv[i, idp + k - 1, 1]
            pl = c.cln * tablenv[i, j, 1]
            inter3(qr, fp, fpd, pl)
            # Derivatives of d ln rho / d ln T.
            qdtp = (tablenv[i, idp, 3] * fpd[1]
                    + tablenv[i, idp + 1, 3] * fpd[2]
                    + tablenv[i, idp + 2, 3] * fpd[3])
            # Question (see the page): no MIN(NPTSX, J) protection here.
            qdtt = (tablenv[idt, j, 3] * ftd[1]
                    + tablenv[idt + 1, j, 3] * ftd[2]
                    + tablenv[idt + 2, j, 3] * ftd[3])
            # Derivatives of ln c_P.
            for k in range(1, 4):
                qr[k] = math.log(tablenv[i, idp + k - 1, 5])
            # Question (see the page): QCPP is not stored.
            qcpp = qr[1] * fpd[1] + qr[2] * fpd[2] + qr[3] * fpd[3]
            for k in range(1, 4):
                ii = idt + k - 1
                jj = min(nptsx[ii], j)
                qr[k] = math.log(tablenv[ii, jj, 5])
            qcpt = qr[1] * ftd[1] + qr[2] * ftd[2] + qr[3] * ftd[3]
            # Derivatives of dU/dT.
            for k in range(1, 4):
                ii = idt + k - 1
                jj = min(nptsx[ii], j)
                qr[k] = tablenv[ii, jj, 6]
                # Question (see the page): this line undoes the one above.
                qr[k] = tablenv[idt + k - 1, j, 6]
            qat = qr[1] * ftd[1] + qr[2] * ftd[2] + qr[3] * ftd[3]
            for k in range(1, 4):
                qr[k] = tablenv[i, idp + k - 1, 6]
            qap = qr[1] * fpd[1] + qr[2] * fpd[2] + qr[3] * fpd[3]
            # d/d ln rho (d ln P / d ln rho).
            for k in range(1, 4):
                qr[k] = c.cln * tablenv[i, idp + k - 1, 2]
            dln = c.cln * tablenv[i, j, 2]
            inter3(qr, fp, fpd, dln)
            qpdd = (fpd[1] / tablenv[i, idp, 4]
                    + fpd[2] / tablenv[i, idp + 1, 4]
                    + fpd[3] / tablenv[i, idp + 2, 4])
            tablenv[i, j, 7] = qdtp
            tablenv[i, j, 8] = qdtt
            tablenv[i, j, 9] = qpdd
            tablenv[i, j, 10] = qcpt
            tablenv[i, j, 11] = qap
            tablenv[i, j, 12] = qat


def _three_point_start(i, n):
    """The first of three grid points centered on i, kept inside 1..n
    (the IF block the Fortran repeats for IDTSCV and IDP)."""
    if i == 1:
        return 1
    elif i == n:
        return i - 2
    else:
        return i - 1
