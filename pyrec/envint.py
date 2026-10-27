"""Python version of envint.f (YREC development, 09c3b9b): the atmosphere
and envelope integrations, from the surface down to the fitting point."""

import math

import numpy as np

from pyrec import acdpth_mod
from pyrec import atmos_mod
from pyrec import atmprt_mod
from pyrec import atmstruct_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import envprt_mod
from pyrec import envstruct_mod
from pyrec import intpar_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec import pualpha_mod
from pyrec import pulse1_mod
from pyrec import pulse2_mod
from pyrec import qenvmix_mod
from pyrec.alsurfp import alsurfp
from pyrec.bsstep import bsstep
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.hra import hra
from pyrec.kcsurfp import kcsurfp
from pyrec.qatm import qatm
from pyrec.qenv import qenv
from pyrec.surfp import surfp
from pyrec.taucal import taucal

MAXSTP = 2000             # "DBG changed MAXSTEP from 200 to 2000"
TINY = 1.0e-30
JSON = dimensions_mod.json

# End-of-data marker for the pulsation files (DATA XYZ/22*99.99D0/).
XYZ = np.full(23, 99.99)

# Locals kept between calls (SAVE) whose old values can be read: the
# pulsation variables (question 5), the step counts, last mass difference
# and errors printed at the end (see update-findings.md), and ENVG, which
# is never set (question 2).
taup = 0.0
op = 0.0
dp = 0.0
nok = 0
nbad = 0
sdif = 0.0
y = np.zeros(4)
errsum = np.zeros(4)
envg = np.zeros(JSON + 1)
envs1 = np.zeros(JSON + 1)
envs2 = np.zeros(JSON + 1)
edel1 = np.zeros(JSON + 1)
edel2 = np.zeros(JSON + 1)


def _pulse_lines(unit, values, per_line):
    """Write values to a pulsation file, per_line to a line (FORMATs 5001:
    5E16.9 and 6001: 5E23.16, as 5, 5, 5, 4; 6003: 6E23.16, as 6, 6, 6,
    4). Python's E format prints d.ddd where Fortran's prints 0.ddd."""
    width, digits = (16, 9) if per_line == 0 else (23, 16)
    per_line = per_line or 5
    lines = []
    for k in range(0, len(values), per_line):
        lines.append("".join(f"{v:{width}.{digits}E}"
                             for v in values[k:k + per_line]))
    iounits_mod.write(unit, "\n".join(lines))


def _fmt260(eg, de):
    """FORMAT 260: one envelope point of the .store table."""
    e = envprt_mod
    q = pulse2_mod
    chrh = 1.0 / q.qqdp
    cht = -chrh * q.qqdt
    cv = q.qqcp - math.exp(constants_mod.cln * (e.ep - e.ed - e.et)) \
        * cht**2 / chrh
    gam1 = chrh * q.qqcp / cv
    return (f" {eg:10.3E}{e.ep:10.7f}{e.et:10.7f}{de:14.7E}{e.es:14.7E}"
            f"{e.ed:10.7f}{e.eo:10.2E}{e.ebeta:7.3f}{e.edel[1]:9.2E}"
            f"{e.edel[2]:6.3f}{e.edel[3]:6.3f}{e.efxion[1]:6.3f}"
            f"{e.efxion[2]:6.3f}{e.efxion[3]:6.3f}{e.evel:9.2E}"
            f"{gam1:12.8f}{q.qqdp:12.8f}")


def envint(b, fpl, ftl, gl, hstot, ie, lprt, lsbc, plim, rl, teffl,
           x, z, estore, istore, katm, kenv, ksaha, ps, rs, ts,
           lpulpt, cmixl):
    """Integrate the atmosphere (gray T-tau, or read from a table) and the
    envelope (in log P) down to the fitting mass; return the fitting
    point's P, T, r, store the envelope profile in envstruct_mod, and find
    the base of the convection zone and its overturn time.

    Translates: envint.f#ENVINT (09c3b9b)

    Arguments: b (L), fpl, ftl (rotational corrections), gl (log g),
    hstot (log M), ie (the point of the envelope triangle), lprt, lsbc
    (compute the boundary condition), plim (the expected end pressure), rl
    (log R), teffl, x, z, lpulpt (pulsation output), cmixl (the mixing
    length, relayed to qenv through qenvmix_mod); istore and the counters
    katm, kenv, ksaha. Arrays changed in place: ps, rs, ts (log P, log T,
    log r at the fitting point, entry ie, set only with lsbc), estore (the
    state where T first reaches TENV).

    Returns (rl, istore, katm, kenv, ksaha): rl as QENV last set it
    (log r at the last envelope point; as passed in if QENV was not
    called).
    """
    global taup, op, dp, nok, nbad, sdif
    c = constants_mod
    nml = nmlpar_mod
    io = iounits_mod
    q = pulse2_mod
    apr = atmprt_mod
    epr = envprt_mod
    env = envstruct_mod
    atm = atmstruct_mod
    acd = acdpth_mod
    dydx = np.zeros(4)
    yscal = np.zeros(4)
    ystart = np.zeros(4)
    err = np.zeros(4)
    fxion = np.zeros(4)

    # QENV is called through BSSTEP with a fixed argument list, so CMIXL
    # reaches it through qenvmix_mod.
    qenvmix_mod.cmixlq = cmixl
    # DBG pulse: turn on the derivative calculator.
    if lpulpt and lprt:
        pulse1_mod.lpumod = True
        lderiv = True
        latmo = False
    else:
        pulse1_mod.lpumod = False
        lderiv = False
        latmo = True
    locond = False
    # JVS 10/07/13: always calculate derivatives.
    lderiv = True

    if lprt and nml.lstatm:
        if not nml.lstch:
            io.write(io.istor, "\n******** ATMOSPHERE BEGIN ********")  # 60

    # Tabulated atmospheres give the photospheric pressure directly.
    tabulated = False
    if nml.kttau == 3:
        # Kurucz atmospheres.
        if nml.lstch:
            nml.lstatm = False
        surfp(teffl, gl, lprt and nml.lstatm)
        tabulated = True
    elif nml.kttau == 5:
        # Kurucz (Castelli-Kurucz 2004) atmospheres.
        if nml.lstch:
            nml.lstatm = False
        kcsurfp(teffl, gl, lprt and nml.lstatm)
        tabulated = True
    elif nml.kttau == 4:
        # Allard & Hauschildt atmospheres.
        if nml.lstch:
            nml.lstatm = False
        lalfail = alsurfp(teffl, gl, lprt and nml.lstatm, False)
        if lalfail:
            # Question (see the page): gray atmospheres for the rest of
            # the run (Teff above the table, or g outside it).
            nml.kttau = 0
            nml.lttau = True
            print(" ENVINT: Change to gray atmosphere (KTTAU=0)")
            io.write(io.ishort, " ENVINT: Change to gray atmosphere (KTTAU=0)")
        else:
            tabulated = True

    if not tabulated:
        # label 2: gray atmosphere. Guess T at an optical depth near zero.
        errsum[1] = 0.0
        if nml.kttau == 0:
            tl = teffl - 0.031235 + 0.25 * math.log10(c.cc23)
        elif nml.kttau == 1:
            tl = teffl - 0.031235 + 0.25 * math.log10(0.550)
        elif nml.kttau == 2:
            tl = teffl + hra(c.cc23) - atmos_mod.hras
        # Very occasionally the start (x0) lies past the end (xlim): the
        # starting density is then divided by 10 and the start retried.
        # Question (see the page): no limit on retries.
        atmdx = nml.atmd0
        while True:                                # label 1998
            t = math.exp(c.cln * tl)
            # P at this T and the chosen starting density.
            pl = math.log10((c.cgas * atmdx + c.ca3 * t**3) * t)
            # The optical depth (x0 = log tau) where the integration begins.
            (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta,
             qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp,
             ksaha) = eqstat(tl, t, pl, 0.0, 0.0, 0.0, x, z, 0.0, 0.0, 0.0,
                             fxion, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                             0.0, 0.0, 0.0, 0.0, 0.0, 0.0, lderiv, latmo,
                             ksaha)
            o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
            x0 = pl - gl + math.log10(o)
            y[1] = pl
            dydx[1] = math.exp(c.cln * (gl + x0 - ol - pl))
            if lprt and nml.lstatm:
                if not nml.lstch:
                    io.write(io.istor,
                             "      TAU         P          T          D"
                             "           O         HII  HEII HEIII   SAHA"
                             "   KATM")                # format 10
                beta = 1.0 - c.ca3 * (t**2)**2 / p
                chrh = 1.0 / qdp
                cht = -chrh * qdt
                # Question (see the page): qcpp where qcp is meant.
                cv = qcpp - math.exp(c.cln * (pl - dl - tl)) * cht**2 / chrh
                gam1 = chrh * qcpp / cv
                if not nml.lstch:
                    io.write(io.istor,
                             f" {x0:11.7f}{pl:11.7f}{tl:11.7f}{dl:11.7f}"
                             f"{o:14.7E}{fxion[1]:6.3f}{fxion[2]:6.3f}"
                             f"{fxion[3]:6.3f}{ksaha:7d}{katm:7d}"
                             f"{gam1:12.8f}{qdp:12.8f}{qdt:12.8f}"
                             f"{beta:12.8f}{qcpp:12.5E}{cv:12.5E}")  # 20
            # DBG pulse: the initial point for pulsation.
            if lpulpt and lprt:
                q.qqed = 0.0
                qesum = 0.0
                q.qqet = 0.0
                q.qfs = 1.0
                taup = math.exp(c.cln * x0)
                op = o
                dp = math.exp(c.cln * dl)
                deltr = 0.0
                pelpf = c.cgas * t * d * emu
                # From the first lines of TPGRAD.
                delr = o * b * math.exp(c.cln * (pl - hstot - 4.0 * tl
                                                 + c.clsunl - c.cgl
                                                 + c.cdelrl)) * ftl / fpl
                del_ = delr if delr - dela <= 1.0e-6 else dela
                if nml.ipver == 1:
                    # MHP 10/02 typo: QED should have been QQED.
                    _pulse_lines(io.iopatm,
                                 [deltr, q.qfs, b, tl, dl, pl, qesum, o, qdp,
                                  q.qqed, q.qqet, qod, qot, del_, dela, qcp,
                                  rmu, qdt, pelpf], 0)        # format 5001
                elif nml.ipver in (2, 3):
                    _pulse_lines(io.iopatm,
                                 [deltr, q.qfs, b, tl, dl, pl, taup, o, qdp,
                                  q.qqed, q.qqet, qod, qot, del_, dela, qcp,
                                  rmu, qdt, pelpf], 5)        # format 6001
            # Integrate dP/dtau from this tau to tau = 2/3.
            nv = 1
            eps = nml.atmerr
            if nml.kttau == 1:
                # Krishna Swamy: tau = 0.312156330 at Teff.
                xlim = -0.505627854
            else:
                # tau = 2/3 at Teff.
                xlim = -0.176091259
            if x0 > xlim:
                io.write(io.ishort, f" ENVINT: X0>XLIM, X0,XLIM: {x0} {xlim}")
                io.write(io.ishort,
                         " ENVINT: get new X0 by dividing ATMD0 by 10")
                atmdx = atmdx / 10.0
                continue                           # go to 1998
            break

        hmax = nml.atmmax
        hmin = nml.atmmin
        h = nml.atmbeg
        stoler = intpar_mod.stolr0
        nok = 0
        nbad = 0
        # The integration loop.
        for nstep in range(1, MAXSTP + 1):
            for i in range(1, nv + 1):
                yscal[i] = abs(y[i]) + abs(h * dydx[i]) + TINY
            # Don't pass the end, nor exceed the maximum step.
            if (x0 - xlim) * (x0 + h - xlim) < 0.0:
                h = xlim - x0
            if h > hmax:
                h = hmax
            # h is the attempted step, hdid the one taken, hnext the next.
            x0, hdid, hnext, rl, katm, ksaha = bsstep(
                y, dydx, nv, x0, h, eps, yscal, 0.0, 0.0, qatm, b, fpl, ftl, gl,
                latmo, lderiv, locond, lprt, rl, teffl, x, z, katm, ksaha, err)
            errsum[1] = errsum[1] + err[1]
            # dP/dtau at the start of the next step.
            katm, ksaha = qatm(x0, y, dydx, b, fpl, ftl, gl, latmo, lderiv,
                               locond, lprt, rl, teffl, x, z, katm, ksaha)
            if lprt and nml.lstatm:
                beta = 1.0 - c.ca3 * math.exp(c.cln * (4.0 * apr.at - apr.ap))
                chrh = 1.0 / q.qqdp
                cht = -chrh * q.qqdt
                cv = q.qqcp - math.exp(c.cln * (apr.ap - apr.ad - apr.at)) \
                    * cht**2 / chrh
                gam1 = chrh * q.qqcp / cv
                if not nml.lstch:
                    io.write(io.istor,
                             f" {apr.taul:11.7f}{apr.ap:11.7f}{apr.at:11.7f}"
                             f"{apr.ad:11.7f}{apr.ao:14.7E}"
                             f"{apr.afxion[1]:6.3f}{apr.afxion[2]:6.3f}"
                             f"{apr.afxion[3]:6.3f}{ksaha:7d}{katm:7d}"
                             f"{gam1:12.8f}{q.qqdp:12.8f}{q.qqdt:12.8f}"
                             f"{beta:12.8f}{q.qqcp:12.5E}{cv:12.5E}")  # 20
                # JvS: save the structure in atmstruct_mod.
                atm.atmop[nstep] = apr.ap
                atm.atmot[nstep] = apr.at
                atm.atmod[nstep] = apr.ad
                atm.abetas[nstep] = beta
                atm.agam1[nstep] = gam1
                atm.aqdt[nstep] = q.qqdt
                atm.afxions[1, nstep] = apr.afxion[1]
                atm.afxions[2, nstep] = apr.afxion[2]
                atm.afxions[3, nstep] = apr.afxion[3]
                atm.atmoo[nstep] = apr.ao
                atm.atmocp[nstep] = q.qqcp
            if hdid == h:
                nok = nok + 1
            else:
                nbad = nbad + 1
            # DBG pulse: atmosphere values for pulsation (JVS 02/11: also
            # for LCLCD).
            # Question (see the page): taup, op, dp are set only with
            # pulsation output; otherwise they are left from an earlier call.
            if (lpulpt and lprt) or acd.lclcd or nml.lstch:
                q.qqed = 0.0
                qesum = 0.0
                q.qqet = 0.0
                q.qfs = 1.0
                on = q.qo
                dn = math.exp(c.cln * q.qdl)
                taun = math.exp(c.cln * apr.taul)
                # See J.P. Cox, Princ. of Stell. Struc., p. 590.
                deltr = (taun - taup) / (((dn * on) + (dp * op)) / 2)
                pelpf = c.cgas * q.qt * q.qd * q.qemu
                op = on
                # From the first lines of TPGRAD.
                delr = q.qo * b * math.exp(c.cln * (q.qpl - hstot
                                                    - 4.0 * q.qtl + c.clsunl
                                                    - c.cgl + c.cdelrl)) \
                    * ftl / fpl
                del_ = delr if delr - q.qdela <= 1.0e-6 else q.qdela
                if lpulpt and lprt:
                    if nml.ipver == 1:
                        _pulse_lines(io.iopatm,
                                     [deltr, q.qfs, b, q.qtl, q.qdl, q.qpl,
                                      qesum, q.qo, q.qqdp, q.qqed, q.qqet,
                                      q.qqod, q.qqot, del_, q.qdela, q.qqcp,
                                      q.qrmu, q.qqdt, pelpf], 0)
                    elif nml.ipver in (2, 3):
                        _pulse_lines(io.iopatm,
                                     [deltr, q.qfs, b, q.qtl, q.qdl, q.qpl,
                                      taun, q.qo, q.qqdp, q.qqed, q.qqet,
                                      q.qqod, q.qqot, del_, q.qdela, q.qqcp,
                                      q.qrmu, q.qqdt, pelpf], 5)
                # JvS: save in atmstruct_mod.
                # Question (see the page): stored on other switches than
                # atmop, atmot, ...
                atm.atmor[nstep] = deltr
                atm.adels[1, nstep] = delr
                atm.adels[2, nstep] = del_
                atm.adels[3, nstep] = q.qdela
                atm.numatm = nstep
                taup = taun
                dp = dn
            # Is the integration complete?
            if abs(x0 - xlim) <= stoler:
                io.write(io.ishort,
                         f" ATMOSPHERE INTEGRATION COMPLETE NUMBER OF STEPS "
                         f"ACCEPTED{nok:5d} REJECTED{nbad:5d}\n"
                         f"     MAXIMUM RELATIVE ERROR IN P "
                         f"{errsum[1]:22.13E}")    # format 35
                break                              # go to 200
            if hnext < hmin:
                hnext = hmin
            h = hnext
        else:
            msg = ("     ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP "
                   "INTEGRATIONS.I QUIT.")         # format 50
            io.write(io.iowr, msg)
            io.write(io.ishort, msg)
            raise SystemExit(msg)

    # label 200 (tables come here directly). G Somers 3/17: if only the
    # photospheric pressure is wanted, stop here.
    ovrtrn_mod.pphot = apr.ap
    if not nml.lcalcenv:
        return rl, istore, katm, kenv, ksaha           # go to 555

    # The envelope: P is the independent variable, m, T and r the
    # dependent ones; integrate from tau = 2/3 to the model's last point.
    if lprt and nml.lstenv:
        if not nml.lstch:
            io.write(io.istor, "\n******** ENVELOPE BEGIN ********")  # 61
    # DBG pulse: end-of-data marker for the atmosphere file.
    if lpulpt and lprt:
        if nml.ipver == 1:
            _pulse_lines(io.iopatm, list(XYZ[1:20]), 0)
        elif nml.ipver in (2, 3):
            _pulse_lines(io.iopatm, list(XYZ[1:20]), 5)

    comp = comp_mod
    if comp.senv > -1.0e-12:
        # The envelope mass is negligible: the photosphere is the fitting
        # point.
        if lsbc:
            ps[ie] = apr.ap
            rs[ie] = rl
            ts[ie] = apr.at
            if lprt:
                if not nml.lstch:
                    io.write(io.istor,
                             f"    {ps[ie]:16.12f}{ts[ie]:16.12f}"
                             f"{rs[ie]:16.12f}        {comp.senv:16.12f}")
    else:
        # Initialize; set the numerical parameters.
        latmo = False
        for i in range(1, 4):
            errsum[i] = 0.0
        nok = 0
        nbad = 0
        x0 = apr.ap
        sl = 0.0
        y[1] = sl
        y[2] = apr.at
        y[3] = rl
        nv = 3
        eps = nml.enverr
        hmax = nml.envmax
        hmin = nml.envmin
        h = nml.envbeg
        stoler = abs(intpar_mod.stolr0 * comp.senv)
        rsurf = math.exp(c.cln * rl)
        if istore == ie:
            istore = 0
        llstore = False
        ltest = not lsbc
        # dy/dx at the start of the step.
        rl, kenv, ksaha = qenv(x0, y, dydx, b, fpl, ftl, gl, latmo,
                               lderiv, locond, lprt, rl, teffl, x, z, kenv,
                               ksaha)
        # DBG pulse: the first point of the envelope.
        pal = pualpha_mod
        if lpulpt and lprt:
            _envelope_pulse_point(rl, b, pal)
        if lprt and nml.lstenv:
            if not nml.lstch:
                io.write(io.istor,
                         " " + "".join(f"{s:>10}" for s in
                                       ("GRAV  ", "P   ", "T   "))
                         + "".join(f"{s:>14}" for s in ("DEPTH    ", "M      "))
                         + "".join(f"{s:>10}" for s in ("D    ", "O   "))
                         + f"{'BETA':>7}{'DELR  ':>9}"
                         + "".join(f"{s:>6}" for s in
                                   ("DELA", "DEL ", "HII ", "HEII", "HEIII"))
                         + f"{'V   ':>9}{'GAM1   ':>12}{'QQDP   ':>12}")  # 240
            srad = math.exp(c.cln * rl)
            de = (rsurf - srad) / rsurf
            sgl = c.cgl + hstot - 2.0 * rl
            eg = math.exp(c.cln * sgl)
            if not nml.lstch:
                io.write(io.istor, _fmt260(eg, de))
        # Store the starting values (07/02: one stored envelope point).
        lsurcz = False
        tauczenv = 0.0
        _store_envelope_point(1, x, z, b)
        # JVS 10/10
        chdelj = epr.edel[2]                       # set, never used
        chdeld = q.qdela                           # set, never used
        if env.lcenv[1]:
            lsurcz = True
        env.numenv = 1
        for nstep in range(1, MAXSTP + 1):
            xstart = x0
            for i in range(1, nv + 1):
                ystart[i] = y[i]
                yscal[i] = abs(y[i]) + abs(h * dydx[i]) + TINY
            sest = y[1] + h * dydx[1]
            if comp.senv - y[1] > 0.0 or comp.senv - sest > 0.0:
                # The integration has overshot the fitting point, or the
                # next step would: shorten the step, or integrate back.
                h = (comp.senv - y[1]) / dydx[1]
            if h < 0.0:
                if h < -hmax:
                    h = -hmax
            else:
                if h > hmax:
                    h = hmax
            # plim estimates the end pressure from the model's last point:
            # the first time the integrator would pass it, limit the step.
            if ltest:
                if x0 + h > plim:
                    h = plim - x0
                    ltest = False
            x0, hdid, hnext, rl, kenv, ksaha = bsstep(
                y, dydx, nv, x0, h, eps, yscal, 0.0, 0.0, qenv, b, fpl, ftl, gl,
                latmo, lderiv, locond, lprt, rl, teffl, x, z, kenv, ksaha, err)
            for k in range(1, 4):
                errsum[k] = errsum[k] + err[k]
            # dy/dx at the start of the next step.
            rl, kenv, ksaha = qenv(x0, y, dydx, b, fpl, ftl, gl, latmo,
                                   lderiv, locond, lprt, rl, teffl, x, z,
                                   kenv, ksaha)
            if lpulpt and lprt:
                _envelope_pulse_point(rl, b, pal)
            if lprt and nml.lstenv:
                srad = math.exp(c.cln * rl)
                de = (rsurf - srad) / rsurf
                sgl = c.cgl + hstot - 2.0 * rl
                eg = math.exp(c.cln * sgl)
                if not nml.lstch:
                    io.write(io.istor, _fmt260(eg, de))
            if hdid == h:
                nok = nok + 1
            else:
                nbad = nbad + 1
            # Is the integration complete? Store the point if it has not
            # overshot the fitting point (07/02).
            sdif = comp.senv - y[1]
            if sdif <= stoler:
                env.numenv = env.numenv + 1
                n = env.numenv
                _store_envelope_point(n, x, z, b)
                if not lsurcz:
                    if env.lcenv[n]:
                        lsurcz = True
                elif epr.evel > 0.0:
                    dr = 10.0**env.envr[n - 1] - 10.0**env.envr[n]
                    # Question (see the page): tauczenv is never kept.
                    tauczenv = tauczenv + dr / epr.evel
            if abs(sdif) <= stoler:
                if lsbc:
                    # Interpolate linearly to the fitting mass.
                    wt0 = sdif / (ystart[1] - y[1])
                    ps[ie] = x0 + wt0 * (xstart - x0)
                    rs[ie] = y[3] + wt0 * (ystart[3] - y[3])
                    ts[ie] = y[2] + wt0 * (ystart[2] - y[2])
                    if lprt:
                        io.write(io.ishort,
                                 f"    {ps[ie]:16.12f}{ts[ie]:16.12f}"
                                 f"{rs[ie]:16.12f}        {comp.senv:16.12f}")
                break                              # go to 300
            elif not llstore:
                # TENV was COMMON/CTLIM/, set by parmin.f; derived here.
                tenv = 0.5 * (nml.tenv0 + nml.tenv1)
                if y[2] >= tenv and lsbc:
                    llstore = True
                    istore = ie
                    estore[1] = x0
                    estore[2] = y[2]
                    estore[3] = y[3]
                    estore[4] = y[1]
            if hnext < hmin:
                hnext = hmin
            h = hnext
        else:
            msg = ("     ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES. "
                   "I QUIT")                       # format 911
            io.write(io.iowr, msg)
            io.write(io.ishort, msg)
            raise SystemExit(msg)

    # label 300 (07/02): reverse the envelope profile, centre outward.
    if comp.senv < -1.0e-12:
        n = env.numenv
        for i in range(1, n + 1):
            i1 = i
            i2 = n - i + 1
            if i1 >= i2:
                break                              # go to 310
            for arr in (env.envd, env.envp, env.envr, env.envs, env.envt,
                        env.envx, env.envz, env.lcenv, env.envo):
                arr[i1], arr[i2] = arr[i2], arr[i1]
            # Question (see the page): ENVL(I2) = DUM1 is missing, so the
            # second half keeps its own values.
            dum1 = env.envl[i1]
            env.envl[i1] = env.envl[i2]
            for arr in (env.eqdt, env.evels, env.ebetas, env.egam1,
                        env.eqcp):
                arr[i1], arr[i2] = arr[i2], arr[i1]
            for k in range(1, 4):
                env.edels[k, i1], env.edels[k, i2] = \
                    env.edels[k, i2], env.edels[k, i1]
                env.efxions[k, i1], env.efxions[k, i2] = \
                    env.efxions[k, i2], env.efxions[k, i1]
    # DBG pulse: end-of-data marker for the envelope file.
    if lpulpt and lprt:
        if nml.ipver == 1:
            _pulse_lines(io.iopenv, list(XYZ[1:20]), 0)
        elif nml.ipver in (2, 3):
            _pulse_lines(io.iopenv, list(XYZ[1:23]), 6)

    # Possible bug (update-findings.md): when the envelope is skipped,
    # nok, nbad, sdif and y are from the atmosphere or an earlier call.
    io.write(io.ishort,
             f" ENVELOPE INTEGRATION COMPLETE NUMBER OF STEPS ACCEPTED"
             f"{nok:5d} REJECTED{nbad:5d}\n"
             f"     SENV-LAST M={sdif:22.13E}  LAST M={y[1]:22.13E}\n"
             f"     MAX RELATIVE ERRORS:M {errsum[1]:14.5E}  T "
             f"{errsum[2]:14.5E}  R {errsum[3]:14.5E}")   # format 215

    # JVS 08/13: if the convection zone's base lies in the envelope, find
    # it, and the overturn time (G Somers 3/17: not with LNEWTCZ).
    if nml.lnewtcz:
        return rl, istore, katm, kenv, ksaha           # go to 555
    jenv = 0
    for i in range(1, env.numenv + 1):
        if jenv == 0 and env.lcenv[i]:
            jenv = i
    if jenv > 1:
        # The base of the surface convection zone in radius, by linear
        # interpolation of del_rad - del.
        dd2 = env.edels[1, jenv - 1] - env.edels[2, jenv - 1]
        dd1 = env.edels[1, jenv] - env.edels[2, jenv]
        fx = dd2 / (dd2 - dd1)
        # Question (see the page): envczl is no longer used.
        envczl = (env.envr[jenv - 1] + fx * (env.envr[jenv]
                                             - env.envr[jenv - 1])
                  - c.crsunl)
        for i in range(1, env.numenv + 1):
            envs1[i] = math.exp(c.cln * (env.envs[i] + hstot))
            if i == env.numenv:
                envs2[i] = 0.0
            else:
                envs2[i] = (math.exp(c.cln * env.envs[i + 1])
                            - math.exp(c.cln * env.envs[i]))
            edel1[i] = env.edels[1, i]
            edel2[i] = env.edels[2, i]
        # Question (see the page): envg is never set.
        taucal(envs2, envs1, env.lcenv, env.envr, env.envp, env.envd, envg,
               env.numenv, env.evels, edel1, edel2)
    return rl, istore, katm, kenv, ksaha               # label 555


def _store_envelope_point(n, x, z, b):
    """Copy the point just integrated (envprt_mod, pulse2_mod) into entry n
    of envstruct_mod: the Fortran writes this block out twice."""
    env = envstruct_mod
    e = envprt_mod
    q = pulse2_mod
    env.envd[n] = e.ed
    env.envp[n] = e.ep
    env.envr[n] = e.er
    env.envs[n] = e.es
    env.envt[n] = e.et
    env.envx[n] = x
    env.envz[n] = z
    env.lcenv[n] = e.evel > 0.0
    for k in range(1, 4):
        env.edels[k, n] = e.edel[k]
    env.evels[n] = e.evel
    env.ebetas[n] = e.ebeta
    # JVS 08/25: always save these.
    chrh = 1.0 / q.qqdp
    cht = -chrh * q.qqdt
    cv = q.qqcp - math.exp(constants_mod.cln * (e.ep - e.ed - e.et)) \
        * cht**2 / chrh
    env.egam1[n] = chrh * q.qqcp / cv
    env.eqcp[n] = q.qqcp
    for k in range(1, 4):
        env.efxions[k, n] = e.efxion[k]
    env.envo[n] = e.eo
    env.envl[n] = b
    env.eqdt[n] = q.qqdt


def _envelope_pulse_point(rl, b, pal):
    """DBG pulse: one envelope point to the pulsation file (FORMATs 5001,
    6001, or 6003 with the mixing-length values appended to the first
    three lines); the Fortran writes this block out twice."""
    q = pulse2_mod
    io = iounits_mod
    q.qqed = 0.0
    qesum = 0.0
    q.qqet = 0.0
    pelpf = constants_mod.cgas * q.qt * q.qd * q.qemu
    vals = [rl, q.qfs, b, q.qtl, q.qdl, q.qpl, qesum, q.qo, q.qqdp, q.qqed,
            q.qqet, q.qqod, q.qqot, q.qdel, q.qdela, q.qqcp, q.qrmu, q.qqdt,
            pelpf]
    if nmlpar_mod.ipver == 1:
        _pulse_lines(io.iopenv, vals, 0)
    elif nmlpar_mod.ipver == 2:
        _pulse_lines(io.iopenv, vals, 5)
    elif nmlpar_mod.ipver == 3:
        _pulse_lines(io.iopenv,
                     vals[0:5] + [pal.alfmlt] + vals[5:10] + [pal.phmlt]
                     + vals[10:15] + [pal.cmxmlt] + vals[15:19], 6)
