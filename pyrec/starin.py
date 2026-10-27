"""Python version of starin.f (YREC 09c3b9b): the starting model of a run,
read from file or kept, and prepared for evolution."""

import math

import numpy as np

from pyrec import cenv_mod
from pyrec import comp2_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import dimensions_mod
from pyrec import dpmix_mod
from pyrec import envstruct_mod
from pyrec import iounits_mod
from pyrec import newopac_mod
from pyrec import nmlbypass_mod  # noqa: F401 (USEd; no member read here)
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldrot_mod
from pyrec import ovrtrn_mod
from pyrec import runloop_mod
from pyrec.envint import envint
from pyrec.eqstat import eqstat
from pyrec.fpft import fpft
from pyrec.getmodel2 import getmodel2
from pyrec.getopac import getopac
from pyrec.gettau import gettau
from pyrec.getyrec7 import getyrec7
from pyrec.momi import momi
from pyrec.ovrot import ovrot
from pyrec.physic import physic
from pyrec.rscale import rscale
from pyrec.setscv import setscv
from pyrec.surfopac import surfopac
from pyrec.tpgrad import tpgrad

JSON = dimensions_mod.json

# DATA ATOMWT: atomic weights of the 12 elements of V (index 0 unused).
ATOMWT = np.array([0.0, 23.0, 26.99, 24.32, 55.86, 28.1, 12.015, 1.008,
                   16.0, 14.01, 39.96, 20.19, 4.004])

# Local arrays kept between calls (SAVE).
mxzone = np.zeros((13, 3), dtype=int)
mrzone = np.zeros((14, 3), dtype=int)
lcz = np.zeros(JSON + 1, dtype=bool)


def _set_debhu(hcomp, m):
    # The composition for the Debye-Hueckel correction, at shell m.
    debhu_mod.xxdh = hcomp[1, m]
    debhu_mod.yydh = hcomp[2, m] + hcomp[4, m]
    debhu_mod.zzdh = hcomp[3, m]
    debhu_mod.zdh[1] = hcomp[5, m] + hcomp[6, m]
    debhu_mod.zdh[2] = hcomp[7, m] + hcomp[8, m]
    debhu_mod.zdh[3] = hcomp[9, m] + hcomp[10, m] + hcomp[11, m]


def _e(x, w, d):
    # Fortran's 0PEw.d: a mantissa in [0.1, 1).
    if x == 0.0:
        return f"{'0.' + '0' * d + 'E+00':>{w}}"
    e = math.floor(math.log10(abs(x))) + 1
    mant = x / 10.0**e
    if round(abs(mant), d) >= 1.0:
        mant /= 10.0
        e += 1
    s = ("-" if x < 0 else "") + "0" + f"{abs(mant):.{d}f}"[1:]
    return f"{s + f'E{e:+03d}':>{w}}"


def starin(bl, cfenv, dage, ddage, delts, deltsh, eta2, fp, ft, ftri, hcomp,
           hd, hi, hjm, hkerot, hl, hp, hr, hs, hs1, hs2, hstot, ht, ikut,
           istore, jenv, large, lc, lnew, m, model, nk, omega, ps, qdp, qdt,
           qiw, r0, rs, sjtot, skerot, smass, teffl, tlumx, tril, trit, ts,
           vel, hg, v, cmixl, compmix):
    """Provide the starting model of run nk: read it from the unit IFIRST
    (YREC7 or MODEL2 format) or keep the one in memory, then prepare it:
    change the mixture of a homogeneous model, extend the core inward,
    rescale, move the envelope's fitting point, compute the shell masses,
    the rotation quantities, the surface composition, and the first
    physics, convection zones and overturn time.

    Translates: starin.f#STARIN (09c3b9b)

    Arguments: the model (hs, hr, hl, hp, ht, hd, hcomp[species, shell],
    lc, omega, m points, index 0 unused) and its rotation arrays (eta2, fp,
    ft, hi, hjm, hkerot, qiw, r0, hg), the envelope triangle (trit, tril,
    ps, ts, rs, cfenv, ftri, istore), bl, teffl, hstot, smass, dage, ddage,
    model, tlumx, the run number nk, jenv (read only, for the spots), large
    (the previous model failed: stop), the totals sjtot, skerot, v (the
    surface mixture's relative mass fractions), cmixl (the run's mixing
    length) and compmix (set by the readers). delts, deltsh, ikut, lnew,
    qdp, qdt and vel are output only. Arrays are changed in place.

    Returns (bl, dage, ddage, delts, deltsh, ftri, hstot, ikut, istore,
    lnew, m, model, qdp, qdt, sjtot, skerot, smass, teffl, vel, compmix).
    """
    c = constants_mod
    nml = nmlpar_mod
    es = envstruct_mod
    ishort = iounits_mod.ishort

    if large:
        iounits_mod.write(ishort, " " + ">" * 39 + "<" * 40 + "\n"
                          "STARIN:        ***** RUN STOPPED *****")
        iounits_mod.write(ishort,
                          "STARIN: ***** MODEL FAILED TO CONVERGE *****")
        raise SystemExit()

    ikut = 0
    iread = iounits_mod.ifirst
    vel = 0.0
    qdt = -1.0
    qdp = 1.0

    if not nml.lfirst[nk]:
        # Keep the model in memory (rescaled if asked).
        if runloop_mod.iresca[nk] != 1:
            hstot, smass = rscale(hl, hcomp, hs, hstot, m, nk, smass, lc)
    else:
        # Read the starting model; its first four characters give the
        # format.
        iounits_mod.rewind(iread)
        atemp = iounits_mod.read(iread)[:4]  # format 10, A4
        # The readers' returns (not yet translated) are taken as every
        # scalar argument except the units, in argument order.
        if atemp == "NMOD":
            iounits_mod.write(ishort, "STARIN:  Input model has YREC7 format")
            (bl, cmixl0, dage, ddage, ftri, hstot, jcore0, jenv0, lexcom0,
             lrot0, m, model, smass, teffl, atm, eos, hik, ldify0, ldifz0,
             ldisk0, linstb0, ljdot00, lok, lovstc0, lovste0, lovstm0,
             lpurez0, lsemic0, compmix, pdisk0, tdisk0, wmax0) = getyrec7(
                bl, cfenv, 0.0, dage, ddage, ftri, hcomp, hd, hl, hp, hr, hs,
                hstot, ht, iread, ishort, 0, 0, lc, False, False, m, model,
                omega, ps, rs, smass, teffl, tlumx, tril, trit, ts, "", "",
                "", False, False, False, False, False, "", False, False,
                False, False, False, compmix, 0.0, 0.0, 0.0)
        elif atemp == "MOD2":
            iounits_mod.write(ishort,
                              "STARIN:  Input model has MODEL2 format")
            (bl, cmixl0, dage, ddage, ftri, hstot, jcore0, jenv0, lexcom0,
             lrot0, m, model, smass, teffl, atm, eos, hik, ldify0, ldifz0,
             ldisk0, linstb0, ljdot00, lok, lovstc0, lovste0, lovstm0,
             lpurez0, lsemic0, compmix, pdisk0, tdisk0, wmax0) = getmodel2(
                bl, cfenv, 0.0, dage, ddage, ftri, hcomp, hd, hl, hp, hr, hs,
                hstot, ht, iread, 0, 0, lc, False, False, m, model, omega,
                ps, rs, smass, teffl, tlumx, tril, trit, ts, "", "", "",
                False, False, False, False, False, "", False, False, False,
                False, False, compmix, 0.0, 0.0, 0.0)
        else:
            iounits_mod.write(ishort, "STARIN: ***** RUN TERMINATED, "
                              "INVALID INPUT MODEL FILE.  *****")
            raise SystemExit()

        # The step in seconds (ddage in years).
        delts = c.csecyr * ddage
        deltsh = abs(delts)
        comp_mod.stotal = hstot

        # An evolved model: warn if its mixing length or composition switch
        # differs from the run's (the run is no longer stopped).
        # Question (see the page): should a mismatch still stop the run?
        if model > 0:
            lmixl = True
            if cmixl0 > 0.0:
                lmixl = abs(cmixl - cmixl0) < 2.0e-3
            if not lmixl or lexcom0 != nml.lexcom:
                # Question (see the page): LEXCP0 is never set; the model's
                # value is LEXCOM0. False stands for the undefined value.
                lexcp0 = False
                text = (" ERROR IN SUBROUTINE STARIN\n USER PARAMETERS"
                        " OF WRONG TYPE FOR INITIAL MODEL\n MIXING LENGTH -"
                        f" USER DESIRES{cmixl:7.3f} MODEL MIX LENGTH"
                        f"{cmixl0:7.3f}\n EXTENDED COMP-USER DESIRES "
                        f"{'T' if nml.lexcom else 'F'} MODEL USED "
                        f"{'T' if lexcp0 else 'F'}")  # format 1040
                iounits_mod.write(ishort, text)
                iounits_mod.write(iounits_mod.iowr, text)

        istore = 0
        ftri = 1.0
        lnew = True  # the envelope is always recomputed
        comp_mod.xnew = hcomp[1, m]
        comp_mod.znew = hcomp[3, m]

        if nml.lrot:
            if nml.lwnew:
                for i in range(1, m + 1):  # loop 540
                    omega[i] = nml.wnew
        else:
            for i in range(1, JSON + 1):  # loop 570
                fp[i] = 1.0
                ft[i] = 1.0
        iounits_mod.rewind(iread)

        # MHP 10/24: change the mixture, only of a homogeneous model.
        testmix = np.zeros(16)
        homogeneous = True
        if nml.lmixture or nml.lisotope:
            for i in range(1, 16):
                testmix[i] = hcomp[i, 1]
            for j in range(2, m + 1):
                for i in range(1, 16):
                    test = abs(hcomp[i, j] - testmix[i])
                    if test > 1.0e-6:
                        text = (f"SPECIES {i:3d} IN SHELL {j:5d} DIFFERS FROM"
                                f" CENTER BY {_e(test, 12, 4)} MIX NOT"
                                f" MODIFIED IN EVOLVED MODEL")  # format 592
                        print(text)
                        iounits_mod.write(ishort, text)
                        homogeneous = False
                        break
                if not homogeneous:
                    break  # GOTO 602
        if homogeneous:
            if nml.lmixture:
                # Scale C, N and O so that each element's share of Z is
                # FRAC_C, FRAC_N, FRAC_O.
                tot_c = (testmix[5] + testmix[6]) / testmix[3]
                tot_n = (testmix[7] + testmix[8]) / testmix[3]
                tot_o = (testmix[9] + testmix[10] + testmix[11]) / testmix[3]
                rat_c = nml.frac_c / tot_c
                rat_n = nml.frac_n / tot_n
                rat_o = nml.frac_o / tot_o
                print(nml.frac_c, nml.frac_n, nml.frac_o)
                for i in range(5, 7):
                    for j in range(1, m + 1):
                        hcomp[i, j] = rat_c * hcomp[i, j]
                for i in range(7, 9):
                    for j in range(1, m + 1):
                        hcomp[i, j] = rat_n * hcomp[i, j]
                for i in range(9, 12):
                    for j in range(1, m + 1):
                        hcomp[i, j] = rat_o * hcomp[i, j]
                text = ("CNO MIX CHANGED IN STARIN. OLD C12 C13 N14"
                        " N15 O16 O17 O18 "
                        + "".join(_e(testmix[k], 12, 4) for k in range(5, 12))
                        + " NEW "
                        + "".join(_e(hcomp[k, 1], 12, 4)
                                  for k in range(5, 12)))  # format 594
                print(text)
                iounits_mod.write(ishort, text)
            if nml.lisotope:
                # Impose 12C/13C and 16O/18O at fixed sums, and the light
                # elements' initial values.
                # Question (see the page): R14_15, R16_17 are not applied.
                sum_c = hcomp[5, 1] + hcomp[6, 1]
                sum_o = hcomp[9, 1] + hcomp[11, 1]
                for j in range(1, m + 1):
                    hcomp[4, j] = nml.xhe3_ini
                    hcomp[5, j] = nml.r12_13 * sum_c / (1.0 + nml.r12_13)
                    hcomp[6, j] = sum_c / (1.0 + nml.r12_13)
                    hcomp[9, j] = nml.r16_18 * sum_o / (1.0 + nml.r16_18)
                    hcomp[11, j] = sum_o / (1.0 + nml.r16_18)
                    hcomp[12, j] = nml.xh2_ini
                    hcomp[13, j] = nml.xli6_ini
                    hcomp[14, j] = nml.xli7_ini
                    hcomp[15, j] = nml.xbe9_ini
                text = ("CNO ISOTOPES AND LIGHT ELEMENTS CHANGED IN "
                        "STARIN. OLD HE3 C12 C13 N14 N15 O16 O17 O18 H2 LI6 "
                        "LI7 BE9"
                        + "".join(_e(testmix[k], 12, 4) for k in range(4, 16))
                        + " NEW "
                        + "".join(_e(hcomp[k, 1], 12, 4)
                                  for k in range(4, 16)))  # format 593
                print(text)
                iounits_mod.write(ishort, text)
        # label 602

    # label 3000
    # MHP 9/98, 9/14: extend the model inward by MCORE points.
    if nml.lcore and nml.lfirst[nk]:
        nml.mcore = int(math.log10(nml.fcore) / nml.hpttol[2]) + 1
        nml.fcore = float(nml.mcore) * nml.hpttol[2]  # now a log offset
        mn = m + nml.mcore
        if mn > JSON:
            # Question (see the page): no stop; the arrays overflow below.
            iounits_mod.write(ishort, "STARIN: Unable to extend core inward "
                              "- JSON too small")
            iounits_mod.write(ishort, f"STARIN: Required size ={mn:8d}"
                              f", JSON = {JSON:8d}")
            iounits_mod.write(ishort, "STARIN:  ***** RUN TERMINATED *****")
        ds = nml.hpttol[2]
        mcore = nml.mcore
        # Shift the old points outward.
        # Question (see the page): hd is not shifted.
        for i in range(m, 0, -1):
            hs[i + mcore] = hs[i]
            hr[i + mcore] = hr[i]
            hl[i + mcore] = hl[i]
            hp[i + mcore] = hp[i]
            ht[i + mcore] = ht[i]
            lc[i + mcore] = lc[i]
            for j in range(1, 16):
                hcomp[j, i + mcore] = hcomp[j, i]
            omega[i + mcore] = omega[i]
        mp1 = mcore + 1
        rhoc = hd[1]
        # Question (see the page): with mcore = 1 this is 0/0.
        dele = (ht[2] - ht[1]) / (hp[2] - hp[1])
        hlc = hl[1]
        facd = hp[mp1] - hd[mp1] - ht[mp1]
        # Inward at constant density, energy generation and gradient.
        for i in range(mcore, 0, -1):
            hs[i] = hs[i + 1] - ds
            hr[i] = c.cc13 * (hs[i] - c.c4pi3l - rhoc)
            hl[i] = math.exp(c.cln * (hs[i] - hs[mp1])) * hlc
            pp = hp[i + 1]
            temp = math.exp(c.cln * (c.cgl + 2.0 * hs[i] - c.c4pil - pp
                                     - 4.0 * hr[i]))
            pp = pp + 0.5 * temp * ds
            temp = math.exp(c.cln * (c.cgl + 2.0 * hs[i] - c.c4pil - pp
                                     - 4.0 * hr[i]))
            hp[i] = hp[i + 1] + temp * ds
            lc[i] = lc[mp1]
            ht[i] = ht[i + 1] + temp * dele * ds
            omega[i] = omega[mp1]
            for j in range(1, 16):
                hcomp[j, i] = hcomp[j, mp1]
            hd[i] = hp[i] - ht[i] - facd
        m = mn

    if runloop_mod.iresca[nk] != 1:
        hstot, smass = rscale(hl, hcomp, hs, hstot, m, nk, smass, lc)

    # Move the envelope's fitting point to SENV0.
    if cenv_mod.lsenv0:
        # (In one branch EQSTAT and TPGRAD set STARIN's qdt, qdp and vel.)
        m, lnew, qdt, qdp, vel = _new_fitting_point(
            hcomp, hd, hjm, hl, hp, hr, hs, ht, lc, omega, m, bl, teffl,
            hstot, jenv, lnew, qdt, qdp, vel, cmixl)

    # Shell masses.
    ds3 = math.exp(c.cln * hs[1])
    ds2 = -ds3
    for i in range(2, m + 1):  # loop 120
        ds1 = ds2
        ds2 = ds3
        ds3 = math.exp(c.cln * hs[i])
        hs1[i - 1] = ds2
        hs2[i - 1] = 0.5 * (ds3 - ds1)
    hs1[m] = ds3
    hs2[m] = math.exp(c.cln * hstot) - 0.5 * (ds2 + ds3)

    # Rotation: distortion, moments of inertia, J and kinetic energy.
    if nml.lrot:
        fpft(hd, hr, hs, m, omega, eta2, fp, ft, hg, r0)
        momi(eta2, hr, hs, hs2, 1, m, omega, r0, hi, qiw)
        sumj = 0.0
        sumke = 0.0
        for i in range(1, m + 1):  # loop 550
            hj = omega[i] * hi[i]
            hjm[i] = hj / hs2[i]
            hkerot[i] = 0.5 * omega[i] * hj
            sumj = sumj + hj
            sumke = sumke + hkerot[i]
        iounits_mod.write(ishort,
                          f" TOTAL J OF STAR - PREVIOUS {sjtot:21.13E} NEW "
                          f"{sumj:21.13E}\n TOTAL ROTATIONAL K.E. OF STAR - "
                          f"PREVIOUS {skerot:21.13E} NEW {sumke:21.13E}")
        sjtot = sumj
        skerot = sumke

    # The surface composition, for the first run only.
    # Question (see the page): skipped for later runs.
    if nk <= 1:
        if not newopac_mod.llaol:
            for i in range(1, 13):
                v[i] = nml.vnew[i]
        comp_mod.xenv = comp_mod.xnew
        comp_mod.zenv = comp_mod.znew
        comp2_mod.yenv = 1.0 - comp_mod.xenv - comp_mod.zenv - hcomp[4, m]
        comp2_mod.y3env = hcomp[4, m]
        vsum = v[1] + v[2] + v[3] + v[4] + v[5] + v[6] + v[8] + v[9] + v[10] \
            + v[11]
        comp_mod.zenvm = comp_mod.zenv * (vsum - v[6] - v[8] - v[9]) / vsum
        fv = comp_mod.zenv / vsum
        v[7] = comp_mod.xenv / fv
        v[12] = (1.0 - comp_mod.xenv - comp_mod.zenv) / fv
        vsum = 0.0
        for i in range(1, 13):  # loop 610
            v[i] = fv * v[i] / ATOMWT[i]
            vsum = vsum + v[i]
        # Question (see the page): the comment calls this the mean atomic
        # weight of the surface Z.
        comp_mod.amuenv = vsum
        fv = 1.0 / comp_mod.amuenv
        for i in range(1, 13):  # loop 620
            comp_mod.fxenv[i] = v[i] * fv
    # label 630
    surfopac(comp_mod.xenv)
    if nml.lscv:
        setscv()

    # Save the old model.
    om = oldmod_mod
    for i in range(1, m + 1):  # loop 710
        om.hpo[i] = hp[i]
        om.hto[i] = ht[i]
        om.hro[i] = hr[i]
        om.hlo[i] = hl[i]
        # Question (see the page): these two are loop-invariant.
        om.tefflo = teffl
        om.mo = m
    if nml.lrot:
        for i in range(1, m + 1):  # loop 720
            oldrot_mod.wold[i] = omega[i]

    # The first physics, convection zones and overturn time.
    physic(fp, ft, hcomp, hd, hg, hl, hp, hr, hs, ht, lc, m, teffl, cmixl)
    nrzone, nzone = ovrot(hcomp, hd, hp, hr, hs, ht, lc, m, lcz, mrzone,
                          mxzone, 0, 0)
    # Possible bug (update-findings, batch 12): the module array ENVR is
    # passed as GETTAU's scalar RBCZ, so ENVR(1) receives the base's radius.
    es.envr[1] = gettau(hcomp, hr, hp, hd, hs1, ht, fp, ft, teffl, hstot,
                        bl, m, lc, es.envr[1], cmixl)
    ovrtrn_mod.taucz0 = ovrtrn_mod.taucz
    ovrtrn_mod.pphot0 = ovrtrn_mod.pphot
    ovrtrn_mod.fracstep = 0.5

    return (bl, dage, ddage, delts, deltsh, ftri, hstot, ikut, istore, lnew,
            m, model, qdp, qdt, sjtot, skerot, smass, teffl, vel, compmix)


def _new_fitting_point(hcomp, hd, hjm, hl, hp, hr, hs, ht, lc, omega, m, bl,
                       teffl, hstot, jenv, lnew, qdt, qdp, vel, cmixl):
    """STARIN's block for LSENV0 (labels 575-599): move the envelope's
    fitting point to SENV0, deeper (drop shells, interpolate the new last
    one) or shallower (integrate an envelope and take its points down to
    the new fitting mass as interior shells). Returns (m, lnew, qdt, qdp,
    vel): STARIN's arguments, which the physics of one branch sets."""
    c = constants_mod
    nml = nmlpar_mod
    ce = cenv_mod
    es = envstruct_mod
    ishort = iounits_mod.ishort

    if ce.senv0 > 0.0:
        ce.senv0 = -ce.senv0
    if ce.senv0 > -1.0e-12:
        ce.senv0 = -1.0e-12
    comp_mod.senv = hs[m] - hstot
    senv1 = comp_mod.senv
    if comp_mod.senv == ce.senv0:
        return m, lnew, qdt, qdp, vel  # GOTO 599
    jend = 11
    if nml.lexcom:
        jend = 15

    if ce.senv0 < comp_mod.senv:
        # A deeper fitting point: drop the shells beyond it.
        hsend = hstot + ce.senv0
        i = 0
        for ii in range(m - 1, 0, -1):  # loop 575
            if hs[ii] < hsend:
                i = ii
                break  # to label 580
        else:
            iounits_mod.write(ishort,
                              f"     ERROR IN SUBROUTINE STARIN\n     DESIRED"
                              f" ENVELOPE MASS{ce.senv0:22.13E} TOO LARGE\n"
                              f"     ENVELOPE MASS NOT CHANGED")  # format 576
            return m, lnew, qdt, qdp, vel  # GOTO 599
        # label 580
        m = i + 1
        comp_mod.senv = ce.senv0
        fs = (hsend - hs[i]) / (hs[i + 1] - hs[i])
        hs[m] = hsend
        hd[m] = hd[i] + fs * (hd[i + 1] - hd[i])
        hl[m] = hl[i] + fs * (hl[i + 1] - hl[i])
        hp[m] = hp[i] + fs * (hp[i + 1] - hp[i])
        hr[m] = hr[i] + fs * (hr[i + 1] - hr[i])
        ht[m] = ht[i] + fs * (ht[i + 1] - ht[i])
        for j in range(1, jend + 1):  # loop 585
            hcomp[j, m] = hcomp[j, i]
        comp_mod.xnew = hcomp[1, m]
        comp_mod.znew = hcomp[3, m]
        if nml.lrot:
            omega[m] = omega[i] + fs * (omega[i + 1] - omega[i])
        if lc[i] and lc[i + 1]:
            lc[m] = True
        elif not lc[i] and not lc[i + 1]:
            lc[m] = False
        else:
            # Between a convective and a radiative shell: let the physics
            # decide. (Since 09c3b9b EQSTAT is always called; it used to be
            # called only with the MHD equation of state.)
            lderiv = False
            latmo = True
            ksaha = 0
            x = hcomp[1, m]
            z = hcomp[3, m]
            pl = hp[m]
            tl = ht[m]
            dl = hd[m]
            b = hl[m]
            rl = hr[m]
            sl = hs[m]
            fpl = 1.0
            ftl = 1.0
            if nml.ldh:
                _set_debhu(hcomp, m)
            fxion = np.zeros(4)
            (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp,
             qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
                tl, 0.0, pl, 0.0, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, lderiv, latmo, ksaha)
            o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
            dpmix_mod.iovim = -1
            (dele, delr, qact, qacp, qacr, vel, lconv) = tpgrad(
                tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, 0.0,
                0.0, dela, qdtt, qdtp, qat, qap, 0.0, 0.0, 0.0, qcpt, qcpp,
                vel, lderiv, False, fpl, ftl, teffl, cmixl)
            hd[m] = dl
            lc[m] = lconv
    else:
        # A shallower fitting point: integrate an envelope from the surface
        # at fixed steps HPTTOL(8), and take its points down to SENV0.
        mm = m
        emax0 = nml.envmax
        emin0 = nml.envmin
        ebeg0 = nml.envbeg
        nml.envmax = nml.hpttol[8]
        nml.envmin = nml.hpttol[8]
        nml.envbeg = nml.hpttol[8]
        lsbc0 = False
        lprt = True
        katm = 0
        kenv = 0
        ksaha = 0
        b = math.exp(c.cln * bl)
        rl = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
        gl = c.cgl + comp_mod.stotal - rl - rl
        x = hcomp[1, m]
        z = hcomp[3, m]
        fpl = 1.0
        ftl = 1.0
        ixx = 0
        plim = hp[m]
        lpulpt = False
        if nml.ldh:
            _set_debhu(hcomp, m)
        idum = 0
        # G Somers 10/14: for spotted stars, the ambient temperature.
        if jenv == m and nml.spotf != 0.0 and nml.spotx != 1.0:
            ateffl = teffl - 0.25 * math.log10(
                nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
        else:
            ateffl = teffl
        dum1 = np.zeros(5)
        dum2 = np.zeros(4)
        dum3 = np.zeros(4)
        dum4 = np.zeros(4)
        rl, idum, katm, kenv, ksaha = envint(
            b, fpl, ftl, gl, hstot, ixx, lprt, lsbc0, plim, rl, ateffl, x, z,
            dum1, idum, katm, kenv, ksaha, dum2, dum3, dum4, lpulpt, cmixl)
        nml.envmax = emax0
        nml.envmin = emin0
        nml.envbeg = ebeg0
        comp_mod.senv = ce.senv0
        if m + es.numenv >= JSON:
            raise SystemExit(9999)
        # Shift the envelope to match the last interior point, dropping its
        # first point.
        dp = hp[m] - es.envp[1]
        dd = hd[m] - es.envd[1]
        dt = ht[m] - es.envt[1]
        dr = hr[m] - es.envr[1]
        for j in range(1, es.numenv):
            es.envd[j] = es.envd[j + 1] + dd
            es.envp[j] = es.envp[j + 1] + dp
            es.envr[j] = es.envr[j + 1] + dr
            es.envs[j] = es.envs[j + 1]
            es.envt[j] = es.envt[j + 1] + dt
            es.envx[j] = es.envx[j + 1]
            es.envz[j] = es.envz[j + 1]
        es.numenv = es.numenv - 1
        senv = comp_mod.senv
        stotal = comp_mod.stotal
        reached = False
        for j in range(m + 1, m + es.numenv + 1):
            jm = j - m
            hl[j] = hl[m]
            if es.envs[jm] <= senv:
                hd[j] = es.envd[jm]
                hp[j] = es.envp[jm]
                hr[j] = es.envr[jm]
                hs[j] = es.envs[jm] + stotal
                ht[j] = es.envt[jm]
                hcomp[1, j] = es.envx[jm]
                hcomp[3, j] = es.envz[jm]
                for k in range(4, jend + 1):
                    hcomp[k, j] = hcomp[k, m]
                hcomp[2, j] = 1.0 - hcomp[1, j] - hcomp[3, j] - hcomp[4, j]
                lc[j] = es.lcenv[jm]
            else:
                # The last new point, interpolated to the fitting mass.
                if jm == 1:
                    x0 = hs[m]
                    x1 = stotal + senv
                    x2 = es.envs[jm] + stotal
                    if x2 - x0 < 1.0e-14:
                        raise SystemExit(9998)
                    fx = (x1 - x0) / (x2 - x0)
                    hd[j] = hd[m] + fx * (es.envd[jm] - hd[m])
                    hp[j] = hp[m] + fx * (es.envp[jm] - hp[m])
                    hr[j] = hr[m] + fx * (es.envr[jm] - hr[m])
                    hs[j] = x1
                    ht[j] = ht[m] + fx * (es.envt[jm] - ht[m])
                    # Question (see the page): the sign is reversed here.
                    hcomp[1, j] = hcomp[1, m] + fx * (hcomp[1, m]
                                                      - es.envx[jm])
                    hcomp[3, j] = hcomp[3, m] + fx * (hcomp[3, m]
                                                      - es.envz[jm])
                    for k in range(4, jend + 1):
                        hcomp[k, j] = hcomp[k, m]
                    hcomp[2, j] = (1.0 - hcomp[1, j] - hcomp[3, j]
                                   - hcomp[4, j])
                    lc[j] = bool(es.lcenv[jm] or lc[m])
                else:
                    x0 = es.envs[jm - 1] + stotal
                    x1 = stotal + senv
                    x2 = es.envs[jm] + stotal
                    if x2 - x0 < 1.0e-14:
                        raise SystemExit(9998)
                    fx = (x1 - x0) / (x2 - x0)
                    hd[j] = es.envd[jm - 1] + fx * (es.envd[jm]
                                                    - es.envd[jm - 1])
                    hp[j] = es.envp[jm - 1] + fx * (es.envp[jm]
                                                    - es.envp[jm - 1])
                    hr[j] = es.envr[jm - 1] + fx * (es.envr[jm]
                                                    - es.envr[jm - 1])
                    hs[j] = x1
                    ht[j] = es.envt[jm - 1] + fx * (es.envt[jm]
                                                    - es.envt[jm - 1])
                    hcomp[1, j] = es.envx[jm - 1] + fx * (es.envx[jm]
                                                          - es.envx[jm - 1])
                    hcomp[3, j] = es.envz[jm - 1] + fx * (es.envz[jm]
                                                          - es.envz[jm - 1])
                    for k in range(4, jend + 1):
                        hcomp[k, j] = hcomp[k, m]
                    hcomp[2, j] = (1.0 - hcomp[1, j] - hcomp[3, j]
                                   - hcomp[4, j])
                    lc[j] = bool(es.lcenv[jm] or es.lcenv[jm - 1])
                m = j
                reached = True
                break  # GOTO 587
        if not reached:
            m = m + es.numenv
        # label 587: the new shells rotate at the old surface rate.
        if nml.lrot:
            for j in range(mm + 1, m + 1):
                omega[j] = omega[mm]
                hjm[j] = c.cc23 * omega[mm] * 10.0**(2.0 * hr[j])
        print(" NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS\n"
              " J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F")
        for j in range(mm, m + 1):  # format 911
            print(f"{j:5d}" + "".join(
                f"{q:16.8E}" for q in (hd[j], hl[j], hp[j], hr[j],
                                       hs[j] - stotal, ht[j]))
                + f" {'T' if lc[j] else 'F'}")
    lnew = True
    iounits_mod.write(ishort,
                      f"     ***** NEW ENVELOPE MASS CALCULATED *****\n"
                      f"        OLD SENV {senv1:22.13E}  NEW SENV"
                      f"{comp_mod.senv:22.13E}")  # format 597
    return m, lnew, qdt, qdp, vel
