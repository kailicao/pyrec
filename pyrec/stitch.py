"""Python version of stitch.f (YREC 09c3b9b): one model's interior,
envelope and atmosphere written to the .store file (JvS 08/25 format)."""

import math

import numpy as np

from pyrec import atmstruct_mod
from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import envstruct_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import pulse1_mod
from pyrec import roten_mod
from pyrec import scrtch_mod
from pyrec import sound_mod
from pyrec import temp2_mod
from pyrec.envint import envint

JSON = dimensions_mod.json

# Locals kept between calls (SAVE). b is set only in the envelope branch but
# written by the atmosphere's lines, so an atmosphere without an envelope
# writes the value of an earlier call.
b = 0.0
# Question (see the page): LPULPT is never set in STITCH; False stands for
# the Fortran's undefined value.
lpulpt = False
envs1 = np.zeros(JSON + 1)


def _e0(x, w, d):
    """Fortran's 0PEw.d: a mantissa in [0.1, 1), e.g. 0.12345678E+01."""
    if x == 0.0:
        return f"{'0.' + '0' * d + 'E+00':>{w}}"
    e = math.floor(math.log10(abs(x))) + 1
    m = x / 10.0**e
    if round(abs(m), d) >= 1.0:  # rounding carried into the next decade
        m /= 10.0
        e += 1
    s = f"{abs(m):.{d}f}"[1:]  # drop the leading zero, as Fortran does
    sign = "-" if x < 0 else ""
    return f"{sign + '0' + s + f'E{e:+03d}':>{w}}"


def _l1(flag):
    return "T" if flag else "F"


def _fmt62(model, i, hs, hr, hl, hp, ht, hd, om, l1, l2, l3, l4, comp):
    # I6,1X,I6,0P2F18.14,1PE24.16,0PF18.14,1X,F18.14,1X,F18.14,1PE24.16,
    # 1X,L1,1X,L1,1X,L1,1X,L1,3(0PF12.9),12(0PE16.8),2X
    s = (f"{model:6d} {i:6d}{hs:18.14f}{hr:18.14f}{hl:24.16E}{hp:18.14f}"
         f" {ht:18.14f} {hd:18.14f}{om:24.16E}"
         f" {_l1(l1)} {_l1(l2)} {_l1(l3)} {_l1(l4)}")
    s += "".join(f"{comp[j]:12.9f}" for j in range(1, 4))
    s += "".join(_e0(comp[j], 16, 8) for j in range(4, 16))
    return s + "  "


def _fmt63(v):
    # 1PE10.4,1PE11.3,E12.4,E12.4,E12.4,1PE12.4,0PF9.5,1X,F9.5,1X,F9.5,
    # 1X,F9.5,1X,F9.5,1X,F9.5,5E12.4,8E13.5 (the last E's with 0P); the
    # 22 values end the record at the fifth E13.5.
    s = f"{v[0]:10.4E}{v[1]:11.3E}{v[2]:12.4E}{v[3]:12.4E}{v[4]:12.4E}"
    s += f"{v[5]:12.4E}{v[6]:9.5f}"
    s += "".join(f" {v[k]:9.5f}" for k in range(7, 12))
    s += "".join(_e0(v[k], 12, 4) for k in range(12, 17))
    s += "".join(_e0(v[k], 13, 5) for k in range(17, 22))
    return s


def _fmt64(v):
    # 7E14.6, 4E11.3 (0P)
    return ("".join(_e0(x, 14, 6) for x in v[:7])
            + "".join(_e0(x, 11, 3) for x in v[7:]))


def stitch(hcomp, hr, hp, hd, hs, ht, hl, hs1, omega, eta2, hi, r0, hjm, fp,
           ft, teffl, hstot, bl, m, lc, model, cmixl):
    """Write one model's full profile to the .store unit ISTOR, one line of
    61 columns per point: the interior shells, then, if asked (LSTENV,
    LSTATM), the envelope integrated again with ENVINT at fixed output
    steps, and the atmosphere.

    Translates: stitch.f#STITCH (09c3b9b)

    Arguments: the model (hcomp[species, shell], hr, hp, hd, hs, ht, hl,
    hs1, omega, lc, index 0 unused), the rotation arrays eta2, hi, r0, hjm,
    fp, ft, log Teff, log M (hstot), log L/Lsun (bl), m, the model number,
    and cmixl (the mixing length, passed to envint). Nothing is changed.

    Returns None.
    """
    global b
    c = constants_mod
    nml = nmlpar_mod
    sc = scrtch_mod
    es = envstruct_mod
    at = atmstruct_mod
    unit = iounits_mod.istor
    zeros11 = [0.0] * 11

    cg = math.exp(c.cln * c.cgl)
    sg = 0.0
    for i in range(1, m + 1):
        sg = math.exp(c.cln * (c.cgl - 2.0 * hr[i])) * hs1[i]
        line = _fmt62(model, i, hs[i], hr[i], hl[i], hp[i], ht[i], hd[i],
                      omega[i], lc[i], True, False, False, hcomp[:, i])
        # Question (see the page): beta of an interior shell is read from
        # the envelope's EBETAS, not SBETA.
        line += _fmt63([sc.so[i], sg, sc.sdel[1, i], sc.sdel[2, i],
                        sc.sdel[3, i], sc.svel[i], sound_mod.gam1[i],
                        sc.sfxion[1, i], sc.sfxion[2, i], sc.sfxion[3, i],
                        es.ebetas[i], sc.seta[i], sc.seg[1, i], sc.seg[2, i],
                        sc.seg[3, i], sc.seg[4, i], sc.seg[5, i],
                        sc.sesum[i], sc.seg[6, i], sc.seg[7, i], sc.scp[i],
                        pulse1_mod.pqdt[i]])
        if nml.lrot:
            fm = math.exp(c.cln * hs[i])
            duma = c.cc13 * omega[i]**2 / (cg * fm) * 5.0 / (2.0 + eta2[i])
            a = duma * r0[i]**3
            rpoleq = (1.0 - a) / (1.0 + 0.5 * a)
            vtot = temp2_mod.ves[i] + temp2_mod.vgsf[i] + temp2_mod.vss[i]
            line += _fmt64([a, rpoleq, fp[i], ft[i], hjm[i], hi[i],
                            roten_mod.derot[i], temp2_mod.ves[i],
                            temp2_mod.vgsf[i], temp2_mod.vss[i], vtot])
        else:
            line += _fmt64(zeros11)
        iounits_mod.write(unit, line)

    if nml.lstenv:
        # Only provide an envelope if asked to do so: fixed output steps
        # ("dropping a sinkline").
        # Question (see the page): the saved steps are never restored.
        abeg0 = nml.atmbeg
        amin0 = nml.atmmin
        amax0 = nml.atmmax
        ebeg0 = nml.envbeg
        emin0 = nml.envmin
        emax0 = nml.envmax
        nml.atmbeg = nml.atmstp
        nml.atmmin = nml.atmstp
        nml.atmmax = nml.atmstp
        nml.envbeg = nml.envstp
        nml.envmin = nml.envstp
        nml.envmax = nml.envstp
        idum = 0
        b = math.exp(c.cln * bl)
        fpl = fp[m]
        ftl = ft[m]
        katm = 0
        kenv = 0
        ksaha = 0
        ixx = 0
        lprt = True
        lsbc0 = False
        x = hcomp[1, m]
        z = hcomp[3, m]
        rl = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
        gl = c.cgl + hstot - rl - rl
        plim = hp[m]
        # Question (see the page): tests lc[m], where STARIN tests jenv == m.
        if lc[m] and nml.spotf != 0.0 and nml.spotx != 1.0:
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

        for i in range(1, es.numenv + 1):
            envs1[i] = math.exp(c.cln * (es.envs[i] + hstot))
        for i in range(m + 1, m + es.numenv + 1):
            k = i - m
            sg = math.exp(c.cln * (c.cgl - 2.0 * es.envr[k])) * envs1[k]
            line = _fmt62(model, i, es.envs[k] + hstot, es.envr[k],
                          es.envl[k], es.envp[k], es.envt[k], es.envd[k],
                          omega[m], es.lcenv[k], False, True, False,
                          hcomp[:, m])
            line += _fmt63([es.envo[k], sg, es.edels[1, k], es.edels[2, k],
                            es.edels[3, k], es.evels[k], es.egam1[k],
                            es.efxions[1, k], es.efxions[2, k],
                            es.efxions[3, k], es.ebetas[k]]
                           + [0.0] * 9 + [es.eqcp[k], es.eqdt[k]])
            line += _fmt64(zeros11)
            iounits_mod.write(unit, line)

    if nml.lstatm:
        # Question (see the page): without LSTENV, numenv, envr and b are
        # those of an earlier integration or call; sg is the last point's.
        for i in range(at.numatm, 0, -1):
            rad = math.log10(math.exp(c.cln * es.envr[es.numenv])
                             + at.atmor[i])
            line = _fmt62(model, at.numatm - i + m + es.numenv, hstot, rad,
                          b, at.atmop[i], at.atmot[i], at.atmod[i], omega[m],
                          False, False, False, True, hcomp[:, m])
            line += _fmt63([at.atmoo[i], sg, at.adels[1, i], at.adels[2, i],
                            at.adels[3, i], 0.0, at.agam1[i],
                            at.afxions[1, i], at.afxions[2, i],
                            at.afxions[3, i], at.abetas[i]]
                           + [0.0] * 9 + [at.atmocp[i], at.aqdt[i]])
            line += _fmt64(zeros11)
            iounits_mod.write(unit, line)
