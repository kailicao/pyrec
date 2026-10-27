"""Python version of wrtmod.f (YREC 09c3b9b): extra profiles of one model
(the sound-speed table, an envelope integrated again for printing, and the
model for pulsation codes)."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import label_mod
from pyrec import nmlpar_mod
from pyrec import pualpha_mod
from pyrec import pulse1_mod
from pyrec import scrtch_mod
from pyrec import sound_mod
from pyrec.envint import envint

JSON = dimensions_mod.json

# Locals kept between calls (SAVE).
# Question (see the page): IDM and ID are never assigned, so the pulsation
# model's shell list is undefined in the Fortran; here it starts empty.
id_ = np.zeros(JSON + 1, dtype=int)
_saved = {'idm': 0,
          # Question (see the page): LPRT is set only if LSTORE; otherwise
          # ENVINT gets its saved (or, on the first call, undefined) value.
          'lprt': False}


def _l1(flag):
    return "T" if flag else "F"


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


def _records(values, per_line, w, d):
    """Values in 0PEw.d, per_line[k] of them on record k (FORMATs 5052,
    6052, 6053); the last record stops when the values run out."""
    lines = []
    k = 0
    for count in per_line:
        chunk = values[k:k + count]
        if not chunk:
            break
        lines.append("".join(_e0(v, w, d) for v in chunk))
        k += count
    return "\n".join(lines)


def wrtmod(m, jenv, hcomp, hs1, hd, hl, hp, hr, ht, model, bl, teffl, fp,
           ft, hs, dage, cmixl, xmsol):
    """Write extra profiles of one model: with LSOUND, a sound-speed and
    structure table (units IMODPT and ISOUND); with LSTATM, LSTENV or
    LPULSE, an envelope integrated again by ENVINT for printing; with
    LPULSE, a model for pulsation codes on IOPMOD (version IPVER).

    Translates: wrtmod.f#WRTMOD (09c3b9b)

    Arguments: the number of shells m, the envelope's base jenv, the model
    (hcomp, hs1, hd, hl, hp, hr, ht, fp, ft, hs, all logarithmic as in the
    model arrays), the model number, log L bl, log Teff teffl, the age
    dage, the mixing length cmixl, and xmsol, the star's mass in solar
    masses as main.f snapshots it. Nothing is returned.
    """
    n = nmlpar_mod
    c = constants_mod
    cln = c.cln
    io = iounits_mod
    gam1 = sound_mod.gam1
    pqdp = pulse1_mod.pqdp
    pqdt = pulse1_mod.pqdt
    sdel = scrtch_mod.sdel

    # Question (see the page): debugging prints at every call.
    print(' wrtmod LSOUND 1: ', _l1(sound_mod.lsound))

    if sound_mod.lsound:
        print(' wrtmod LSOUND 2: ', _l1(sound_mod.lsound))
        # CFD 10/09: the sound-speed profile, on ISOUND (opened in PARMIN).
        # Question (see the page): the header is written at every call.
        io.write(io.isound, " Shell#     R/Rsun       M/Msun       Cs"
                            "               Rho")  # FORMAT 9124
        for i in range(1, m + 1):
            fr = math.exp(cln * hr[i]) / c.crsun
            fm = math.exp(cln * hs[i]) / c.cmsun
            xx1 = fm / fr**3
            xxx = -math.exp(cln * (c.cgl + hs[i] + hd[i] - hp[i] - hr[i]))
            xx2 = -xxx / gam1[i]
            xx3 = gam1[i]
            xx4 = -xx2 - xxx * (pqdp[i] + sdel[2, i] * pqdt[i])
            xx5 = math.exp(cln * (c.c4pil + hd[i] + 3.0 * hr[i] - hs[i]))
            sv = 1.0e-5 * math.sqrt(gam1[i] * math.exp(cln * (hp[i] - hd[i])))
            io.write(io.imodpt,  # FORMAT 123
                     f" {fr:12.8f}{fm:12.8f}{xx1:16.8E}{xx2:16.8E}"
                     f"{xx3:16.8E}{xx4:16.8E}{xx5:16.8E}{sv:16.8E}")
            io.write(io.isound,  # FORMAT 9123
                     f" {i:6d} {fr:12.8f}{fm:12.8f}{sv:16.8E}{hd[i]:16.8E}"
                     f"{ht[i]:16.8E}")

            # Derivatives of dP/dr, drho/dr.
            if i < m:
                rmid = 0.5 * (math.exp(cln * hr[i]) + math.exp(cln * hr[i + 1]))
                dr = math.exp(cln * hr[i + 1]) - math.exp(cln * hr[i])
                divp = 0.5 * (gam1[i] * math.exp(cln * hp[i])
                              + gam1[i + 1] * math.exp(cln * hp[i + 1])) * dr
                divr = 0.5 * (math.exp(cln * hd[i])
                              + math.exp(cln * hd[i + 1])) * dr
                qpr1 = math.exp(cln * (c.cgl + hs[i] + hd[i] - 2.0 * hr[i]))
                qdr1 = math.exp(cln * (hd[i] - hp[i])) * pqdp[i] * qpr1
                qpr2 = math.exp(cln * (c.cgl + hs[i + 1] + hd[i + 1]
                                       - 2.0 * hr[i + 1]))
                # Question (see the page): QPR1 where QPR2 is expected.
                qdr2 = math.exp(cln * (hd[i + 1] - hp[i + 1])) * pqdp[i + 1] \
                    * qpr1
                qqpr = (qpr1 - qpr2) / divp
                qqdr = (qdr1 - qdr2) / divr
                io.write(io.imodpt,  # FORMAT 124
                         f" {rmid / c.crsun:11.7f}{qpr1:16.8E}{qdr1:16.8E}"
                         f"{qqpr:16.8E}{qqdr:16.8E}")

    # Integrate an envelope from the surface to the converged model,
    # printing the results (ENVINT prints them).
    if n.lstatm or n.lstenv or n.lpulse:
        lsbc0 = False
        if n.lstore:
            _saved['lprt'] = True
        lprt = _saved['lprt']
        katm = 0
        kenv = 0
        ksaha = 0
        # Save the integration step parameters and enforce the spacing
        # requested for printout.
        abeg0 = n.atmbeg
        amin0 = n.atmmin
        amax0 = n.atmmax
        ebeg0 = n.envbeg
        emin0 = n.envmin
        emax0 = n.envmax
        n.atmbeg = n.atmstp
        n.atmmin = n.atmstp
        n.atmmax = n.atmstp
        n.envbeg = n.envstp
        n.envmin = n.envstp
        n.envmax = n.envstp
        b = math.exp(cln * bl)
        rl = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
        gl = c.cgl + comp_mod.stotal - rl - rl
        x = hcomp[1, m]
        z = hcomp[3, m]
        fpl = fp[m]
        ftl = ft[m]
        ixx = 0
        hstot = comp_mod.stotal
        plim = hp[m]
        lpulpt = n.lpulse
        if n.ldh:
            debhu_mod.xxdh = hcomp[1, m]
            debhu_mod.yydh = hcomp[2, m] + hcomp[4, m]
            debhu_mod.zzdh = hcomp[3, m]
            debhu_mod.zdh[1] = hcomp[5, m] + hcomp[6, m]
            debhu_mod.zdh[2] = hcomp[7, m] + hcomp[8, m]
            debhu_mod.zdh[3] = hcomp[9, m] + hcomp[10, m] + hcomp[11, m]
        idum = 0
        # G. Somers 10/14: for a spotted star, the ambient temperature.
        if jenv == m and n.spotf != 0.0 and n.spotx != 1.0:
            ateffl = teffl - 0.25 * math.log10(
                n.spotf * n.spotx**4.0 + 1.0 - n.spotf)
        else:
            ateffl = teffl
        dum1 = np.zeros(5)
        dum2 = np.zeros(4)
        dum3 = np.zeros(4)
        dum4 = np.zeros(4)
        rl, idum, katm, kenv, ksaha = envint(
            b, fpl, ftl, gl, hstot, ixx, lprt, lsbc0, plim, rl, ateffl, x,
            z, dum1, idum, katm, kenv, ksaha, dum2, dum3, dum4, lpulpt,
            cmixl)
        n.atmbeg = abeg0
        n.atmmin = amin0
        n.atmmax = amax0
        n.envbeg = ebeg0
        n.envmin = emin0
        n.envmax = emax0

    fsi = math.exp(-cln * comp_mod.stotal)
    idm = _saved['idm']
    # The pulsation model's header.
    if n.lpulse:
        rsurfl = 0.5 * (bl - c.c4pil - c.csigl - 4.0 * teffl + c.clsunl)
        tempr = rsurfl - c.crsunl
        qsmass = xmsol
        io.write(io.iopmod,  # FORMAT 5001
                 f" MODEL#={model:5d}  NUMBER OF SHELLS IN MODEL={idm:5d}"
                 f" VER={n.ipver:2d}\n"
                 f" MASS={qsmass:8.5f}  LOG(TEFF)={teffl:8.5f}\n"
                 f" LOG(L/LSUN)={bl:16.10f}  LOG(R/RSUN)={tempr:16.10f}\n"
                 f" AGE={dage:12.5E} GYR\n"
                 f" MIXING LENGTH PARAMETER={cmixl:16.10f}\n"
                 f" ZAMS (X,Z)={label_mod.xenv0:18.10f}"
                 f"{label_mod.zenv0:18.10f}")

    # The pulsation model, shell by shell, on the list ID(1..IDM).
    p = pulse1_mod
    s = scrtch_mod
    v = pualpha_mod
    for j in range(1, idm + 1):  # loop 220
        i = id_[j]
        fs = fsi * hs1[i]
        if not (n.lpulse and n.lstore):
            continue
        if j == 2 and i == 1:  # GOTO 5003
            continue
        if n.ipver == 1:
            pelpf = c.cgas * math.exp(cln * (ht[i] + hd[i])) * p.pemu[i]
            vals = [hr[i], fs, hl[i], ht[i], hd[i],
                    hp[i], s.sesum[i], s.so[i], p.pqdp[i], p.pqed[i],
                    p.pqet[i], p.pqod[i], p.pqot[i], s.sdel[2, i],
                    s.sdel[3, i],
                    p.pqcp[i], p.prmu[i], p.pqdt[i], pelpf]
            io.write(io.iopmod, _records(vals, (5, 5, 5, 5), 16, 9))
        elif n.ipver == 2:
            vals = [hr[i], fs, hl[i], ht[i], hd[i],
                    hp[i], s.sesum[i], s.so[i], p.pqdp[i], p.pqed[i],
                    p.pqet[i], p.pqod[i], p.pqot[i], s.sdel[2, i],
                    s.sdel[3, i],
                    p.pqcp[i], p.prmu[i], p.pqdt[i], hcomp[1, i],
                    hcomp[3, i]]
            io.write(io.iopmod, _records(vals, (5, 5, 5, 5), 23, 16))
        elif n.ipver == 3:
            # DBG 7/95: with the mixing-length variables.
            vals = [hr[i], fs, hl[i], ht[i], hd[i], v.valfmlt[i],
                    hp[i], s.sesum[i], s.so[i], p.pqdp[i], p.pqed[i],
                    v.vphmlt[i],
                    p.pqet[i], p.pqod[i], p.pqot[i], s.sdel[2, i],
                    s.sdel[3, i], v.vcmxmlt[i],
                    p.pqcp[i], p.prmu[i], p.pqdt[i], hcomp[1, i],
                    hcomp[3, i]]
            io.write(io.iopmod, _records(vals, (6, 6, 6, 5), 23, 16))
