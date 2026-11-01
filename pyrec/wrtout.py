"""Python version of wrtout.f (YREC 09c3b9b): the output of each converged
model (summary block, track lines, last model, detailed profiles), with
the summary quantities computed on the way."""

import math

import numpy as np

from pyrec import acdpth_mod
from pyrec import cent_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import fluxes_mod
from pyrec import iounits_mod
from pyrec import label_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec import scrtch_mod
from pyrec.calcad import calcad
from pyrec.eqstat import eqstat
from pyrec.gettau import gettau
from pyrec.putstore import putstore
from pyrec.wrtlst import wrtlst
from pyrec.wrtmil import wrtmil
from pyrec.wrtmod import wrtmod

# DATA: the track's cross sections for the Cl and Ga rates (SNU), for the
# fluxes pp, pep, hep, Be7, B8, N13, O15, F17.
GASNUF = np.array([0.0, 1.18e1, 2.15e2, 7.14e4, 7.17e1, 2.40e4, 6.04e1,
                   1.137e2, 1.139e2])
CLSNUF = np.array([0.0, 0.0, 1.6e1, 4.26e4, 2.4, 1.14e4, 1.7, 6.8, 6.9])

# SAVE: locals that a call can use before setting them, so their values
# come from an earlier model (see the page's questions 1, 3 and 5).
icheck = 0  # DATA ICHECK/0/
fx = 0.0
czsurf_rl = 0.0
rotp = 0.0
rotv = 0.0
hicz = 0.0
sss1 = sss2 = sss3 = 0.0
rrr1 = rrr2 = rrr3 = 0.0
ksaha = 0  # Question (see the page): never set.


# --- Fortran edit descriptors -------------------------------------------

def _e(v, w, d, e=2, scale=1):
    """Fortran Ew.d (Ee) under the scale factor 1P (scale=1: d.ddd) or 0P
    (scale=0: 0.ddd). With e = 2 an exponent beyond 99 drops the 'E', as
    Fortran does; a field too narrow is filled with '*'."""
    digits = d + 1 if scale == 1 else d
    m, x = f"{v:.{digits - 1}E}".split('E')
    x = int(x)
    if scale == 0:
        sign = '-' if m.startswith('-') else ''
        m = sign + '0.' + m.lstrip('-').replace('.', '')
        x = 0 if v == 0 else x + 1
    if e == 2 and abs(x) > 99:
        text = f"{m}{x:+04d}"
    else:
        text = f"{m}E{x:+0{e + 1}d}"
    return text.rjust(w) if len(text) <= w else '*' * w


def _es(values, w, d, e=2, scale=1):
    """A repeated E descriptor."""
    return ''.join(_e(v, w, d, e, scale) for v in values)


def _f(v, w, d):
    """Fortran Fw.d (Fw.0 keeps its decimal point)."""
    text = f"{v:.{d}f}" + ('.' if d == 0 else '')
    return text.rjust(w) if len(text) <= w else '*' * w


def _l(flag):
    return 'T' if flag else 'F'


# --- The track lines ----------------------------------------------------

def _model_block(model, m, dage, bl, rl, gl, teffl, corem, czsurf_m,
                 w, d, e):
    """The common start of FORMATs 1499, 1500, 1800: 1X,2I8, then seven E
    fields (DAGE ... CZsurf_M)."""
    return (f" {model:8d}{m:8d}"
            + _es([dage, bl, rl, gl, teffl, corem, czsurf_m], w, d, e))


def _hcomp_columns(hcomp, m):
    """HCOMP(4..11,1), HCOMP(4..15,M), HCOMP(1..3,M), Z/X at the surface."""
    return ([hcomp[i, 1] for i in range(4, 12)]
            + [hcomp[i, m] for i in range(4, 16)]
            + [hcomp[i, m] for i in range(1, 4)]
            + [hcomp[3, m] / hcomp[1, m]])


def wrtout(hcomp, hd, hl, hp, hr, hs, hs1, ht, lc, m, model, dage,
           ddage, smass, teffl, bl, gl, lshell, jxbeg, jxmid, jxend, jcore,
           jenv, tlumx, trit, tril, ps, ts, rs, cfenv, ftri, hstot, omega,
           lpunch, fp, ft, eta2, r0, hjm, hi, sjtot, skerot, hs2, cmixl,
           xmsol, compmix):
    """Write the output of a converged model: the summary block to ISHORT,
    the track line(s) to ITRACK (by ITRVER), the last model to ILAST, and at
    intervals the detailed profiles (PUTSTORE, or WRTMIL and WRTMOD) and the
    acoustic depths (CALCAD). On the way it rescales the luminosity terms,
    finds the central conditions and the base of the surface convection
    zone (stored in CENT_MOD), and the neutrino capture rates.

    Translates: wrtout.f#WRTOUT (09c3b9b)

    Arguments: the model (hcomp, hd, hl, hp, hr, hs, hs1, ht, lc, m, the
    rotation arrays omega, fp, ft, eta2, r0, hjm, hi, hs2), its run data
    (model, dage, ddage, smass, teffl, bl, hstot, sjtot, skerot), the shell
    and zone indices (lshell, jxbeg, jxmid, jxend, jcore, jenv), the
    luminosity terms tlumx (rescaled in place), the envelope data passed to
    WRTLST/PUTSTORE (trit, tril, ps, ts, rs, cfenv, ftri), and cmixl, xmsol,
    compmix, passed on. Returns (gl, lpunch): log g, set here, and lpunch,
    cleared when the model is stored.
    """
    global icheck, fx, czsurf_rl, rotp, rotv, hicz
    global sss1, sss2, sss3, rrr1, rrr2, rrr3
    n = nmlpar_mod
    c = constants_mod
    sc = scrtch_mod
    fl = fluxes_mod
    ac = acdpth_mod
    ce = cent_mod
    ishort = iounits_mod.ishort
    itrack = iounits_mod.itrack
    w = iounits_mod.write

    # Question (see the page): LTSCAL is set and never used.
    ltscal = False

    # Rescale the luminosity terms to the surface luminosity.
    if not n.lkuthe:
        thl = 0.0
        for i in range(1, 9):
            thl = thl + tlumx[i]
        temp = hl[m] / thl
        for i in range(1, 9):  # label 10
            tlumx[i] = tlumx[i] * temp

    # The summary block on ISHORT: global properties and the center.
    w(ishort, "\n " + '*' * 127)  # 21
    if not n.lkuthe:
        step = _f(ddage, 12, 0)  # 30
    else:
        step = _e(ddage, 12, 4)  # 40
    w(ishort, f" MODEL NO.{model:5d}  MASS{_f(smass, 13, 7)}  (X,Z)=("
              f"{_f(comp_mod.xnew, 11, 9)},{_f(comp_mod.znew, 11, 9)})  "
              f"AGE(GYRS){_f(dage, 14, 8)} STEP(YRS)={step}")

    # The convective core's mass.
    if jcore > 1:
        corem = hs1[jcore] / c.cmsun
    else:
        corem = 0.0
    # Its radius (JVS 10/11).
    if jcore > 1:
        # Question (see the page): FX2 is 1 (or 0/0), and FX is used.
        fx2 = ((sc.sdel[3, jcore + 1] - sc.sdel[1, jcore])
               / (sc.sdel[3, jcore + 1] - sc.sdel[1, jcore]))
        corerl = hr[jcore] + fx * (hr[jcore + 1] - hr[jcore]) - c.crsunl
        corer = math.exp(c.cln * corerl)
    else:
        corer = 0.0

    # Bolometric magnitude, radius and gravity.
    bolm = c.cmbol - 2.5 * bl
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    gl = c.cgl + comp_mod.stotal - rl - rl

    # Central conditions, extrapolated from shell 1.
    temp = 0.5 * math.exp(c.cln * (c.cc13 * (c.c4pi3l + hd[1] - hs[1])
                                   + hd[1] + c.cgl + hs[1]))
    p = math.exp(c.cln * hp[1])
    pl = math.log10(p + temp)
    tl = ht[1] + math.log10(1.0 + temp * sc.sdel[2, 1] / p)
    dl = hd[1]
    x = hcomp[1, 1]
    z = hcomp[3, 1]
    latmo = True
    lderiv = False
    if n.ldh:
        dh = debhu_mod
        dh.xxdh = hcomp[1, 1]
        dh.yydh = hcomp[2, 1] + hcomp[4, 1]
        dh.zzdh = hcomp[3, 1]
        dh.zdh[1] = hcomp[5, 1] + hcomp[6, 1]
        dh.zdh[2] = hcomp[7, 1] + hcomp[8, 1]
        dh.zdh[3] = hcomp[9, 1] + hcomp[10, 1] + hcomp[11, 1]
    fxion = np.zeros(4)
    # Question (see the page): KSAHA is never set here.
    (t, p, dl, d, beta, betai, beta14, rmu, amu, emu, eta, qdt, qdp, qcp,
     dela, qdtt, qdtp, qat, qap, qcpt, qcpp, _) = eqstat(
        tl, 0.0, pl, p, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        lderiv, latmo, ksaha)
    ce.pcenter = pl
    ce.tcenter = tl
    ce.dcenter = dl

    # The base of the surface convection zone, where DEL_rad - DEL_ad
    # changes sign between JENV-1 and JENV.
    if jenv < m:
        if jenv > 1:
            dd2 = sc.sdel[1, jenv - 1] - sc.sdel[3, jenv - 1]
            dd1 = sc.sdel[1, jenv] - sc.sdel[3, jenv]
            fx = dd2 / (dd2 - dd1)
            hsb = hs1[jenv - 1] + fx * (hs1[jenv] - hs1[jenv - 1])
            ce.czsurf_m = (math.exp(c.cln * hstot) - hsb) / c.cmsun
            czsurf_rl = hr[jenv - 1] + fx * (hr[jenv] - hr[jenv - 1]) - c.crsunl
            ce.czsurf_r = math.exp(c.cln * czsurf_rl)
            czsurf_o = sc.so[jenv - 1] + fx * (sc.so[jenv] - sc.so[jenv - 1])
            czsurf_tl = ht[jenv - 1] + fx * (ht[jenv] - ht[jenv - 1])
            czsurf_dl = hd[jenv - 1] + fx * (hd[jenv] - hd[jenv - 1])
            czsurf_pl = hp[jenv - 1] + fx * (hp[jenv] - hp[jenv - 1])
            czsurf_t = math.exp(c.cln * czsurf_tl)
            czsurf_d = math.exp(c.cln * czsurf_dl)
            czsurf_p = math.exp(c.cln * czsurf_pl)
        else:
            # Fully convective.
            ce.czsurf_m = smass
            ce.czsurf_r = 0.0
            czsurf_t = 10.0 ** ce.tcenter
            czsurf_d = 10.0 ** ce.dcenter
            czsurf_p = 10.0 ** ce.pcenter
            czsurf_o = sc.so[1]
    else:
        # No surface convection zone.
        ce.czsurf_m = 0.0
        ce.czsurf_r = 0.0
        czsurf_t = 0.0
        czsurf_d = 0.0
        czsurf_p = 0.0
        czsurf_o = 0.0
    w(ishort, f" SHELLS={m:5d}  (X0,Z0)=({_f(label_mod.xenv0, 9, 7)},"
              f"{_f(label_mod.zenv0, 9, 7)})  CONV. ZONE MASSES(MSUN): CORE"
              f"{_f(corem, 10, 7)} ENV.{_f(ce.czsurf_m, 10, 7)}"
              f" RAD. FRAC.{_f(ce.czsurf_r, 10, 7)}")  # 50
    rl = rl - c.crsunl
    w(ishort, f" LOG(TEFF)={_f(teffl, 11, 8)}  M(BOL)={_f(bolm, 11, 7)}"
              f"  LOG(L/LSUN)={_f(bl, 12, 8)}  LOG(R/RSUN)={_f(rl, 12, 8)}"
              f"  LOG(G) ={_f(gl, 12, 8)}")  # 60
    w(ishort, f" CENTER: LOG P={_f(pl, 10, 7)} LOG T={_f(tl, 10, 8)} LOG D="
              f"{_f(dl, 10, 6)} BETA={_f(beta, 9, 7)} ETA={_f(eta, 10, 5)}"
              f"  X={_f(x, 9, 7)} Z={_f(z, 9, 7)} O16={_f(hcomp[9, 1], 9, 7)}")
    # FORMAT 80: PPI, PPII, PPIII, CNO; triple alpha, He-C, neutrinos, grav.
    w(ishort, f" ENERGY: PPI{_e(tlumx[1], 13, 6)}  PPII{_e(tlumx[2], 13, 6)}"
              f"  PPIII{_e(tlumx[3], 13, 6)}  CNO{_e(tlumx[4], 13, 6)}\n"
              f"{'':9s}TRIPLE ALPHA{_e(tlumx[5], 13, 6)}"
              f"  HE-C{_e(tlumx[8], 13, 6)}  NEUTRINOS{_e(tlumx[6], 13, 6)}"
              f"  GRAV{_e(tlumx[7], 13, 6)}")

    # The hydrogen shell and the temperature maximum.
    ss3 = 0.0
    ss4 = 0.0
    if lshell:
        ss1 = hs1[jxmid] / c.cmsun
        ss2 = (hs1[jxend] - hs1[jxbeg - 1]) / c.cmsun
        ss3 = hs1[jxbeg - 1] / c.cmsun
        tlmax = tl
        line120 = None
        if not n.lkuthe:
            # The first temperature maximum from the center.
            for i in range(2, m + 1):  # label 100
                if ht[i] < ht[i - 1]:
                    break
            else:
                i = m + 1
            itmax = i - 1  # label 110
            if itmax > 1:
                ss4 = hs1[itmax] / c.cmsun
                tlmax = ht[itmax]
            else:
                ss4 = 0.0
                tlmax = ht[1]
            line120 = True
        else:
            # Helium flash: the first maximum above log T = 7.98.
            # Question (see the page): the threshold 7.98.
            for i in range(2, m + 1):  # label 130
                if ht[i] < ht[i - 1] and ht[i - 1] > 7.98:
                    break
            else:
                i = m + 1
            itmax = i - 1  # label 140
            if itmax > 1:
                ss4 = hs1[itmax] / c.cmsun
                tlmax = ht[itmax]
                rlmax = hr[itmax]
                lcmax = lc[itmax]
                line120 = True
        if line120:
            w(ishort, f" H-SHELL MID-PT={_f(ss1, 10, 7)} MASS TOTAL="
                      f"{_f(ss2, 10, 7)}  HE-CORE MASS={_f(ss3, 10, 7)} MAX-T="
                      f"{_f(tlmax, 10, 7)} (MASS={_f(ss4, 9, 7)})")  # 120
            if n.lkuthe:
                w(ishort, f" CONVECTION = {_l(lcmax)}     LOG(R) MAX-T ="
                          f"{_f(rlmax, 8, 5)}")  # 150

    # Neutrino fluxes and detector capture rates.
    fluxtot = fl.fluxtot
    w(ishort, " NEUTRINOS 1E10ERG/CM^2 PP,PEP,HEP,BE7,B8,N13,O15,F17:"
              + _es(fluxtot[1:9], 9, 2))  # 160
    # Question (see the page): these cross sections differ from CLSNUF and
    # GASNUF, and the line calls the rates "energies".
    f = fluxtot
    fl7li = (0.0 * f[1] + 665.0 * f[2] + 8.4e4 * f[3] + 9.6 * f[4]
             + 3.9e4 * f[5] + 42.4 * f[6] + 246.0 * f[7] + 249.0 * f[8])
    fl37cl = (0.0 * f[1] + 16.0 * f[2] + 4.26e4 * f[3] + 2.4 * f[4]
              + 1.09e4 * f[5] + 1.7 * f[6] + 6.8 * f[7] + 6.9 * f[8])
    fl71ga = (11.8 * f[1] + 215.0 * f[2] + 7.3e4 * f[3] + 73.2 * f[4]
              + 2.43e4 * f[5] + 61.8 * f[6] + 116.0 * f[7] + 117.0 * f[8])
    fl81br = (0.0 * f[1] + 75.0 * f[2] + 9.0e4 * f[3] + 18.3 * f[4]
              + 2.7e4 * f[5] + 14.5 * f[6] + 36.7 * f[7] + 37.0 * f[8])
    fl98mo = (0.0 * f[1] + 0.0 * f[2] + 10.0e4 * f[3] + 0.0 * f[4]
              + 3.0e4 * f[5] + 0.0 * f[6] + 0.0 * f[7] + 0.0 * f[8])
    fl115in = (78.0 * f[1] + 576.0 * f[2] + 6.1e4 * f[3] + 248.0 * f[4]
               + 2.5e4 * f[5] + 224.0 * f[6] + 355.0 * f[7] + 356.0 * f[8])
    w(ishort, f" NEUTRINO ENERGIES (1.E-36ERG): 7Li={_e(fl7li, 9, 2)}"
              f" 37Cl={_e(fl37cl, 9, 2)} 71Ga={_e(fl71ga, 9, 2)}"
              f" 81Br={_e(fl81br, 9, 2)} 98Mo={_e(fl98mo, 9, 2)}"
              f" 115In={_e(fl115in, 9, 2)}")  # 2160
    fitmas = hs1[m] / c.cmsun
    w(ishort, f" FIT-POINT    M/MSUN={_f(fitmas, 16, 12)}     (P,T,R) ="
              f"{_f(hp[m], 12, 7)}{_f(ht[m], 12, 7)}{_f(hr[m], 12, 7)}")
    w(ishort, " " + '*' * 127)  # 20

    if n.ltrack:
        # The total moment of inertia.
        hitot = 0.0
        if not n.lrot:
            for i in range(1, m + 1):
                hitot = hitot + c.cc23 * hs2[i] * math.exp(2.0 * c.cln * hr[i])
        else:
            for i in range(1, m + 1):
                hitot = hitot + hi[i]

        if n.itrver == 0:
            # The Cl and Ga rates in SNU.
            fl.clsnu = 0.0
            fl.gasnu = 0.0
            if n.lsnu:
                for i in range(1, 9):
                    fl.clsnu = fl.clsnu + CLSNUF[i] * fluxtot[i]
                    fl.gasnu = fl.gasnu + GASNUF[i] * fluxtot[i]
            else:
                # Question (see the page): FLUXTOT zeroed for later readers.
                for i in range(1, 11):
                    fluxtot[i] = 0.0
            # Rotation period (days), velocity (km/s), the zone's inertia.
            hicz = 0.0
            if n.lrot:
                # Question (see the page): the old cap of 9999 days.
                rotp = min(9999.0, 0.5 * c.c4pi / omega[m] / 8.64e4)
                rotv = omega[m] * math.exp(c.cln * (rl + c.crsunl)) * 1.0e-5
                if jenv < m:
                    for k in range(jenv, m + 1):
                        hicz = hicz + hi[k]
            else:
                rotp = 0.0
                rotv = 0.0
                if jenv < m:
                    for k in range(jenv, m + 1):
                        hicz = hicz + c.cc23 * hs2[k] * math.exp(
                            2.0 * c.cln * hr[k])
            _shell_location(lshell, jxbeg, jxmid, jxend, hs1, hr, rl)
            # tau_cz and P_phot; GETTAU also sets CZsurf_R.
            ce.czsurf_r = gettau(hcomp, hr, hp, hd, hs1, ht, fp, ft, teffl,
                                 hstot, bl, m, lc, ce.czsurf_r, cmixl)
            ovrtrn_mod.taucz0 = ovrtrn_mod.taucz
            ovrtrn_mod.pphot0 = ovrtrn_mod.pphot
            # FORMAT 1499 (E17.8E3, for three-digit exponents) (MCR 12/25).
            w(itrack, _model_block(model, m, dage, bl, rl, gl, teffl, corem,
                                   ce.czsurf_m, 17, 8, 3)
                      + _es([ce.czsurf_r, czsurf_t, czsurf_d, czsurf_p,
                             czsurf_o], 12, 4)
                      + _es([tl, dl, pl, beta, eta, hcomp[1, 1], hcomp[2, 1],
                             hcomp[3, 1]] + list(tlumx[1:6])
                            + [tlumx[8], tlumx[7], tlumx[6]], 17, 8, 3)
                      + _es([fl.clsnu, fl.gasnu] + list(fluxtot[1:11]), 10, 3)
                      + _es(_hcomp_columns(hcomp, m)
                            + [sjtot, skerot, hitot, hicz, omega[m], omega[1],
                               rotp, rotv, ovrtrn_mod.taucz,
                               sss1, sss2, sss3, rrr1, rrr2, rrr3,
                               ovrtrn_mod.pphot, smass], 17, 8, 3))
        elif n.itrver == 1 or n.itrver == 2:
            w(itrack, _model_block(model, m, dage, bl, rl, gl, teffl, corem,
                                   ce.czsurf_m, 16, 8, 2)
                      + _f(ce.czsurf_r, 8, 4)
                      + _es([czsurf_t, czsurf_d, czsurf_p, czsurf_o],
                            12, 4))  # 1500
            w(itrack, " " + _es([tl, dl, pl, beta, eta, hcomp[1, 1],
                                 hcomp[2, 1], hcomp[3, 1], hitot], 16, 8))
            w(itrack, " " + _es(list(tlumx[1:6])
                                + [tlumx[8], tlumx[7], tlumx[6]], 16, 8))
            fl.clsnu = 0.0
            fl.gasnu = 0.0
            for i in range(1, 9):
                fl.clsnu = fl.clsnu + CLSNUF[i] * fluxtot[i]
                fl.gasnu = fl.gasnu + GASNUF[i] * fluxtot[i]
            w(itrack, f" {_f(fl.clsnu, 8, 3)}{_f(fl.gasnu, 8, 3)}"
                      + _es(fluxtot[1:11], 10, 3))  # 1515
            w(itrack, " " + _es([hcomp[i, 1] for i in range(4, 12)], 16, 8))
            surf = _hcomp_columns(hcomp, m)[8:]
            w(itrack, " " + _es(surf[:8], 16, 8) + "\n "
                      + _es(surf[8:], 16, 8))  # 1520
            # Question (see the page): without rotation ROTP, ROTV and HICZ
            # keep earlier values, and version 2 writes them.
            if n.lrot:
                rotp = 0.5 * c.c4pi / omega[m] / 8.64e4
                rotv = omega[m] * math.exp(c.cln * (rl + c.crsunl)) * 1.0e-5
                hicz = 0.0
                if lc[m]:
                    for k in range(jenv, m + 1):
                        hicz = hicz + hi[k]
                else:
                    hicz = 0.0
            _shell_location(lshell, jxbeg, jxmid, jxend, hs1, hr, rl)
            w(itrack, " " + _es([sss1, sss2, sss3, rrr1, rrr2, rrr3],
                                16, 8))  # 1530
        if n.itrver == 2:
            # FORMAT 1540, (1X,1P6E13.5,0P,2F11.5,1E13.5): the last field
            # is under 0P, so it prints as 0.ddddd.
            w(itrack, " " + _es([sjtot, skerot, hitot, hicz, omega[m],
                                 omega[1]], 13, 5)
                      + _f(rotp, 11, 5) + _f(rotv, 11, 5)
                      + _e(ovrtrn_mod.taucz, 13, 5, scale=0))
        if n.itrver == 3:
            w(itrack, f" {model:8d}{m:8d}"
                      + _es([dage, bl, rl, gl, teffl], 13, 5)
                      + _es([corem, ce.czsurf_m], 11, 3)
                      + _f(ce.czsurf_r, 8, 4)
                      + _e(comp_mod.xnew, 13, 5))  # 1501

    # The last model, every time.
    iwrite = iounits_mod.ilast
    wrtlst(iwrite, hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps,
           ts, rs, cfenv, ftri, tlumx, jcore, jenv, model, m, smass, teffl,
           bl, hstot, dage, ddage, omega, cmixl, compmix)

    # The detailed output every NPRTMOD models (PUTSTORE also writes the
    # stitched profile); otherwise the pulsation files.
    stored = n.lstore and model % n.nprtmod == 0
    if stored:
        putstore(hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps,
                 ts, rs, cfenv, ftri, tlumx, jcore, jenv, model, m, smass,
                 teffl, bl, hstot, dage, ddage, omega, hs1, eta2, r0, fp, ft,
                 hjm, hi, cmixl, xmsol, compmix)
        lpunch = False
    if not stored and n.lpulse:
        if n.lmilne:
            wrtmil(hcomp, hd, hl, hp, hr, hs1, m, model)
        wrtmod(m, jenv, hcomp, hs1, hd, hl, hp, hr, ht, model, bl, teffl,
               fp, ft, hs, dage, cmixl, xmsol)

    # Acoustic depths, and a second track line (JVS).
    if ac.lacout:
        if jenv > 1 and ac.ladon:
            # Question (see the page): CZsurf_RL may be from an earlier model.
            calcad(hr, czsurf_rl, m, hd, hp, ht, bl, fp, ft, hstot, teffl,
                   hcomp, dage, jenv, cmixl)
        elif jenv == 1:
            ac.tauczn = 0.0
        if lc[m]:
            icheck = 1
        elif not lc[m]:
            icheck = 0
        if ac.ljlast:
            iwrite = iounits_mod.ijlast
            wrtlst(iwrite, hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps,
                   ts, rs, cfenv, ftri, tlumx, jcore, jenv, model, m, smass,
                   teffl, bl, hstot, dage, ddage, omega, cmixl, compmix)
        # FORMAT 1800. Question (see the page): a second line format in the
        # track file, with values that may be from earlier models.
        w(itrack, _model_block(model, m, dage, bl, rl, gl, teffl, corem,
                               ce.czsurf_m, 16, 8, 2)
                  + _f(ce.czsurf_r, 8, 4)
                  + _es([czsurf_t, czsurf_d, czsurf_p, czsurf_o], 12, 4)
                  + _es([tl, dl, pl, beta, eta, hcomp[1, 1], hcomp[2, 1],
                         hcomp[3, 1]] + list(tlumx[1:6])
                        + [tlumx[8], tlumx[7], tlumx[6]], 16, 8)
                  + _es([fl.clsnu, fl.gasnu] + list(fluxtot[1:11]), 10, 3)
                  + _es(_hcomp_columns(hcomp, m)
                        + [sjtot, ac.tauczn, ac.tcz, ac.tnorm, ac.tauhe,
                           ac.whe, ac.tatmos, rotv, ovrtrn_mod.taucz,
                           sss1, sss2, sss3, rrr1, rrr2, rrr3], 16, 8)
                  + f"{icheck:8d}")
    return gl, lpunch


def _shell_location(lshell, jxbeg, jxmid, jxend, hs1, hr, rl):
    """The hydrogen shell's base, middle and top, in mass (Msun) and in
    radius (fraction of R), for the track (the block both track versions
    repeat)."""
    global sss1, sss2, sss3, rrr1, rrr2, rrr3
    c = constants_mod
    if lshell:
        sss1 = hs1[jxbeg] / c.cmsun
        sss2 = hs1[jxmid] / c.cmsun
        sss3 = hs1[jxend] / c.cmsun
        rrr1 = math.exp(c.cln * (hr[jxbeg] - rl - c.crsunl))
        rrr2 = math.exp(c.cln * (hr[jxmid] - rl - c.crsunl))
        rrr3 = math.exp(c.cln * (hr[jxend] - rl - c.crsunl))
    else:
        sss1 = 0.0
        sss2 = sss1
        sss3 = sss1
        rrr1 = sss1
        rrr2 = sss1
        rrr3 = sss1
