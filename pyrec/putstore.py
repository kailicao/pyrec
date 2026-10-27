"""Python version of putstore.f (YREC 09c3b9b): writes a verbose model to
the .store file."""

import math

import numpy as np

from pyrec import ccout_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import roten_mod
from pyrec import scrtch_mod
from pyrec import sound_mod
from pyrec import temp2_mod
from pyrec.pindex import pindex
from pyrec.stitch import stitch
from pyrec.wrtmil import wrtmil
from pyrec.wrtmod import wrtmod

# Locals kept between calls (SAVE). ATM keeps its last value when KTTAU is
# not 0-4; JXBEG, JXEND, LSHELL are never set (see question 1).
atm = None
jxbeg = 0
jxend = 0
lshell = False

_ATM_LABELS = {0: 'EDD ', 1: 'KS  ', 2: 'HRA ', 3: 'KUR ', 4: 'ALL '}

# Output formats are approximate: Fortran's 0PE (0.ddd) is written as
# Python's E (d.ddd).
_FMT1013 = ('# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX'
            '  LEXCOM  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM'
            '  LROT  LINSTB  LJDOT0  LDISK  TDISK  PDISK  WMAX  LSTORE'
            '  LSTATM  LSTENV  LSTMOD  LSTPHYS  LSTROT')
_FMT1014 = (
    'MODEL SHELL MASS RADIUS LUMINOSITY PRESSURE TEMPERATURE DENSITY OMEGA '
    'CONVECTIVE INTERIOR_PT ENV_PT ATM_PT H1 He4 METALS He3 C12 C13 N14 N15 O16 '
    'O17 O18 H2 Li6 Li7 Be9 OPACITY GRAV DELR DEL DELAD V_CONV GAM1 HII HEII HEIII '
    'BETA ETA PPI PPII PPIII CNO TRIPLE_ALPHA E_NUC E_NEU E_GRAV CP DLNRHODLNT A '
    'RP/RE FP FT J/M MOMENT DEL_KE V_ES V_GSF V_SS VTOT ')
_FMT55 = (
    '\n'
    ' SHELL       MASS             RADIUS             LUMINOSITY            '
    'PRESSURE         TEMPERATURE         DENSITY               OMEGA      '
    '    C     H1          He4        METALS         He3             C12   '
    '          C13             N14             N15             O16         '
    '    O17             O18             H2              Li6             Li7'
    '             Be9           OPAC       GRAV        DELR        DEL      '
    '   DELA       V_CONV     GAM1      HII     HEII     HEIII    BETA      '
    'ETA       PPI         PPII       PPIII        CNO         3HE         '
    'E_NUC        E_NEU       E_GRAV          A           RP/RE           FP'
    '            FT           J/M          MOMENT        DEL_KE       V_ES '
    '      V_GSF      V_SS       VTOT   \n')


def _l(flag):
    """Fortran's L1."""
    return 'T' if flag else 'F'


def _flags_line(jcore, jenv, cmixl, eos, alok, hik, compmix):
    """FORMATs 29 and 30 (29 adds a leading '#'): the physics flags."""
    n = nmlpar_mod
    return (f"{jcore:8d}{jenv:8d}{cmixl:16.10f} {eos:6s} {atm:4s} "
            f"{alok:4s} {hik:4s} {_l(n.lpurez)} {compmix:4s} "
            + "".join(_l(v) + " " for v in (
                n.lexcom, n.ldify, n.ldifz, n.lsemic, n.lovstc, n.lovste,
                n.lovstm, n.lrot, n.linstb, n.ljdot0, n.ldisk))
            + f"{n.tdisk:18.10E}{n.pdisk:18.10E}{n.wmax:18.10E} "
            + "".join(_l(v) + " " for v in (
                n.lstore, n.lstatm, n.lstenv, n.lstmod, n.lstphys,
                n.lstrot)))


def putstore(hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps, ts, rs,
             cfenv, ftri, tlumx, jcore, jenv, model, m, smass, teffl, bl,
             hstot, dage, ddage, omega, hs1, eta2, r0, fp, ft, hjm, hi,
             cmixl, xmsol, compmix):
    """Write the model to the .store file (unit ISTOR): either the
    "stitched" format through stitch (LSTCH), or the older one: a MOD2
    header, the flags, TLUMX, the envelope lines and, with LSTMOD, a line
    per selected shell with optional physics (LSTPHYS) and rotation
    (LSTROT) columns; then, with LSTATM or LSTENV, the envelope and
    atmosphere through wrtmil (LMILNE) and wrtmod.

    Translates: putstore.f#PUTSTORE (09c3b9b)

    Arguments: the model's arrays (hcomp ... lc, omega, hs1, eta2, r0, fp,
    ft, hjm, hi; index 0 unused), the envelope triangle (trit, tril, ps,
    ts, rs, cfenv, ftri), the luminosities tlumx, and the scalars; the
    mixing length cmixl, xmsol (passed to wrtmod) and compmix, the initial
    mixture's code. Arrays changed in place: tlumx (converted to erg/s if
    given in solar units). Nothing is returned.
    """
    global atm
    n = nmlpar_mod
    c = constants_mod
    istor = iounits_mod.istor
    sc = scrtch_mod
    t2 = temp2_mod

    # Physics flags. Atmosphere flag ATM.
    # Question (see the page): KTTAU = 5 has no label.
    if n.kttau in _ATM_LABELS:
        atm = _ATM_LABELS[n.kttau]
    # Equation of state flag EOS.
    eos = 'SAHA  '
    if n.ldh:
        eos = 'SAH+DH'
    if n.lscv:
        eos = 'SCV   '
        if n.ldh:
            # Question (see the page): overwritten on the next line.
            if n.lopale06:
                eos = 'SCV+O6'
            eos = 'SCV+DH'
            if n.lopale06:
                eos = 'SCDHO6'
    else:
        if n.lopale06:
            eos = 'OPAL06'
            if n.ldh:
                eos = 'OP6+DH'
    # Low-temperature opacity flag (question 7: always NONE).
    alok = 'NONE'
    # High-temperature opacity flag.
    hik = 'NONE'
    if n.lopal95:
        hik = 'OP95'

    if atm == ' ? ':
        iounits_mod.write(iounits_mod.ishort,
                          '*** YREC7 input file, flags, etc., have been '
                          'defaulted.  ***')

    # JvS 09/25: the stitched format, or the old one, not both.
    if n.lstch:
        if ccout_mod.lphhd:
            iounits_mod.write(istor, _FMT1013)
            iounits_mod.write(istor, '#' + _flags_line(
                jcore, jenv, cmixl, eos, alok, hik, compmix))
            iounits_mod.write(istor, _FMT1014)
            ccout_mod.lphhd = False   # the physics header only once
        stitch(hcomp, hr, hp, hd, hs, ht, hl, hs1, omega, eta2, hi, r0, hjm,
               fp, ft, teffl, hstot, bl, m, lc, model, cmixl)
        return

    # Header record (FORMATs 10-13): the decimals depend on DAGE's size.
    if dage < 1.0e3:
        head = "".join(f"{v:16.11f}" for v in (smass, teffl, bl, hstot, dage))
    else:
        if dage < 1.0e4:
            nd = 10
        elif dage < 1.0e5:
            nd = 9
        else:
            nd = 8
        head = ("".join(f"{v:16.12f}" for v in (smass, teffl, bl, hstot))
                + f"{dage:16.{nd}f}")
    iounits_mod.write(istor, f"MOD2 {model:8d}{m:8d}" + head
                      + f"{ddage:18.10E}{hs[1]:16.12f}{hs[m]:16.12f}")

    # Physics flags (FORMAT 30).
    iounits_mod.write(istor, _flags_line(jcore, jenv, cmixl, eos, alok,
                                         hik, compmix))

    # Luminosities: if TLUMX are in solar units, convert to erg/s (decided
    # by comparing with 1e20).
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax <= 1.0e20:
        for j in range(1, 8):
            tlumx[j] = tlumx[j] * c.clsun
    iounits_mod.write(istor, 'TLUMX     '
                      + "".join(f"{tlumx[j]:17.9E}" for j in range(1, 8)))

    # Envelope data (FORMAT 50).
    for i in range(1, 4):
        i0 = i
        if ftri < 0.0:
            i0 = -i
        iounits_mod.write(
            istor, f"ENV{i0:2d}{trit[i]:16.12f}{tril[i]:16.12f}"
            f"{ps[i]:16.12f}{ts[i]:16.12f}{rs[i]:16.12f}"
            + "".join(f"{cfenv[3 * i - 3 + j]:20.12E}" for j in range(1, 4)))

    # Which shells to print. Question (see the page): jxbeg, jxend and
    # lshell are never set.
    id_ = np.zeros(len(hs), dtype=int)
    idm = pindex(jxbeg, jxend, lshell, m, id_, 0)

    if n.lstmod:
        iounits_mod.write(istor, _FMT55)
        cg = math.exp(c.cln * c.cgl)
        for ii in range(1, idm + 1):
            i = id_[ii]
            # The basic columns, led by MODEL (FORMAT 62).
            line = (f"{model:6d} {i:6d}{hs[i]:18.14f}{hr[i]:18.14f}"
                    f"{hl[i]:24.16E}{hp[i]:18.14f} {ht[i]:18.14f} "
                    f"{hd[i]:18.14f}{omega[i]:24.16E} {_l(lc[i])}"
                    + "".join(f"{hcomp[j, i]:12.9f}" for j in range(1, 4))
                    + "".join(f"{hcomp[j, i]:16.8E}" for j in range(4, 16))
                    + "  ")
            # Additional physics (FORMAT 63).
            if n.lstphys:
                sg = math.exp(c.cln * (c.cgl - 2.0 * hr[i])) * hs1[i]
                # Question (see the page): the ionization fractions are
                # written as zeros.
                vals = [sc.so[i], sg, sc.sdel[1, i], sc.sdel[2, i],
                        sc.sdel[3, i], sc.svel[i], sound_mod.gam1[i],
                        0.0, 0.0, 0.0, sc.sbeta[i], sc.seta[i]]
                vals += [sc.seg[k, i] for k in range(1, 6)]
                vals += [sc.sesum[i], sc.seg[6, i], sc.seg[7, i]]
            else:
                vals = [0.0] * 20
            line += (f"{vals[0]:10.4E}{vals[1]:11.3E}{vals[2]:12.4E}"
                     f"{vals[3]:12.4E}{vals[4]:12.4E}{vals[5]:12.4E}"
                     f"{vals[6]:9.5f}"
                     + "".join(f" {v:9.5f}" for v in vals[7:12])
                     + "".join(f"{v:12.4E}" for v in vals[12:17])
                     + "".join(f"{v:13.5E}" for v in vals[17:20]))
            # Additional rotation info (FORMAT 64).
            if n.lstrot and n.lrot:
                fm = math.exp(c.cln * hs[i])
                duma = c.cc13 * omega[i]**2 / (cg * fm) * 5.0 / (2.0 + eta2[i])
                a = duma * r0[i]**3
                rpoleq = (1.0 - a) / (1.0 + 0.5 * a)
                vtot = t2.ves[i] + t2.vgsf[i] + t2.vss[i]
                rot = [a, rpoleq, fp[i], ft[i], hjm[i], hi[i],
                       roten_mod.derot[i], t2.ves[i], t2.vgsf[i], t2.vss[i],
                       vtot]
            else:
                rot = [0.0] * 11
            line += ("".join(f"{v:14.6E}" for v in rot[:7])
                     + "".join(f"{v:11.3E}" for v in rot[7:]))
            iounits_mod.write(istor, line)

    # The envelope and atmosphere.
    if n.lstatm or n.lstenv:
        if n.lmilne:
            wrtmil(hcomp, hd, hl, hp, hr, hs1, m, model)
        wrtmod(m, jenv, hcomp, hs1, hd, hl, hp, hr, ht, model, bl, teffl,
               fp, ft, hs, dage, cmixl, xmsol)

    iounits_mod.write(istor, '\n\n')   # FORMAT 65, (/,/)
