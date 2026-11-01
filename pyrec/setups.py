"""Python version of setups.f (YREC 09c3b9b): the one-time setup of a run's
constants and tables."""

import math

from pyrec import atmos2_mod
from pyrec import atmos2c_mod
from pyrec import atmos_mod
from pyrec import ccr_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import fac_mod
from pyrec import iounits_mod  # unit numbers; Python opens files by name
from pyrec import nmlpar_mod
from pyrec import scveos_mod
from pyrec.alfilein import alfilein
from pyrec.esac06 import seteos06
from pyrec.hra import hra
from pyrec.setupopac import setupopac

NT, NG = 57, 11      # the atmosphere grid (KTTAU = 3, 4)
NTC, NGC = 76, 11    # the Kurucz/Castelli grid (KTTAU = 5)
NTS = 63             # temperatures of the SCV tables


def _fields(line, widths):
    """Split a fixed-format record into floats (a blank field reads as 0,
    as in Fortran)."""
    out, pos = [], 0
    for w in widths:
        text = line[pos:pos + w].strip()
        out.append(float(text.replace('D', 'E')) if text else 0.0)
        pos += w
    return out


def _read_e16(f, n, rows):
    """Read n values in E16.8 fields, as many per record as `rows` says
    (formats 201, 202, 203, 206, 208)."""
    values = []
    for count in rows:
        values += _fields(f.readline(), [16] * count)
    return values[:n]


def setups(falex06, fallard, fatm, ffermi, fliv95, fpurez, fscvh, fscvhe,
           fscvz):
    """Set the physical constants, then read every table the run needs:
    opacities (through setupopac), the Fermi-Dirac tables, the
    model-atmosphere pressure tables, the OPAL 2006 EOS tables and the SCV
    EOS tables. All results go into modules; nothing is returned.

    Translates: setups.f#SETUPS (09c3b9b)

    Arguments: the file names of the tables (strings).
    """
    c = constants_mod

    # Setup constants (cgs).
    c.clndp = math.log(10.0)
    c.cln = c.clndp
    c.clni = 1.0 / c.cln
    c.clsunl = math.log10(c.clsun)          # clsun as already set
    c.clnsun = c.cln / c.clsun
    c.cmsun = 1.9891e33
    c.cmsunl = math.log10(c.cmsun)
    c.crsunl = math.log10(c.crsun)          # crsun as already set
    c.cmbol = 4.79
    c.csecyr = 3.1558e7
    c.cc13 = 1.0 / 3.0
    c.cc23 = c.cc13 + c.cc13
    c.cpi = 3.1415926535898
    c.c4pi = 4.0 * c.cpi
    c.c4pil = math.log10(c.c4pi)
    c.c4pi3l = math.log10(c.cc13 * c.c4pi)
    cc = 2.99792458e10                      # speed of light (local)
    c.csig = 5.67051e-5
    c.csigl = math.log10(c.csig)
    c.ca3 = c.csig * 4.0 / (3.0 * cc)       # radiation constant / 3
    c.ca3l = math.log10(c.ca3)
    c.cgas = 8.314510e7
    c.cgl = -7.17571                        # log G, G = 6.67259e-8
    cmelec = 9.1093897e-28
    cboltz = 1.380658e-16
    cplnck = 6.6260755e-27
    c.cmkh = (1.50 * math.log10(2.0 * c.cpi * cmelec)
              + 2.5 * math.log10(cboltz) - 3.0 * math.log10(cplnck))
    c.cdelrl = -c.c4pi3l - c.csigl - math.log10(16.0)
    c.cmixl2 = c.cc13
    c.cmixl3 = 16.0 * math.sqrt(2.0) * c.csig
    # The Debye-Hueckel constant CDH (DBG).
    hmass = 1.0 / 6.0222137e23
    elec = 4.802e-10
    debhu_mod.cdh = (-math.sqrt(c.cpi / (cboltz * hmass * hmass * hmass))
                     * elec * elec * elec / 3.0)
    # Charges of C, N, O, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, Ca, Ti, Cr, Mn,
    # Fe, Ni.
    for k, z in enumerate([6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18,
                           20, 22, 24, 25, 26, 28], start=1):
        debhu_mod.dhnue[k] = float(z)

    # Temperature at tau = 2/3 of the HRA.
    atmos_mod.hras = hra(c.cc23)

    # Opacity tables.
    setupopac(comp_mod.xenv, falex06, fliv95, fpurez)

    # F-tables for the degenerate equation of state.
    with open(ffermi) as f:
        for arr in (ccr_mod.dval, ccr_mod.etat8):
            vals = []
            while len(vals) < 43:                 # format 150: 8X,9F8.4
                line = f.readline()
                vals += _fields(line[8:], [8] * min(9, 43 - len(vals)))
            arr[1:44] = vals
        for ncards in range(1, 861):              # format 160: 2I5,5F12.7
            line = f.readline()
            # Question (see the page): I and J come from the record and are
            # used as subscripts at once.
            i, j = int(line[0:5]), int(line[5:10])
            ccr_mod.dtable[1:6, i, j] = _fields(line[10:], [12] * 5)
    k1 = 1
    for i in range(1, 42):
        # The table's density grid must increase ("should be INT(DVAL, etc.)") (MHP 10/02).
        if ccr_mod.dval[i + 1] <= ccr_mod.dval[i]:
            print(' ' + '>' * 39 + '<' * 40)      # format 1000, to ISHORT
            print(' ERROR IN SUBROUTINE SETUPS')
            print(f' GLITCH IN FERMI TABLE ELEMENT{i:4d}')
            print(' RUN STOPPED')
            raise SystemExit
        k = int((ccr_mod.dval[i + 1] - ccr_mod.dval[i]) * 20.0 + 0.10)
        k2 = k1 + k - 1
        if i == 41:
            k2 = 261
        for j in range(k1, k2 + 1):
            ccr_mod.iden[j] = i
        k1 = k2 + 1

    # Pressure tables for the surface boundary condition.
    kttau = nmlpar_mod.kttau
    if kttau == 3 or kttau == 4:
        a = atmos2_mod
        with open(fatm) as f:
            a.atmz = _fields(f.readline(), [11, 16])[1]   # format 200
            f.readline()
            a.atmtl[1:NT + 1] = _read_e16(f, NT, [4] * 14 + [1])  # 201
            f.readline()
            a.atmgl[1:NG + 1] = _read_e16(f, NG, [4, 4, 3])       # 202
            f.readline()
            for j in range(1, NT + 1):
                a.atmpl[j, 1:NG + 1] = _read_e16(f, NG, [4, 4, 3])  # 203
        _set_imin(a.atmpl, NT, NG, fac_mod.imin, fac_mod.iminmax, NG)
        # Allard model atmospheres, below 10,000 K (MHP 6/97).
        if kttau == 4:
            alfilein(fallard)
    elif kttau == 5:
        # Kurucz/Castelli atmospheres (JNT 6/14).
        a = atmos2c_mod
        with open(fatm) as f:
            atmos2_mod.atmz = _fields(f.readline(), [11, 16])[1]  # 200
            f.readline()
            a.atmtlc[1:NTC + 1] = _read_e16(f, NTC, [4] * 19)     # 206
            f.readline()
            a.atmglc[1:NGC + 1] = _read_e16(f, NGC, [4, 4, 3])    # 208
            f.readline()
            for j in range(1, NTC + 1):
                a.atmplc[j, 1:NGC + 1] = _read_e16(f, NGC, [4, 4, 3])
        # Question (see the page): an empty row sets IMIN2 to NG, not NGC.
        _set_imin(a.atmplc, NTC, NGC, fac_mod.imin2, fac_mod.iminmax2, NG)

    # OPAL 2006 EOS table.
    if nmlpar_mod.lopale06:
        seteos06()

    # SCV equation-of-state tables (MHP 5/97).
    if nmlpar_mod.lscv:
        s = scveos_mod
        with open(fscvh) as fh, open(fscvhe) as fhe, open(fscvz) as fz:
            for i in range(1, NTS + 1):
                line = fh.readline()                      # format 1
                s.tlogx[i] = _fields(line, [5])[0]
                s.nptsx[i] = int(line[5:9])
                fhe.readline()                            # dum1, dum2
                fz.readline()                             # dum3, dum4
                widths2 = [6, 13, 13] + [9] * 8           # format 2
                for j in range(1, s.nptsx[i] + 1):
                    s.tablex[i, j, 1:12] = _fields(fh.readline(), widths2)
                    s.tabley[i, j, 1:12] = _fields(fhe.readline(), widths2)
                # Metal table (Prather EOS, old Yale code), in reverse order.
                for j in range(s.nptsx[i], 0, -1):        # format 3
                    s.tablez[i, j, 1:14] = _fields(fz.readline(),
                                                   [6] + [9] * 12)


def _set_imin(atmpl, nt, ng, imin, iminmax, ng_empty):
    """G Somers 5/15: for each temperature j, the range of gravities with a
    valid pressure: iminmax[j], the highest, and imin[j], the lowest of the
    contiguous valid run below it (or ng_empty if there is none). Written
    out twice in the Fortran, for each grid."""
    for j in range(1, nt + 1):
        lcatch = False
        for k in range(ng, 0, -1):
            if atmpl[j, k] <= 0.0:
                if lcatch:
                    imin[j] = k + 1
                    break  # to label 5 (6)
            else:
                if not lcatch:
                    iminmax[j] = k
                lcatch = True
        else:
            imin[j] = 1
        if not lcatch:
            imin[j] = ng_empty
