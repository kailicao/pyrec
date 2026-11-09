"""Python version of getmodel2.f (YREC 7be4133): read a stellar model in
the Model2 text format."""

from pyrec import constants_mod
import numpy as np
from pyrec import iounits_mod


def _fields(line, spec):
    """Read one record by a FORMAT, given as (kind, width[, decimals])
    items: 'a' text, 'i' integer, 'f' real (F or E), 'l' logical, 'x'
    skipped. Blanks inside numbers are ignored and a blank field reads as
    0, as Fortran does by default; a real without a decimal point has its
    last `decimals` digits taken as decimals. (The 1P scale factor of the
    FORMATs only matters for fields without an exponent, which PUTMODEL2
    does not write, and is not applied.)"""

    out, pos = [], 0
    for item in spec:
        kind, w = item[0], item[1]
        text = line[pos:pos + w] if line is not None else ''
        pos += w
        if kind == 'x':
            continue
        if kind == 'a':
            out.append(text.ljust(w))
            continue
        t = text.replace(' ', '')
        if kind == 'l':
            t = t.lstrip('.')
            out.append(t[:1].upper() == 'T')
        elif kind == 'i':
            out.append(int(t) if t else 0)
        else:
            t = t.upper().replace('D', 'E')
            if not t:
                out.append(0.0)
                continue
            mant, sep, expo = t.partition('E')
            value = float(mant)
            if '.' not in mant:
                value = value / 10.0 ** item[2]
            if sep:
                value = value * 10.0 ** int(expo)
            out.append(value)
    return out


def getmodel2(bl: float, cfenv: np.ndarray, cmixl: float, dage: float, ddage: float, ftri: float, hcomp: np.ndarray,
              hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray, hs: np.ndarray, hstot: float,
              ht: np.ndarray, iread: int, jcore: int, jenv: int, lc: np.ndarray, lexcom: bool, lrot: bool, m: int,
              model: int, omega: np.ndarray, ps: np.ndarray, rs: np.ndarray, smass: float, teffl: float,
              tlumx: np.ndarray, tril: np.ndarray, trit: np.ndarray, ts: np.ndarray, atm: str, eos: str, hik: str,
              ldify: bool, ldifz: bool, ldisk: bool, linstb: bool, ljdot0: bool, alok: str, lovstc: bool, lovste: bool,
              lovstm: bool, lpurez: bool, lsemic: bool, compmix: str, pdisk: float, tdisk: float,
              wmax: float) -> tuple[float, float, float, float, float, float, int, int, bool, bool, int, int, float,
              float, str, str, str, bool, bool, bool, bool, bool, str, bool, bool, bool, bool, bool, str, float, float,
              float]:
    """Read a Model2-format model from unit iread (rewound first): the
    header, the physics flags of the run that wrote it, the luminosities
    by source (converted to solar units if written in erg/s), the
    envelope triangle with its orientation ftri, and the shells.

    Translates: getmodel2.f#GETMODEL2 (7be4133)

    Arguments: iread (the unit); every other scalar is only set.
    Arrays changed in place: cfenv, hcomp, hd, hl, hp, hr, hs, ht, lc,
    omega, ps, rs, tlumx, tril, trit, ts.

    Returns (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv, lexcom,
    lrot, m, model, smass, teffl, atm, eos, hik, ldify, ldifz, ldisk,
    linstb, ljdot0, alok, lovstc, lovste, lovstm, lpurez, lsemic,
    compmix, pdisk, tdisk, wmax).

    """

    # Fortran line 54.
    iounits_mod.rewind(iread)

    # Header (FORMAT 10): A5,2I8,5F16.12,1PE18.10,0P2F16.12 (Fortran lines
    # 56–59).
    # Suspected bug (see the page): m is not checked against 1..json
    # before hs[m] is set.
    (atemp, model, m, smass, teffl, bl, hstot, dage, ddage, hs1,
     hsm) = _fields(iounits_mod.read(iread),
                    [('a', 5), ('i', 8), ('i', 8)] + [('f', 16, 12)] * 5
                    + [('f', 18, 10), ('f', 16, 12), ('f', 16, 12)])
    hs[1] = hs1
    hs[m] = hsm

    # The run's physics (FORMAT 30) (Fortran lines 61–67).
    vals = _fields(iounits_mod.read(iread),
                   [('i', 8), ('i', 8), ('f', 16, 10), ('x', 1), ('a', 6),
                    ('x', 1)] + [('a', 4), ('x', 1)] * 3
                   + [('l', 1), ('x', 1), ('a', 4), ('x', 1)]
                   + [('l', 1), ('x', 1)] * 11 + [('f', 18, 10)] * 3)
    (jcore, jenv, cmixl, eos, atm, alok, hik, lpurez, compmix, lexcom,
     ldify, ldifz, lsemic, lovstc, lovste, lovstm, lrot, linstb, ljdot0,
     ldisk, tdisk, pdisk, wmax) = vals

    # Luminosities by source (FORMAT 40); in erg/s if any exceeds 1e20
    # (Fortran lines 69–83).
    tlumx[1:8] = _fields(iounits_mod.read(iread),
                         [('x', 10)] + [('f', 17, 9)] * 7)
    # Question (see the page): the unit is guessed from the size.
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax > 1.0e20:
        for j in range(1, 7 + 1):
            tlumx[j] = tlumx[j] / constants_mod.clsun

    # The envelope triangle (FORMAT 50); a negative record number means
    # its orientation is reversed (Fortran lines 85–95).
    ftri = 1.0
    for i in range(1, 3 + 1):
        vals = _fields(iounits_mod.read(iread),
                       [('x', 4), ('i', 2)] + [('f', 16, 12)] * 5
                       + [('f', 20, 12)] * 3)
        io, trit[i], tril[i], ps[i], ts[i], rs[i] = vals[:6]
        for j in range(1, 3 + 1):
            cfenv[3 * i - 3 + j] = vals[5 + j]
        if io < 0:
            ftri = -1.0

    # Three header records of the shell table (Fortran lines 97–100).
    iounits_mod.read(iread)
    iounits_mod.read(iread)
    iounits_mod.read(iread)

    # The shells (FORMAT 60) (Fortran lines 102–108).
    spec60 = ([('i', 6), ('f', 18, 14), ('f', 18, 14), ('f', 24, 16)]
              + [('f', 18, 14)] * 3 + [('f', 24, 16), ('x', 1), ('l', 1)]
              + [('f', 12, 9)] * 3 + [('f', 16, 8)] * 12)
    for i in range(1, m + 1):
        vals = _fields(iounits_mod.read(iread), spec60)
        (k, hs[i], hr[i], hl[i], hp[i], ht[i], hd[i], omega[i],
         lc[i]) = vals[:9]
        for j in range(1, 15 + 1):
            hcomp[j, i] = vals[8 + j]

    return (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv, lexcom, lrot,
            m, model, smass, teffl, atm, eos, hik, ldify, ldifz, ldisk,
            linstb, ljdot0, alok, lovstc, lovste, lovstm, lpurez, lsemic,
            compmix, pdisk, tdisk, wmax)
