"""Python version of getmodel2.f (YREC 09c3b9b): read a stellar model in
the Model2 text format."""

from pyrec import constants_mod
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


def getmodel2(bl, cfenv, cmixl, dage, ddage, ftri, hcomp, hd, hl, hp, hr,
              hs, hstot, ht, iread, jcore, jenv, lc, lexcom, lrot, m, model,
              omega, ps, rs, smass, teffl, tlumx, tril, trit, ts, atm, eos,
              hik, ldify, ldifz, ldisk, linstb, ljdot0, alok, lovstc, lovste,
              lovstm, lpurez, lsemic, compmix, pdisk, tdisk, wmax):
    """Read a Model2-format model from unit iread (rewound first): the
    header, the physics flags of the run that wrote it, the luminosities
    by source (converted to solar units if written in erg/s), the
    envelope triangle with its orientation ftri, and the shells.

    Translates: getmodel2.f#GETMODEL2 (09c3b9b)

    Arguments: iread (the unit); every other scalar is only set.
    Arrays changed in place: cfenv, hcomp, hd, hl, hp, hr, hs, ht, lc,
    omega, ps, rs, tlumx, tril, trit, ts.

    Returns (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv, lexcom,
    lrot, m, model, smass, teffl, atm, eos, hik, ldify, ldifz, ldisk,
    linstb, ljdot0, alok, lovstc, lovste, lovstm, lpurez, lsemic,
    compmix, pdisk, tdisk, wmax).
    """
    iounits_mod.rewind(iread)

    # Header (FORMAT 10): A5,2I8,5F16.12,1PE18.10,0P2F16.12.
    (atemp, model, m, smass, teffl, bl, hstot, dage, ddage, hs1,
     hsm) = _fields(iounits_mod.read(iread),
                    [('a', 5), ('i', 8), ('i', 8)] + [('f', 16, 12)] * 5
                    + [('f', 18, 10), ('f', 16, 12), ('f', 16, 12)])
    hs[1] = hs1
    hs[m] = hsm

    # The run's physics (FORMAT 30).
    vals = _fields(iounits_mod.read(iread),
                   [('i', 8), ('i', 8), ('f', 16, 10), ('x', 1), ('a', 6),
                    ('x', 1)] + [('a', 4), ('x', 1)] * 3
                   + [('l', 1), ('x', 1), ('a', 4), ('x', 1)]
                   + [('l', 1), ('x', 1)] * 11 + [('f', 18, 10)] * 3)
    (jcore, jenv, cmixl, eos, atm, alok, hik, lpurez, compmix, lexcom,
     ldify, ldifz, lsemic, lovstc, lovste, lovstm, lrot, linstb, ljdot0,
     ldisk, tdisk, pdisk, wmax) = vals

    # Luminosities by source (FORMAT 40); in erg/s if any exceeds 1e20.
    tlumx[1:8] = _fields(iounits_mod.read(iread),
                         [('x', 10)] + [('f', 17, 9)] * 7)
    # Question (see the page): the unit is guessed from the size.
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax > 1.0e20:
        for j in range(1, 7 + 1):
            tlumx[j] = tlumx[j] / constants_mod.clsun

    # The envelope triangle (FORMAT 50); a negative record number means
    # its orientation is reversed.
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

    # Three header records of the shell table.
    iounits_mod.read(iread)
    iounits_mod.read(iread)
    iounits_mod.read(iread)

    # The shells (FORMAT 60).
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
