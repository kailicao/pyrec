"""Python version of getyrec7.f (YREC 09c3b9b): read a stellar model in
the older YREC7 text format."""

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod


def _fields(line, spec):
    """Read one record by a FORMAT, given as (kind, width[, decimals])
    items: 'a' text, 'i' integer, 'f' real (F or E), 'l' logical, 'x'
    skipped. Blanks inside numbers are ignored and a blank field reads as
    0, as Fortran does by default; a real without a decimal point has its
    last `decimals` digits taken as decimals. (The 1P scale factor of the
    FORMATs only matters for fields without an exponent, and is not
    applied.)"""
    out, pos = [], 0
    for item in spec:
        kind, w = item[0], item[1]
        text = line[pos:pos + w]
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


def _read_values(unit, n, spec):
    """Read n values with a FORMAT that is reused for further records
    (Fortran's format reversion); None at the end of the file."""
    values = []
    while len(values) < n:
        line = iounits_mod.read(unit)
        if line is None:
            return None
        values += _fields(line, spec)
    return values[:n]


def getyrec7(bl, cfenv, cmixl, dage, ddage, ftri, hcomp, hd, hl, hp, hr,
             hs, hstot, ht, iread, ishort, jcore, jenv, lc, lexcom, lrot, m,
             model, omega, ps, rs, smass, teffl, tlumx, tril, trit, ts, atm,
             eos, hik, ldify, ldifz, ldisk, linstb, ljdot0, alok, lovstc,
             lovste, lovstm, lpurez, lsemic, compmix, pdisk, tdisk, wmax):
    """Read a YREC7-format model from unit iread (rewound before and
    after): the header, the luminosities, the envelope triangle, the
    shells with X and Z in units of 1e-6, the other species stored only
    between the convective core and envelope, and Omega (as -log Omega)
    with rotation. The Model2-only items get defaults.

    Translates: getyrec7.f#GETYREC7 (09c3b9b)

    Arguments: iread (the unit), ishort (for the error messages); every
    other scalar is only set.
    Arrays changed in place: cfenv, hcomp, hd, hl, hp, hr, hs, ht, lc,
    omega, ps, rs, tlumx, tril, trit, ts.

    Returns (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv, lexcom,
    lrot, m, model, smass, teffl, atm, eos, hik, ldify, ldifz, ldisk,
    linstb, ljdot0, alok, lovstc, lovste, lovstm, lpurez, lsemic,
    compmix, pdisk, tdisk, wmax).
    """
    json = dimensions_mod.json
    for i in range(1, json + 1):
        omega[i] = 0.0
        for j in range(1, 15 + 1):
            hcomp[j, i] = 0.0

    # The Model2-only items, which this format lacks.
    atm = "  ? "
    eos = "    ? "
    hik = "  ? "
    alok = "  ? "
    compmix = "AG93"
    ldify = False
    ldifz = False
    ldisk = False
    linstb = False
    ljdot0 = False
    lovstc = False
    lovste = False
    lovstm = False
    lpurez = False
    lsemic = False
    # Question (see the page): defaults that change the physics.
    pdisk = 7.2722e-6
    tdisk = 0.0
    wmax = 3e-4

    iounits_mod.rewind(iread)

    # Header (FORMAT 10).
    (atemp, model, m, smass, teffl, bl, hstot, dage, ddage) = _fields(
        iounits_mod.read(iread),
        [('a', 4), ('x', 1), ('i', 5), ('i', 5)] + [('f', 8, 5)] * 3
        + [('x', 1), ('f', 14, 10), ('f', 13, 9), ('f', 13, 0)])

    # Flags (FORMAT 30).
    # Question (see the page): 5X skips one column more than PUTYREC7
    # writes.
    lrot, lexcom, jcore, jenv, cmixl = _fields(
        iounits_mod.read(iread),
        [('l', 1), ('x', 1), ('l', 1), ('x', 1), ('x', 5), ('i', 4),
         ('i', 4), ('f', 7, 3)])

    # Luminosities by source (FORMAT 40); in erg/s if any exceeds 1e20.
    tlumx[1:8] = _fields(iounits_mod.read(iread),
                         [('x', 10)] + [('f', 10, 3)] * 7)
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax > 1.0e20:
        for j in range(1, 7 + 1):
            tlumx[j] = tlumx[j] / constants_mod.clsun

    # The envelope triangle (FORMAT 70).
    ftri = 1.0
    for i in range(1, 3 + 1):  # label 80
        vals = _fields(iounits_mod.read(iread),
                       [('x', 3), ('i', 2), ('f', 7, 5)] + [('f', 8, 5)] * 4
                       + [('f', 12, 5)] * 3)
        io, trit[i], tril[i], ps[i], ts[i], rs[i] = vals[:6]
        for j in range(1, 3 + 1):
            cfenv[i + i + i - 3 + j] = vals[5 + j]
        if io < 0:
            ftri = -1.0

    # The shells (FORMAT 100), X and Z in units of 1e-6.
    spec100 = ([('f', 13, 10), ('f', 10, 7), ('f', 14, 7), ('f', 10, 7),
                ('f', 10, 7), ('f', 10, 7), ('l', 1), ('i', 6), ('i', 6)])
    for i in range(1, m + 1):  # label 110
        (hs[i], hr[i], hl[i], hp[i], ht[i], hd[i], lc[i], ix,
         iz) = _fields(iounits_mod.read(iread), spec100)
        if hs[i] < 0.0 or hs[i] > hstot:
            iounits_mod.write(ishort, ' ' + '>' * 39 + '<' * 40)
            iounits_mod.write(ishort, ' RUN STOPPED')
            iounits_mod.write(ishort, ' ERROR IN SUBROUTINE GETY7')
            iounits_mod.write(ishort, f' GLITCH IN SHELL{i:3d}, SHELL MASS '
                              'LESS THAN ZERO OR GREATER THAN STAR MASS')
            raise SystemExit
        hcomp[1, i] = 1.0e-6 * float(ix)
        hcomp[3, i] = 1.0e-6 * float(iz)

    # Species 4-11, stored from jcore to jenv only (FORMAT 200).
    spec200 = [('f', 9, 3), ('x', 1)] * 8
    if jenv <= jcore:
        hcomp[4:12, jcore] = _fields(iounits_mod.read(iread), spec200)
        hcomp[4:12, jenv] = _fields(iounits_mod.read(iread), spec200)
    else:
        for j in range(jcore, jenv + 1):  # label 210
            hcomp[4:12, j] = _fields(iounits_mod.read(iread), spec200)
    if jcore > 1:
        # The homogeneous core.
        for j in range(1, jcore):  # label 230
            for i in range(4, 11 + 1):  # label 220
                hcomp[i, j] = hcomp[i, jcore]
    if jenv < m:
        # The homogeneous envelope.
        for j in range(jenv + 1, m + 1):  # label 250
            for i in range(4, 11 + 1):  # label 240
                hcomp[i, j] = hcomp[i, jenv]
    for i in range(1, m + 1):  # label 260
        hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # The light elements (species 12-15), two shells per record.
    if lexcom:
        spec300 = [('f', 9, 3), ('x', 1)] * 4
        if jenv == 1:
            # Fully convective: one record for the whole star.
            hcomp[12:16, 1] = _fields(iounits_mod.read(iread), spec300)
            for j in range(1, m + 1):  # label 310
                for i in range(12, 15 + 1):  # label 305
                    hcomp[i, j] = hcomp[i, 1]
        else:
            for k in range(jcore, jenv, 2):  # label 315
                vals = _fields(iounits_mod.read(iread), spec200)
                hcomp[12:16, k] = vals[0:4]
                hcomp[12:16, k + 1] = vals[4:8]
            # Fortran's K after the loop: one step past its last value.
            kt = jenv - 1 - jcore
            k = jcore + 2 * max(0, int((kt + 2) / 2))
            if kt % 2 != 0:
                # An odd number of shells: the last one alone.
                hcomp[12:16, k] = _fields(iounits_mod.read(iread), spec300)
            if jcore > 1:
                for j in range(1, jcore):  # label 330
                    for i in range(12, 15 + 1):  # label 320
                        hcomp[i, j] = hcomp[i, jcore]
            if jenv < m:
                for j in range(jenv + 1, m + 1):  # label 350
                    for i in range(12, 15 + 1):  # label 340
                        hcomp[i, j] = hcomp[i, jenv]

    # Rotation (FORMAT 500): -log Omega, 58.9 or more for Omega = 0.
    if lrot:
        omegal = _read_values(iread, m, [('f', 10, 7)] * 8)
        if omegal is None:  # END=9999
            iounits_mod.write(ishort,
                              "GETYREC7: OMEGA records are missing from "
                              " input model file - OMEGA array zeroed")
            for i in range(1, json + 1):
                omega[i] = 0.0
            iounits_mod.rewind(iread)
            return (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv,
                    lexcom, lrot, m, model, smass, teffl, atm, eos, hik,
                    ldify, ldifz, ldisk, linstb, ljdot0, alok, lovstc,
                    lovste, lovstm, lpurez, lsemic, compmix, pdisk, tdisk,
                    wmax)
        for i in range(1, m + 1):  # label 510
            if omegal[i - 1] < 58.9:
                omega[i] = 10.0 ** (-omegal[i - 1])
            else:
                omega[i] = 0.0

    iounits_mod.rewind(iread)
    return (bl, cmixl, dage, ddage, ftri, hstot, jcore, jenv, lexcom, lrot,
            m, model, smass, teffl, atm, eos, hik, ldify, ldifz, ldisk,
            linstb, ljdot0, alok, lovstc, lovste, lovstm, lpurez, lsemic,
            compmix, pdisk, tdisk, wmax)
