"""Python version of alfilein.f (YREC 09c3b9b): reads a table of Allard
model atmospheres (the 1999 NextGen format or the newer "Allard" format)
into ALATM01_MOD, then has ALTABINIT check it."""

import math

import numpy as np

from pyrec import alatm01_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.altabinit import altabinit

NTA = 250
NGA = 25

# SAVE: locals kept between calls.
teffs = np.zeros(NTA + 1)
latmtptau100 = False  # Question (see the page): never set in the Fortran.


class _ReadError(Exception):
    """A Fortran READ that would take its ERR= branch."""


def _fields(line, specs):
    """Read a fixed-format record: specs is a list of (width, decimals).
    A blank field reads as 0; a field without a decimal point or exponent
    gets the format's implied decimals, as in Fortran; D exponents read as
    E. A field that cannot be read raises _ReadError."""
    out, pos = [], 0
    for width, decimals in specs:
        text = line[pos:pos + width].strip().upper().replace('D', 'E')
        pos += width
        if not text:
            out.append(0.0)
            continue
        try:
            value = float(text)
        except ValueError:
            raise _ReadError(text) from None
        if '.' not in text and 'E' not in text:
            value /= 10.0 ** decimals
        out.append(value)
    return out


def _read(unit):
    """One record; None at end of file (END=)."""
    return iounits_mod.read(unit)


def _error_exit():
    """label 9999: the error exit."""
    msg = " **************** PROGRAM ALFilein TERMINATED @ 9999 "
    print("")
    print("")
    print(msg + "**************")
    print("")
    print("")
    print("")
    w = iounits_mod.write
    w(iounits_mod.ishort, "")
    w(iounits_mod.ishort, msg + "*************")
    w(iounits_mod.ishort, "")
    w(iounits_mod.ishort, "")
    raise SystemExit(9999)


def _differs(feh, alpha):
    """The passes' test for a record to skip: [Fe/H] off by more than 1e-5,
    or alpha off at all."""
    n = nmlpar_mod
    # Question (see the page): two tolerances for [Fe/H] and alpha.
    return abs(feh - n.alatm_feh) > 1e-5 or abs(alpha - n.alatm_alpha) > 0.0


def _read_old(ioatma):
    """The old 1999 NextGen format (label 100): fixed sizes, log P only."""
    a = alatm01_mod
    ishort = iounits_mod.ishort
    a.loldnextgen = True
    iounits_mod.write(ishort, " ALFileIn: File Description: 1999 NEXTGEN"
                              " (Old NEXTGEN)")
    # Ensure that PL, TL at tau = 100 are not requested: they are not in
    # the old NextGen files.
    if latmtptau100:  # Question (see the page): the misspelled flag.
        text = (" ALFileIN: Invalid old Allard input file for requested"
                " PT,TL at Tau=100")
        print("")
        print(text)
        iounits_mod.write(ishort, "")
        iounits_mod.write(ishort, text)
        _error_exit()
    a.nteff = 54
    a.ngl = 5
    # FORMAT 905, (5F7.2/): the gravities, then one record skipped.
    a.gls[1:a.ngl + 1] = _fields(_read(ioatma), [(7, 2)] * 5)
    _read(ioatma)
    # FORMAT 907, (10(1P5E16.8,/),1P4E16.8,/): the 54 temperatures, as
    # given (Question (see the page): log Teff or Teff?).
    values = []
    for _ in range(10):
        values += _fields(_read(ioatma), [(16, 8)] * 5)
    values += _fields(_read(ioatma), [(16, 8)] * 4)
    a.teffls[1:a.nteff + 1] = values
    _read(ioatma)
    # FORMAT 909, (1P5E16.8): one record of nGL values per Teff.
    for j in range(1, a.nteff + 1):
        a.pls[j, 1:a.ngl + 1] = _fields(_read(ioatma), [(16, 8)] * a.ngl)
    iounits_mod.close(ioatma)


def _first_pass(ioatma):
    """Pass 1 of the new format (labels 299 to 316): collect the distinct
    Teff, log g, [Fe/H] and alpha of the records with the requested [Fe/H]
    and alpha, up to the end of the file. Raises _ReadError on a bad
    record."""
    a = alatm01_mod
    spec915 = [(6, 0), (6, 2), (6, 2), (6, 2)]  # FORMAT 915, (F6.0,3F6.2)
    # label 299: up to the first acceptable record.
    while True:
        rec = _read(ioatma)
        if rec is None:
            return
        teff, gl, feh, alpha = _fields(rec, spec915)
        if not _differs(feh, alpha):
            break
    a.nteff = 1
    a.ngl = 1
    a.nfeh = 1
    a.nalpha = 1
    teffs[a.nteff] = teff
    a.gls[a.ngl] = gl
    a.fehs[a.nfeh] = feh
    a.alphas[a.nalpha] = alpha
    # label 300: the rest of the file.
    while True:
        rec = _read(ioatma)
        if rec is None:
            return
        teff, gl, feh, alpha = _fields(rec, spec915)
        if _differs(feh, alpha):
            continue  # label 316
        if not any(abs(teff - teffs[i]) < 1e-6
                   for i in range(1, a.nteff + 1)):
            a.nteff += 1
            teffs[a.nteff] = teff
        # label 310
        if not any(abs(gl - a.gls[i]) < 1e-6 for i in range(1, a.ngl + 1)):
            a.ngl += 1
            a.gls[a.ngl] = gl
        # label 312. Question (see the page): the loop runs to nGL, not nFeH.
        if not any(abs(feh - a.fehs[i]) < 1e-6 for i in range(1, a.ngl + 1)):
            a.nfeh += 1
            a.fehs[a.nfeh] = feh
        # label 314. Question (see the page): to nGL, not nALPHA.
        if not any(abs(alpha - a.alphas[i]) < 1e-6
                   for i in range(1, a.ngl + 1)):
            a.nalpha += 1
            a.alphas[a.nalpha] = alpha
        # label 316


def _second_pass(ioatma, nhdr):
    """Pass 2 of the new format (labels 410 to 440): place each accepted
    record's log P at Teff, log P and log T at tau = 100 at its indices.
    Raises _ReadError on a bad record."""
    a = alatm01_mod
    n = nmlpar_mod
    iounits_mod.rewind(ioatma)
    for _ in range(nhdr):
        _read(ioatma)
    # FORMAT 920, (F6.0,3F6.2,1P4D16.8).
    spec920 = [(6, 0), (6, 2), (6, 2), (6, 2), (16, 8), (16, 8), (16, 8)]
    while True:  # label 410
        rec = _read(ioatma)
        if rec is None:
            return
        teff, gl, feh, alpha, pl, p100l, t100 = _fields(rec, spec920)
        if _differs(feh, alpha):
            continue  # label 440
        for i in range(1, a.nteff + 1):
            if abs(teff - teffs[i]) < 1e-6:
                i1 = i
                break
        else:
            print(" ALFilein: Impossible failure #1")
            print(f" Teff of  {teff} not in Teff table.")
            _error_exit()
        # label 420
        for j in range(1, a.ngl + 1):
            if abs(gl - a.gls[j]) < 1e-6:
                j1 = j
                break
        else:
            print(" ALFilein: Impossible failure #2")
            print(f" GL of  {gl} not in GL table.")
            _error_exit()
        # label 430: wrong [Fe/H] or alpha, bypass the record.
        if abs(feh - n.alatm_feh) >= 1e-6:
            continue
        if abs(alpha - n.alatm_alpha) >= 1e-6:
            continue
        a.teffls[i1] = math.log10(teff)  # log Teff from here on
        a.gls[j1] = gl
        a.pls[i1, j1] = pl
        a.p100ls[i1, j1] = p100l
        a.t100ls[i1, j1] = t100
        # label 440


def alfilein(fallard):
    """Read the Allard atmosphere file FALLARD (old or new format) on unit
    IOATMA, keep the records with the requested [Fe/H] and alpha
    (NMLPAR_MOD), fill ALATM01_MOD and call ALTABINIT. Entries absent from
    the file keep -999. Stops the run on any error.

    Translates: alfilein.f#ALFILEIN (09c3b9b)

    Argument: the file name fallard. Nothing is returned.
    """
    a = alatm01_mod
    ioatma = iounits_mod.ioatma
    ishort = iounits_mod.ishort

    # Every entry starts invalid (-999).
    for i in range(1, NTA + 1):
        a.teffls[i] = -999.0
        for j in range(1, NGA + 1):
            a.gls[j] = -999.0
            a.pls[i, j] = -999.0
            a.p100ls[i, j] = -999.0
            a.t100ls[i, j] = -999.0

    try:
        iounits_mod.open(ioatma, fallard, status='OLD')
    except OSError:
        # label 899
        print("")
        print(" ALFILEIN: Allard File OPEN Failure")
        print(f"      Unit-Number, FIle-Name are:  {ioatma} {fallard}")
        _error_exit()

    # label 900: the first record decides the format (FORMAT 901, (a),
    # into a CHARACTER*35).
    inrec = (_read(ioatma) or '')[:35].ljust(35)
    iounits_mod.rewind(ioatma)
    if inrec[3:7] == '3.50':
        _read_old(ioatma)  # old NextGen file type
    elif inrec[3:9] == 'ALLARD':
        # label 200: the new Allard file type.
        a.loldnextgen = False
        # FORMAT 911, (i2,X,a): the number of header records to skip.
        rec = _read(ioatma)
        nhdr = int(rec[0:2])
        inrec2 = rec[3:259].ljust(256)
        iounits_mod.write(ishort, "ALFilein: New Allard Atm: File Description: "
                          + inrec2[0:47])  # 912
        if nhdr > 1:
            for _ in range(2, nhdr + 1):
                _read(ioatma)
        a.ngl = 0
        a.nteff = 0
        a.nfeh = 0
        a.nalpha = 0
        try:
            _first_pass(ioatma)
        except _ReadError:
            # label 399
            print("")
            print(" AtTabInit Pass 1 File Error ")
            print("")
            _error_exit()
        # label 400: sort the values found.
        sort_shell(a.nteff, teffs)
        sort_shell(a.ngl, a.gls)
        sort_shell(a.nfeh, a.fehs)
        sort_shell(a.nalpha, a.alphas)
        try:
            _second_pass(ioatma, nhdr)
        except _ReadError:
            # label 499
            print("")
            print(" AtTabInit File Read error")
            print("")
            _error_exit()
    else:
        for show in (lambda t: iounits_mod.write(ishort, t), print):
            show(" *** Invalid Allard Atmosphere input file ***")
            show(f" First record was '{inrec}'")
        _error_exit()

    # label 500
    altabinit()


def sort_shell(nn, a):
    """Shell's sort (Numerical Recipes), ascending, in place, on a[1..nn],
    with the increments 1, 4, 13, 40, ...

    Translates: alfilein.f#SORT_SHELL (09c3b9b)

    Arguments: the count nn and the array a (changed in place). Nothing is
    returned.
    """
    if nn == 1:
        return
    n = nn
    inc = 1
    while True:  # label 1: the starting increment
        inc = 3 * inc + 1
        if inc > n:
            break
    while True:  # label 2: loop over the partial sorts
        inc = inc // 3
        for i in range(inc + 1, n + 1):
            v = a[i]  # straight insertion
            j = i
            while a[j - inc] > v:  # label 3
                a[j] = a[j - inc]
                j = j - inc
                if j <= inc:
                    break
            a[j] = v  # label 4
        if inc <= 1:
            break
