"""Python version of readalex06.f (YREC 09c3b9b): reads the Ferguson et
al. (2005) low-temperature opacity tables ("Alexander 2006")."""

import numpy as np

from pyrec import alot06_mod
from pyrec import alot06all_mod
from pyrec import galot06_mod
from pyrec import iounits_mod

NUMX = 9
NUMZ = 16
NUMT = 85
NUMD = 19


def _fields(line, widths, start=0):
    """Split a fixed-format record into floats (a blank field reads as 0,
    as in Fortran)."""
    out, pos = [], start
    for w in widths:
        text = line[pos:pos + w].strip()
        out.append(float(text) if text else 0.0)
        pos += w
    return out


def _read_header(unit):
    """FORMAT 10, (30X,F9.6,7X,F9.6//): X and Z, then two records skipped."""
    xx, _, zz = _fields(iounits_mod.read(unit), [9, 7, 9], start=30)
    iounits_mod.read(unit)
    iounits_mod.read(unit)
    return xx, zz


def _read_table(unit, jj):
    """Check the log R row (only its first 16 values) and read the NUMT
    rows of table JJ, in descending order of log T, checking each log T."""
    g = galot06_mod
    # FORMAT 20, (6X,19F7.3).
    r0grtab = np.zeros(NUMD + 1)
    r0grtab[1:] = _fields(iounits_mod.read(unit), [7] * NUMD, start=6)
    # Question (see the page): only 16 of the 19 densities are checked.
    for kk in range(1, 17):
        if r0grtab[kk] != g.r0gr06[kk]:
            print(f" ERROR IN ALEXANDER OPACITY TABLES:\n"
                  f" EXPECTED AND ACTUAL R{kk:3d}{g.r0gr06[kk]:7.3f}"
                  f"{r0grtab[kk]:7.3f} RUN STOPPED")
            raise SystemExit
    # Opacity at each log T (FORMAT 30, (F5.3,1X,19F7.3)).
    for j in range(NUMT, 0, -1):
        vals = _fields(iounits_mod.read(unit), [5, 1] + [7] * NUMD)
        t = vals[0]
        alot06all_mod.captalex06[jj, j, 1:NUMD + 1] = vals[2:]
        if t != g.tgr06[j]:
            print(f" ERROR IN ALEXANDER OPACITY TABLES:\n"
                  f" EXPECTED AND ACTUAL T{j:3d}{g.tgr06[j]:7.3f}{t:7.3f}"
                  f" RUN STOPPED")
            raise SystemExit


def _fmt15(xexp, xx, zexp, zz):
    return (" ERROR IN ALEXANDER OPACITY TABLES:\n"
            f" EXPECTED AND ACTUAL X,Z{xexp:7.2f}{xx:7.2f}{zexp:7.2f}"
            f"{zz:7.2f} RUN STOPPED")


def readalex06(falex06):
    """Read all 143 compositions of the low-temperature opacity tables into
    alot06all_mod.captalex06, checking each table's X, Z, log R and log T
    against the grids of galot06_mod, then set alot06_mod.capalex06 to the
    first table (X = 0, Z = 0).

    Translates: readalex06.f#READALEX06 (09c3b9b)

    Arguments: falex06, the file name (opened on unit IALEX06). Nothing is
    returned; the tables are stored in alot06all_mod and alot06_mod.
    """
    g = galot06_mod
    unit = iounits_mod.ialex06
    # Question (see the page): opened without STATUS='OLD'.
    iounits_mod.open(unit, falex06)
    # The eight X values below 1-Z, each with all 16 values of Z.
    for i in range(1, NUMX):
        for ii in range(1, NUMZ + 1):
            xx, zz = _read_header(unit)
            if xx != g.xxg06[i] or zz != g.zzg06[ii]:
                print(_fmt15(g.xxg06[i], xx, g.zzg06[ii], zz))
                # Question (see the page): the STOP here is commented out.
            jj = (i - 1) * NUMZ + ii
            _read_table(unit, jj)
    # The final tables, at X = 1-Z; none for Z = 0.1 (read at X = 0.9).
    for ii in range(1, NUMZ):
        xx, zz = _read_header(unit)
        xtest = 1.0 - g.zzg06[ii]
        if xx != xtest or zz != g.zzg06[ii]:
            print(_fmt15(xtest, xx, g.zzg06[ii], zz))
            raise SystemExit
        jj = (NUMX - 1) * NUMZ + ii
        _read_table(unit, jj)
    iounits_mod.close(unit)
    # Initialize the fixed (X, Z) table.
    for i in range(1, NUMT + 1):
        for j in range(1, NUMD + 1):
            alot06_mod.capalex06[i, j] = alot06all_mod.captalex06[1, i, j]
