"""Python version of ll95tbl.f (YREC 7be4133): read the OPAL 1995 opacity
tables for all metallicities and build the table at the run's Z."""

from pyrec import iounits_mod  # unit numbers; Python opens files by name
from pyrec import llot95a_mod
from pyrec import nmlpar_mod
from pyrec.op95ztab import op95ztab

NUMD = 19     # densities, log R = log(rho/T6**3)
NUMX = 10     # hydrogen abundances
NUMZ = 13     # metallicities
NUMXZ = 126   # tables in the file


class _EndOfFile(Exception):
    """An END= branch of a Fortran READ."""


def _record(f):
    """One record, as a READ consumes it; _EndOfFile at the end."""

    line = f.readline()
    if line == '':
        raise _EndOfFile
    return line.rstrip('\n')


def _fields(line, widths, start=0):
    """Split a fixed-format record into floats (a blank field reads as 0,
    as in Fortran)."""

    out, pos = [], start
    for w in widths:
        text = line[pos:pos + w].strip()
        out.append(float(text.replace('D', 'E')) if text else 0.0)
        pos += w
    return out


def _xz(line):
    """X and Z from a table's title line, format (36X,F7.4,11X,F7.4): the
    11 columns between them ("Y=... Z=") are skipped, not read."""

    return _fields(line, [7], start=36)[0], _fields(line, [7], start=54)[0]


def _read_row(f, n, i, ncol):
    """A row in log T: F4.2 for TGR(I), then NCOL values F7.3 into
    CAPPA2(N,I,1..NCOL) (formats 30 to 35)."""

    a = llot95a_mod
    vals = _fields(_record(f), [4] + [7] * ncol)
    a.tgr[i] = vals[0]
    a.cappa2[n, i, 1:ncol + 1] = vals[1:]


def _stop_eof(fliv95, nn):
    """Label 9999: an unexpected end of file; format 9000 to the terminal
    and the .short file, then STOP 3, the input-error code."""

    text = (f" OPAL95: UNEXPECTED END OF FILE IN {fliv95.strip()}\n"
            f" WHILE READING TABLE{nn:4d} OF{NUMXZ:4d}. RUN STOPPED")
    print(text)
    iounits_mod.write(iounits_mod.ishort, text)
    raise SystemExit(3)


def ll95tbl(fliv95: str) -> None:
    """Read the full set of OPAL95 opacity tables (log kappa against Z, X,
    log T and log R) into llot95a_mod, checking each table's X and Z
    against the grids, filling the corner without data with 9.999; then
    build the fixed-Z table with OP95ZTAB at Z = ZOPAL951.

    Translates: ll95tbl.f#LL95TBL (7be4133)

    Arguments: fliv95, the file name. Nothing is returned.

    """

    a = llot95a_mod

    # Read in the opacity tables, skipping the header. "Hey, it works."
    # (Fortran: OPEN(ILIV95,...,STATUS='OLD').) Every END= goes to label
    # 9999, which stops the run (Fortran lines 104–111).
    with open(fliv95) as f:
        nfmt = 1
        nn = 0
        try:
            while True:                                   # label 10
                dummy = _record(f)
                if dummy[nfmt - 1:nfmt + 4] == 'TABLE':
                    break
        except _EndOfFile:
            _stop_eof(fliv95, nn)                         # label 9999
        xx, zz = _xz(dummy)
        iz = 1
        ix = 1
        nn = 1
        try:
            while True:
                # Entry point for getting new tables (label 15): the slot,
                # and the check of X and Z (Fortran lines 116–127).
                # Suspected bug (see the page): the X = 1-Z tables at Z = 0.06
                # and 0.08 overwrite the X = 0 tables at Z = 0.08 and 0.10.
                n = a.nz[iz] + ix
                if ix < NUMX:
                    xxt = a.xxg[ix]
                else:
                    xxt = 1.0 - a.zzg[iz]
                if a.zzg[iz] != zz or xxt != xx:
                    iounits_mod.write(iounits_mod.ishort,
                                      ' OPAL95: Z ERROR INCOMPATIBLE TABLE')
                    raise SystemExit(3)    # STOP 3, the input-error code

                # Header: the grid in rho/T6**3 (format 20: ///,4X,19F7.1,/)
                # (Fortran lines 132–133).
                # An end of file here goes to label 9999 (the outer except).
                for _ in range(3):
                    _record(f)
                a.r0gr[1:NUMD + 1] = _fields(_record(f), [7] * NUMD,
                                             start=4)
                _record(f)

                # Full table: log kappa against log T and log R (lines
                # 136–171).
                # Question (see the page): hard-coded table shape.
                for i in range(1, 58):
                    _read_row(f, n, i, NUMD)              # format 30
                # The corner of the table without data (MHP 12/97).
                _read_row(f, n, 58, NUMD - 1)             # format 31
                a.cappa2[n, 58, 19] = 9.999
                for i in range(59, 61):
                    _read_row(f, n, i, NUMD - 2)          # format 32
                    a.cappa2[n, i, 18:20] = 9.999
                for i in range(61, 65):
                    _read_row(f, n, i, NUMD - 3)          # format 33
                    a.cappa2[n, i, 17:20] = 9.999
                for i in range(65, 70):
                    _read_row(f, n, i, NUMD - 4)          # format 34
                    a.cappa2[n, i, 16:20] = 9.999
                i = 70
                _read_row(f, n, i, NUMD - 5)              # format 35
                a.cappa2[n, i, 15:20] = 9.999

                # Exit if the correct number of tables was read in (lines
                # 174–175).
                if nn >= NUMXZ:
                    break                                 # go to 1000
                nn = nn + 1

                # Fewer Z values at the highest X; the next metallicity at
                # fixed X, or the first at the next X (Fortran lines 177–191).
                if ix <= 8:
                    nmax = NUMZ
                elif ix == 9:
                    nmax = NUMZ - 3
                else:
                    nmax = NUMZ - 1
                if iz < nmax:
                    iz = iz + 1
                else:
                    iz = 1
                    ix = ix + 1

                # The next table's X and Z (format 900: /36X,F7.4,11X,F7.4;
                # no END=, so an end of file here is a Fortran run-time
                # error, not label 9999) (Fortran lines 194–196).
                try:
                    _record(f)
                    line = _record(f)
                except _EndOfFile:
                    raise EOFError('LL95TBL: end of file in READ 900 '
                                   '(Fortran run-time error)') from None
                xx, zz = _xz(line)
        except _EndOfFile:
            _stop_eof(fliv95, nn)                         # label 9999

    # Label 1000 (Fortran lines 197–206): CLOSE(ILIV95), then a table at a fixed value of Z. (For
    # metal diffusion a 4-D interpolation in X, Z, T, rho is done instead.)
    z = nmlpar_mod.zopal951
    op95ztab(z)
