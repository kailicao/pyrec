"""Python version of rdzlaol.f (YREC 09c3b9b): read the LAOL89 opacity
table for pure metals (the "pure Z" table)."""

from pyrec import iounits_mod
from pyrec import zlaol_mod


def _e_records(unit, n, per=6, w=12):
    """Read n reals as FORMAT (1P6E12.5) does: per to a record, so
    ceil(n/per) records. (1P would scale a field without an exponent by
    1/10; the tables' fields have exponents, and it is not applied.)"""
    out = []
    while len(out) < n:
        line = iounits_mod.read(unit)
        for k in range(min(per, n - len(out))):
            text = line[k * w:(k + 1) * w].replace(' ', '')
            out.append(float(text.replace('D', 'E')) if text else 0.0)
    return out


def _skip(unit, nrec):
    for _ in range(nrec):
        iounits_mod.read(unit)


def rdzlaol(fpurez):
    """Read the pure-metal opacity table from the file fpurez (unit
    IOPUREZ) into ZLAOL_MOD: the sizes NUMRHOZ, NUMTZ, the grids ZORHO,
    ZOT and the opacities ZOLAOL(rho, T). The header (mix, hydrogen
    fraction) is read and dropped. Stops the run if the file has more
    than one composition or the grids exceed 104 x 52.

    Translates: rdzlaol.f#RDZLAOL (09c3b9b)

    Arguments: fpurez, the file name. Returns nothing.
    """
    z = zlaol_mod
    iopurez = iounits_mod.iopurez
    iounits_mod.open(iopurez, fpurez, status='OLD')
    # FORMAT 100: (/,18X,I2,9X,I3,14X,I3), the second record.
    _skip(iopurez, 1)
    line = iounits_mod.read(iopurez)
    n = int(line[18:20])
    z.numrhoz = int(line[29:32])
    z.numtz = int(line[46:49])
    if n != 1 or z.numrhoz > 104 or z.numtz > 52:
        iounits_mod.write(iounits_mod.ishort, '  Z OPACITY INPUT ERROR.')
        raise SystemExit()
    # FORMAT 120: (47X,F8.5,//,1P6E12.5,/,1P4E12.5), the mix: 4 records.
    _skip(iopurez, 4)
    # FORMAT 131: (54X,F8.5,///////), the hydrogen fraction: 8 records.
    _skip(iopurez, 8)
    # FORMAT 140: (/,(1P6E12.5)), one dummy value: 2 records.
    _skip(iopurez, 2)
    # FORMATs 150 and 160: (/,(1P6E12.5)), the grids.
    _skip(iopurez, 1)
    vals = _e_records(iopurez, z.numrhoz)
    for ii in range(1, z.numrhoz + 1):
        z.zorho[ii] = vals[ii - 1]
    _skip(iopurez, 1)
    vals = _e_records(iopurez, z.numtz)
    for ii in range(1, z.numtz + 1):
        z.zot[ii] = vals[ii - 1]
    _skip(iopurez, 1)                    # FORMAT 170
    for ir in range(1, z.numrhoz + 1):
        _skip(iopurez, 1)                # FORMAT 200
        vals = _e_records(iopurez, z.numtz)     # FORMAT 210
        for it in range(1, z.numtz + 1):
            z.zolaol[ir, it] = vals[it - 1]
    iounits_mod.close(iopurez)
