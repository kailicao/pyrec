"""Python version of putyrec7.f (YREC 7be4133): write a model in the older
YREC7 text format."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod

# The local array OMEGAL (SAVE).
omegal = np.zeros(dimensions_mod.json + 1)


def _l1(flag: bool) -> str:
    """Fortran's L1."""

    return "T" if flag else "F"


def _e1(x: float, w: int, d: int) -> str:
    """Fortran's 1PEw.d: one digit before the point, e.g. 1.234E+05."""

    return f"{x:{w}.{d}E}"


def putyrec7(bl: float, cfenv: np.ndarray, cmixl: float, dage: float, ddage: float, ftri: float, hcomp: np.ndarray,
             hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray, hs: np.ndarray, hstot: float,
             ht: np.ndarray, iwrite: int, jcore: int, jenv: int, lc: np.ndarray, lexcom: bool, lrot: bool, m: int,
             model: int, omega: np.ndarray, ps: np.ndarray, rs: np.ndarray, smass: float, teffl: float,
             tlumx: np.ndarray, tril: np.ndarray, trit: np.ndarray, ts: np.ndarray) -> None:
    """Write the model to unit iwrite in the YREC7 format: header, flags,
    TLUMX (erg/s), three ENV lines, one line per shell (X and Z as
    integers in units of 1e-6), the abundances 4-11 from jcore to jenv,
    with lexcom the abundances 12-15, with lrot |log10 Omega|.

    Translates: putyrec7.f#PUTYREC7 (7be4133)

    Arguments: as in the Fortran (arrays with index 0 unused; trit, tril,
    ps, ts, rs [1..3], cfenv [1..9], tlumx [1..8]). tlumx (converted from
    solar units) and omega (floored at 1e-59) are changed in place.
    Returns nothing.

    """

    w = iounits_mod.write

    # The header record (Fortran lines 56–62).
    head = (f"NMOD {model:5d}{m:5d}{smass:8.5f}{teffl:8.5f}{bl:8.5f} "
            f"{hstot:14.10f}{dage:13.9f}")
    if ddage >= 1e6:
        w(iwrite, head + f"{ddage:#13.0f}")     # F13.0 keeps the point
    else:
        w(iwrite, head + _e1(ddage, 13, 6))

    # The flags (Fortran lines 65–66). Suspected bug (see the page):
    # GETYREC7 reads this record with 5X, not 4.
    w(iwrite, f"{_l1(lrot)} {_l1(lexcom)} F F {jcore:4d}{jenv:4d}"
      f"{cmixl:7.3f}")

    # The luminosities, converted to erg/s in the caller's array if in
    # solar units (Fortran lines 71–79).
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax <= 1.0e20:
        for j in range(1, 8):
            tlumx[j] = tlumx[j] * constants_mod.clsun
    w(iwrite, "TLUMX" + " " * 5
      + "".join(_e1(tlumx[j], 10, 3) for j in range(1, 8)))

    # The envelope data (Fortran lines 82–88).
    for i in range(1, 4):
        i0 = i
        if ftri < 0.0:
            i0 = -i
        w(iwrite, f"ENV{i0:2d}{trit[i]:7.5f}{tril[i]:8.5f}{ps[i]:8.5f}"
          f"{ts[i]:8.5f}{rs[i]:8.5f}"
          + "".join(_e1(cfenv[3 * i - 3 + j], 12, 5) for j in range(1, 4)))

    # One line per shell; X and Z as integers in units of 1e-6
    # (Fortran lines 91–97). Question (see the page): the precision.
    for i in range(1, m + 1):
        ix = int(1.0e6 * hcomp[1, i] + 0.50)
        iz = int(1.0e6 * hcomp[3, i] + 0.50)
        w(iwrite, f"{hs[i]:13.10f}{hr[i]:10.7f}" + _e1(hl[i], 14, 7)
          + f"{hp[i]:10.7f}{ht[i]:10.7f}{hd[i]:10.7f}{_l1(lc[i])}"
          f"{ix:6d}{iz:6d}")

    def fmt200(values: list) -> str:
        # FORMAT(8(1PE9.3,1X)); more than 8 values continue on a new line.
        vals = [_e1(v, 9, 3) + " " for v in values]
        return "\n".join("".join(vals[n:n + 8]) for n in range(0, len(vals), 8))

    # The abundances 4-11 from jcore to jenv; two lines for a fully
    # convective star, as GETYREC7 reads (Fortran lines 102–110).
    if jenv <= jcore:
        w(iwrite, fmt200([hcomp[i, jcore] for i in range(4, 12)]))
        w(iwrite, fmt200([hcomp[i, jenv] for i in range(4, 12)]))
    else:
        for j in range(jcore, jenv + 1):
            w(iwrite, fmt200([hcomp[i, j] for i in range(4, 12)]))

    # With lexcom, the abundances 12-15, two shells per line (Fortran
    # lines 113–126).
    if lexcom:
        if jenv <= jcore:
            w(iwrite, fmt200([hcomp[i, jenv] for i in range(12, 16)]))
        else:
            for k in range(jcore, jenv, 2):
                w(iwrite, fmt200([hcomp[i, j] for j in (k, k + 1)
                                  for i in range(12, 16)]))
            # Fortran's K after the loop: one step past the last value.
            ntrip = max(0, (jenv - 1 - jcore) // 2 + 1)
            k = jcore + 2 * ntrip
            ktest = jenv - 1 - jcore
            if ktest % 2 != 0:
                w(iwrite, fmt200([hcomp[i, k] for i in range(12, 16)]))

    # With lrot, |log10 Omega|, eight per line; omega is floored in the
    # caller's array (Fortran lines 129–143). Question (see the page): the
    # sign of log Omega is lost.
    if lrot:
        for i in range(1, m + 1):
            if omega[i] <= 1.0e-59:
                omega[i] = 1.0e-59
            omegal[i] = abs(math.log10(omega[i]))
        n1 = m // 8
        n2 = m - n1 * 8
        j1 = 1
        for i in range(1, n1 + 1):
            w(iwrite, "".join(f"{omegal[j]:10.7f}" for j in range(j1, j1 + 8)))
            j1 = j1 + 8
        if n2 > 0:
            w(iwrite, "".join(f"{omegal[j]:10.7f}"
                              for j in range(j1, j1 + n2)))
