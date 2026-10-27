"""Python version of putmodel2.f (YREC 09c3b9b): write a model in the
MODEL2 text format."""

import math

from pyrec import constants_mod
from pyrec import iounits_mod


def _l1(flag):
    """Fortran's L1."""
    return "T" if flag else "F"


def _e1(x, w, d):
    """Fortran's 1PEw.d: one digit before the point, e.g. 1.2345E+05."""
    return f"{x:{w}.{d}E}"


def _e0(x, w, d):
    """Fortran's 0PEw.d: a mantissa in [0.1, 1), e.g. 0.12345678E+01."""
    if x == 0.0:
        return f"{'0.' + '0' * d + 'E+00':>{w}}"
    e = math.floor(math.log10(abs(x))) + 1
    m = x / 10.0**e
    if round(abs(m), d) >= 1.0:  # rounding carried into the next decade
        m /= 10.0
        e += 1
    s = f"{abs(m):.{d}f}"[1:]
    sign = "-" if x < 0 else ""
    return f"{sign + '0' + s + f'E{e:+03d}':>{w}}"


def putmodel2(bl, cfenv, cmixl, dage, ddage, ftri, hcomp, hd, hl,
              hp, hr, hs, hstot, ht, iwrite, ishort, jcore, jenv, lc, lexcom,
              lrot, m, model, omega, ps, rs, smass, teffl, tlumx, tril, trit,
              ts, atm, eos, hik, ldify, ldifz, ldisk, linstb, ljdot0, alok,
              lovstc, lovste, lovstm, lpurez, lsemic, compmix, pdisk, tdisk,
              wmax):
    """Write the model to unit iwrite in the MODEL2 format: a header line
    (its age format chosen by the age's size), a line of physics flags,
    the luminosities TLUMX(1..7) in erg/s, three ENV lines, column titles
    and one line per shell.

    Translates: putmodel2.f#PUTMODEL2 (09c3b9b)

    Arguments: as in the Fortran (arrays with index 0 unused; trit, tril,
    ps, ts, rs [1..3], cfenv [1..9], tlumx [1..8]). tlumx is changed in
    place when converted from solar units. Returns nothing.
    """
    w = iounits_mod.write
    if atm == ' ? ':
        w(ishort, "*** YREC7 input file, flags, etc., have been "
                  "defaulted.  ***")

    # The header line; DAGE's format by its size (formats 10-13).
    if dage < 1e3:
        fages = [f"{smass:16.11f}", f"{teffl:16.11f}", f"{bl:16.11f}",
                 f"{hstot:16.11f}", f"{dage:16.11f}"]
    else:
        if dage < 1e4:
            fage = f"{dage:16.10f}"
        elif dage < 1e5:
            fage = f"{dage:16.9f}"
        else:
            fage = f"{dage:16.8f}"
        fages = [f"{smass:16.12f}", f"{teffl:16.12f}", f"{bl:16.12f}",
                 f"{hstot:16.12f}", fage]
    w(iwrite, "MOD2 " + f"{model:8d}{m:8d}" + "".join(fages)
      + _e1(ddage, 18, 10) + f"{hs[1]:16.12f}{hs[m]:16.12f}")

    # Format 30: the physics flags.
    w(iwrite, f"{jcore:8d}{jenv:8d}{cmixl:16.10f} {eos:6.6s} "
      f"{atm:4.4s} {alok:4.4s} {hik:4.4s} {_l1(lpurez)} {compmix:4.4s} "
      + "".join(_l1(f) + " " for f in (lexcom, ldify, ldifz, lsemic,
                                       lovstc, lovste, lovstm, lrot,
                                       linstb, ljdot0, ldisk))
      + _e1(tdisk, 18, 10) + _e1(pdisk, 18, 10) + _e1(wmax, 18, 10))

    # Luminosities taken to be solar if all are at most 1e20, converted in
    # place. Question (see the page): units guessed from size; TLUMX(8)
    # is not written.
    cccmax = max(tlumx[1], tlumx[2], tlumx[3], tlumx[4], tlumx[5],
                 abs(tlumx[6]), tlumx[7])
    if cccmax <= 1.0e20:
        for j in range(1, 8):
            tlumx[j] = tlumx[j] * constants_mod.clsun
    w(iwrite, "TLUMX" + " " * 5
      + "".join(_e1(tlumx[j], 17, 9) for j in range(1, 8)))

    for i in range(1, 4):
        i0 = i
        if ftri < 0.0:
            i0 = -i
        w(iwrite, f"ENV{i0:2d}{trit[i]:16.12f}{tril[i]:16.12f}"
          f"{ps[i]:16.12f}{ts[i]:16.12f}{rs[i]:16.12f}"
          + "".join(_e1(cfenv[3 * i - 3 + j], 20, 12) for j in range(1, 4)))

    w(iwrite, "")
    # Format 55: A6,2A18,A24,3A18,A24,A2,3A12,12A16.
    titles = [' SHELL', 'MASS       ', 'RADIUS      ', 'LUMINOSITY       ',
              'PRESSURE     ', 'TEMPERATURE   ', 'DENSITY     ',
              'OMEGA         ', ' C', '  H1     ', '   He4    ',
              '   METALS  ', 'He3      ', 'C12      ', 'C13      ',
              'N14      ', 'N15      ', 'O16      ', 'O17      ',
              'O18      ', 'H2       ', 'Li6      ', 'Li7      ',
              'Be9      ']
    widths = [6, 18, 18, 24, 18, 18, 18, 24, 2] + [12] * 3 + [16] * 12
    w(iwrite, "".join(f"{t:>{n}}" for t, n in zip(titles, widths)))
    w(iwrite, "")

    # Format 60: one line per shell.
    for i in range(1, m + 1):
        w(iwrite, f"{i:6d}{hs[i]:18.14f}{hr[i]:18.14f}" + _e1(hl[i], 24, 16)
          + f"{hp[i]:18.14f}{ht[i]:18.14f}{hd[i]:18.14f}"
          + _e1(omega[i], 24, 16) + " " + _l1(lc[i])
          + "".join(f"{hcomp[j, i]:12.9f}" for j in range(1, 4))
          + "".join(_e0(hcomp[j, i], 16, 8) for j in range(4, 16)))
