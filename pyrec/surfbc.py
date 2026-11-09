"""Python version of surfbc.f (YREC 7be4133): the outer boundary condition
from the envelope triangle."""

import math

import numpy as np

from pyrec import alatm05_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.envint import envint

# DATA NUMENV/0/ with SAVE: a count of the envelopes computed, kept between
# calls (not the envelope's point count in envstruct_mod) (Fortran line 77).
numenv = 0


def surfbc(trit: np.ndarray, tril: np.ndarray, cfenv: np.ndarray, ps: np.ndarray, ts: np.ndarray, rs: np.ndarray,
           ftri: float, istore: int, estore: np.ndarray, lnew: bool, lreset: bool, ksaha: int, kenv: int, katm: int,
           hstot: float, bl: float, teffl: float, x: float, z: float, fpl: float, ftl: float, lflp: bool, plim: float,
           lc: np.ndarray, m: int, cmixl: float) -> tuple[float, int, bool, bool, int, int, int, bool]:
    """Keep three envelope integrations around (log Teff, log L), flip the
    triangle when the model leaves it, and fit the plane cfenv that gives
    log R, log L and log Teff from log P and log T at the fitting point;
    switch between atmosphere tables and the gray atmosphere as Teff
    changes.

    Translates: surfbc.f#SURFBC (7be4133)

    Arguments: the triangle trit, tril (log Teff, log L of its vertices,
    index 1-3) and its envelopes' fitting-point values ps, ts, rs (log P,
    log T, log R), cfenv (the plane's 9 coefficients), ftri (orientation),
    istore, estore (passed to envint), lnew (start a new triangle), lreset
    (recompute all three envelopes), the counters ksaha, kenv, katm, hstot
    (log M), bl (log L/Lsun), teffl (log Teff), the surface x, z, the
    rotation factors fpl, ftl, plim, lc (convective flags; lc[m] for the
    spots), m, and cmixl (the mixing length, passed to envint). Arrays are
    changed in place; lflp is output only. Stops (SystemExit) when the
    triangle's reflection walk exceeds 1000 reflections or the plane fit is
    degenerate.

    Returns (ftri, istore, lnew, lreset, ksaha, kenv, katm, lflp): lflp is
    True if an envelope was recomputed. teffl is no longer changed (7be4133:
    each vertex is integrated at a local vteffl).

    """

    global numenv
    c = constants_mod
    nml = nmlpar_mod
    ltri = [None, False, False, False]

    # Above the atmosphere table's range, switch to the gray atmosphere and
    # start a new triangle (MHP 9/01) (Fortran lines 82–108).
    # Question (see the page): the limit 3.95 is hard-coded; KTTAU = 5 has
    # no switch.
    if nml.kttau == 3:
        if teffl >= 3.95:
            print()
            print(f"LOG TEFF OF {teffl:7.3f} ABOVE 3.95 - SWITCH"
                  f"TO GRAY ATMOSPHERE BOUNDARY CONDITION")  # format 5
            print()
            nml.kttau = 0
            lnew = True
            nml.lttau = True  # remember that the flag is switched (MHP 06/13)
    if nml.kttau == 4:
        if teffl >= alatm05_mod.al_tefflmax:
            print()
            print(f"LOG TEFF OF {teffl:7.3f} ABOVE Allard Table max "
                  f"{alatm05_mod.al_tefflmax:7.3f}"
                  f"  - SWITCH TO GRAY ATMOSPHERE BOUNDARY CONDITION")
            print()
            nml.kttau = 0
            lnew = True
            nml.lttau = True

    # Switch back below the limit, with a new triangle, since the old
    # vertices were integrated with the gray T(tau) (MHP 10/30/26)
    # (Fortran lines 109–129).
    if nml.lttau:
        if nml.kttau0 == 3 and teffl < 3.95:
            nml.kttau = nml.kttau0
            nml.lttau = False
            lnew = True
            print(f"LOG TEFF OF {teffl:7.3f} BELOW 3.95 - SWITCH"
                  f" BACK TO KURUCZ ATMOSPHERE BOUNDARY CONDITION")
        elif nml.kttau0 == 4 and teffl < alatm05_mod.al_tefflmax:
            nml.kttau = nml.kttau0
            nml.lttau = False
            lnew = True
            print(f"LOG TEFF OF {teffl:7.3f} below Allard Table max "
                  f"{alatm05_mod.al_tefflmax:7.3f}"
                  f"  - SWITCH BACK TO ALLARD ATMOSPHERE BOUNDARY CONDITION")

    # Which envelopes to redo (Fortran lines 130–145). Question (see the
    # page): trierr is always zero.
    lflp = False
    if lnew or lreset:
        # Redo all three envelopes.
        ltri[1] = ltri[2] = ltri[3] = False
        istore = 0
        trierr = 0.0
        lreset = False
    else:
        ltri[1] = ltri[2] = ltri[3] = True
        trierr = 0.0

    if lnew:
        # A new triangle about (teffl, bl) (Fortran lines 146–155).
        ftri = 1.0
        trit[3] = teffl
        trit[1] = trit[3] - 0.5 * nml.tridt
        trit[2] = trit[1] + nml.tridt
        tril[3] = bl + 0.5 * nml.tridl
        tril[1] = tril[3] - nml.tridl
        tril[2] = tril[1]
        lnew = False
    else:
        # Check that (teffl, bl) is inside the triangle; flip a vertex
        # across the opposite side until it is, at most 1000 times
        # (MHP 10/30/26) (Fortran lines 156–185).
        nrefl = 0
        flipped = True
        while flipped:  # label 10
            flipped = False
            for i1 in range(1, 4):  # loop 20
                i2 = i1 % 3 + 1
                i3 = i2 % 3 + 1
                temp = ftri * ((tril[i2] - tril[i3]) * (teffl - trit[i2])
                               + (trit[i3] - trit[i2]) * (bl - tril[i2]))
                if temp < -trierr:
                    ftri = -ftri
                    trit[i1] = trit[i2] + trit[i3] - trit[i1]
                    tril[i1] = tril[i2] + tril[i3] - tril[i1]
                    ltri[i1] = False
                    trierr = 0.0
                    nrefl = nrefl + 1
                    if nrefl > 1000:
                        text = (" SURFBC: ENVELOPE TRIANGLE WALK FAILED AFTER"
                                f" 1000 REFLECTIONS AT LOG TEFF{teffl:9.5f}"
                                f" LOG L{bl:9.5f}\n RUN STOPPED")  # format 60
                        print(text)
                        iounits_mod.write(iounits_mod.ishort, text)
                        raise SystemExit(1)
                    flipped = True
                    break  # GOTO 10

    # Compute new envelopes where needed (Fortran lines 186–229).
    lrenv = False
    for i in range(1, 4):  # loop 40
        if not ltri[i]:
            lflp = True
            # The vertex's log Teff, kept local (MHP 10/29/26).
            vteffl = trit[i]
            b = math.exp(c.cln * tril[i])
            rl = 0.5 * (tril[i] + c.clsunl - 4.0 * vteffl - c.c4pil
                        - c.csigl)
            gl = c.cgl + comp_mod.stotal - rl - rl
            numenv = numenv + 1
            lpenv = False  # part of LSTENV (G Somers 11/14)
            ixx = i
            lsbc0 = True
            lpulpt = False  # no pulsation output
            # For spotted stars, the ambient temperature (G Somers 10/14).
            if lc[m] and nml.spotf != 0.0 and nml.spotx != 1.0:
                ateffl = vteffl - 0.25 * math.log10(
                    nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
            else:
                ateffl = vteffl
            rl, istore, katm, kenv, ksaha = envint(
                b, fpl, ftl, gl, hstot, ixx, lpenv, lsbc0, plim, rl, ateffl,
                x, z, estore, istore, katm, kenv, ksaha, ps, rs, ts, lpulpt,
                cmixl)
            lrenv = True

    if lrenv:
        # The plane through the three envelopes, for log R, log L, log Teff;
        # a degenerate fit stops the run (MHP 10/30/26) (Fortran lines
        # 230–262).
        if abs(ts[2] - ts[1]) <= 1.0e-10 * abs(ts[3] - ts[1]):
            _degenerate_fit(ps, ts)
        temp = 1.0 / (ts[2] - ts[1])
        temp1 = (ts[3] - ts[1]) * temp
        den = ps[3] - ps[1] - temp1 * (ps[2] - ps[1])
        if abs(den) <= 1.0e-10 * (abs(ps[3] - ps[1])
                                  + abs(temp1 * (ps[2] - ps[1]))):
            _degenerate_fit(ps, ts)
        temp2 = 1.0 / den
        cfenv[1] = (rs[3] - rs[1] - temp1 * (rs[2] - rs[1])) * temp2
        cfenv[4] = (tril[3] - tril[1] - temp1 * (tril[2] - tril[1])) * temp2
        cfenv[7] = (trit[3] - trit[1] - temp1 * (trit[2] - trit[1])) * temp2
        cfenv[2] = (rs[2] - rs[1] - cfenv[1] * (ps[2] - ps[1])) * temp
        cfenv[5] = (tril[2] - tril[1] - cfenv[4] * (ps[2] - ps[1])) * temp
        cfenv[8] = (trit[2] - trit[1] - cfenv[7] * (ps[2] - ps[1])) * temp
        cfenv[3] = rs[1] - cfenv[1] * ps[1] - cfenv[2] * ts[1]
        cfenv[6] = tril[1] - cfenv[4] * ps[1] - cfenv[5] * ts[1]
        cfenv[9] = trit[1] - cfenv[7] * ps[1] - cfenv[8] * ts[1]

        # Format 50 (Fortran lines 263–268).
        labels = ["", "        LOG(R)", "     LOG(L/L0)", "       LOG(TE)"]
        text = (" ENVELOPE TRIANGLE  LOG(TE)  LOG(L)   LOG(P)   LOG(T)"
                "   LOG(R)   COEFFICIENTS")
        for i in range(1, 4):
            row = "".join(f"{v:9.5f}"
                          for v in (trit[i], tril[i], ps[i], ts[i], rs[i]))
            coef = "".join(f"{cfenv[3 * i - 3 + j]:10.5f}"
                           for j in range(1, 4))
            text += (f"\n{i:14d}{'T' if ltri[i] else 'F'}  {row}"
                     f"{labels[i]}{coef}")
        iounits_mod.write(iounits_mod.ishort, text)
    else:
        # Fortran lines 270–273.
        iounits_mod.write(iounits_mod.ishort,
                          " ENVELOPE TRIANGLE NOT CHANGED")  # format 70

    return ftri, istore, lnew, lreset, ksaha, kenv, katm, lflp


def _degenerate_fit(ps, ts):
    """Write format 61 to the screen and ISHORT, then STOP 1."""

    # Fortran lines 237–252.
    text = (" SURFBC: DEGENERATE ENVELOPE TRIANGLE FIT; FITTING POINTS"
            " (LOG T, LOG P):\n"
            + "".join(f"{ts[i]:12.7f}{ps[i]:12.7f}\n" for i in range(1, 4))
            + " RUN STOPPED")
    print(text)
    iounits_mod.write(iounits_mod.ishort, text)
    raise SystemExit(1)
