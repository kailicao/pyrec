"""Python version of surfbc.f (YREC 09c3b9b): the outer boundary condition
from the envelope triangle."""

import math

from pyrec import alatm05_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.envint import envint

# DATA NUMENV/0/ with SAVE: a count of the envelopes computed, kept between
# calls (not the envelope's point count in envstruct_mod).
numenv = 0


def surfbc(trit, tril, cfenv, ps, ts, rs, ftri, istore, estore, lnew,
           lreset, ksaha, kenv, katm, hstot, bl, teffl, x, z, fpl, ftl,
           lflp, plim, lc, m, cmixl):
    """Keep three envelope integrations around (log Teff, log L), flip the
    triangle when the model leaves it, and fit the plane cfenv that gives
    log R, log L and log Teff from log P and log T at the fitting point;
    switch between atmosphere tables and the gray atmosphere as Teff
    changes.

    Translates: surfbc.f#SURFBC (09c3b9b)

    Arguments: the triangle trit, tril (log Teff, log L of its vertices,
    index 1-3) and its envelopes' fitting-point values ps, ts, rs (log P,
    log T, log R), cfenv (the plane's 9 coefficients), ftri (orientation),
    istore, estore (passed to envint), lnew (start a new triangle), lreset
    (recompute all three envelopes), the counters ksaha, kenv, katm, hstot
    (log M), bl (log L/Lsun), teffl (log Teff), the surface x, z, the
    rotation factors fpl, ftl, plim, lc (convective flags; lc[m] for the
    spots), m, and cmixl (the mixing length, passed to envint). Arrays are
    changed in place; lflp is output only.

    Returns (ftri, istore, lnew, lreset, ksaha, kenv, katm, teffl, lflp):
    teffl is the last recomputed vertex's temperature if any envelope was
    recomputed; lflp is True if one was.
    """
    global numenv
    c = constants_mod
    nml = nmlpar_mod
    ltri = [None, False, False, False]

    # Above the atmosphere table's range, switch to the gray atmosphere and
    # start a new triangle (MHP 9/01).
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
    # Switch back below the limit.
    # Question (see the page): switching back does not set lnew, so the
    # gray envelopes are kept.
    if nml.lttau:
        if nml.kttau0 == 3 and teffl < 3.95:
            nml.kttau = nml.kttau0
            nml.lttau = False
            print(f"LOG TEFF OF {teffl:7.3f} BELOW 3.95 - SWITCH"
                  f" BACK TO KURUCZ ATMOSPHERE BOUNDARY CONDITION")
        elif nml.kttau0 == 4 and teffl < alatm05_mod.al_tefflmax:
            nml.kttau = nml.kttau0
            nml.lttau = False
            print(f"LOG TEFF OF {teffl:7.3f} below Allard Table max "
                  f"{alatm05_mod.al_tefflmax:7.3f}"
                  f"  - SWITCH BACK TO ALLARD ATMOSPHERE BOUNDARY CONDITION")

    lflp = False
    # Question (see the page): trierr is always zero.
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
        # A new triangle about (teffl, bl).
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
        # across the opposite side until it is.
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
                    flipped = True
                    break  # GOTO 10

    # Compute new envelopes where needed.
    lrenv = False
    for i in range(1, 4):  # loop 40
        if not ltri[i]:
            lflp = True
            # Question (see the page): teffl, an argument, is overwritten.
            teffl = trit[i]
            b = math.exp(c.cln * tril[i])
            rl = 0.5 * (tril[i] + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
            gl = c.cgl + comp_mod.stotal - rl - rl
            numenv = numenv + 1
            lpenv = False  # part of LSTENV (G Somers 11/14)
            ixx = i
            lsbc0 = True
            lpulpt = False  # no pulsation output
            # For spotted stars, the ambient temperature (G Somers 10/14).
            if lc[m] and nml.spotf != 0.0 and nml.spotx != 1.0:
                ateffl = teffl - 0.25 * math.log10(
                    nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
            else:
                ateffl = teffl
            rl, istore, katm, kenv, ksaha = envint(
                b, fpl, ftl, gl, hstot, ixx, lpenv, lsbc0, plim, rl, ateffl,
                x, z, estore, istore, katm, kenv, ksaha, ps, rs, ts, lpulpt,
                cmixl)
            lrenv = True

    if lrenv:
        # The plane through the three envelopes, for log R, log L, log Teff.
        temp = 1.0 / (ts[2] - ts[1])
        temp1 = (ts[3] - ts[1]) * temp
        temp2 = 1.0 / (ps[3] - ps[1] - temp1 * (ps[2] - ps[1]))
        cfenv[1] = (rs[3] - rs[1] - temp1 * (rs[2] - rs[1])) * temp2
        cfenv[4] = (tril[3] - tril[1] - temp1 * (tril[2] - tril[1])) * temp2
        cfenv[7] = (trit[3] - trit[1] - temp1 * (trit[2] - trit[1])) * temp2
        cfenv[2] = (rs[2] - rs[1] - cfenv[1] * (ps[2] - ps[1])) * temp
        cfenv[5] = (tril[2] - tril[1] - cfenv[4] * (ps[2] - ps[1])) * temp
        cfenv[8] = (trit[2] - trit[1] - cfenv[7] * (ps[2] - ps[1])) * temp
        cfenv[3] = rs[1] - cfenv[1] * ps[1] - cfenv[2] * ts[1]
        cfenv[6] = tril[1] - cfenv[4] * ps[1] - cfenv[5] * ts[1]
        cfenv[9] = trit[1] - cfenv[7] * ps[1] - cfenv[8] * ts[1]
        # Format 50.
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
        iounits_mod.write(iounits_mod.ishort,
                          " ENVELOPE TRIANGLE NOT CHANGED")  # format 70
    return ftri, istore, lnew, lreset, ksaha, kenv, katm, teffl, lflp
