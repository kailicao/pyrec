"""Python version of momi.f (YREC 09c3b9b): the moment of inertia of each
shell of a rotating, distorted model, and its derivative dI/dOmega."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec.intmom import intmom


def momi(eta2, hr, hs, hs2, jstart, jend, omega, r0, hi, qiw):
    """Compute the moment of inertia hi and dI/dOmega qiw of the shells
    jstart..jend (Law's thesis, Yale 1980, p. 61), from the shapes found by
    SHAPE; thin spherical shells when WALPCZ is not 0.

    Translates: momi.f#MOMI (09c3b9b)

    Arguments: the shape parameters eta2 and r0; log r (hr), log m (hs),
    the shell masses hs2 and the angular velocities omega; the range
    jstart..jend. Arrays have index 0 unused.

    Sets in place: hi, qiw. Returns nothing.
    """
    c = constants_mod
    # Question (see the page): r0p, rphip (and h0 below) are never used.
    if jstart == 1:
        r0p = 0.0
        rphip = 0.0
    else:
        r0p = math.log(r0[jstart - 1])
        rphip = c.cln * hr[jstart - 1]
    con = c.cc13 * 5.0 / math.exp(c.cln * c.cgl)
    # Question (see the page): WALPCZ as a switch for every shell.
    if nmlpar_mod.walpcz != 0.0:
        for i in range(jstart, jend + 1):
            hi[i] = c.cc23 * hs2[i] * math.exp(c.cln * 2.0 * hr[i])
            qiw[i] = 0.0
        if jstart == 1:
            rad = 0.5 * (hr[1] + hr[2])
            hi[1] = 0.4 * hs2[1] * math.exp(c.cln * 2.0 * rad)
        return
    for i in range(jstart, jend + 1):                    # loop 20
        # QR0R = d ln r0 / d ln r.
        qr0r = 1.0
        h0 = c.cc23 * hs2[i] * math.exp(c.cln * 2.0 * hr[i])
        r03 = r0[i]**3
        rphi3 = math.exp(c.cln * 3.0 * hr[i])
        r0r = (r03 / rphi3) * r0[i]**2
        a = con * omega[i]**2 * r03 / (math.exp(c.cln * hs[i])
                                       * (2.0 + eta2[i]))
        eta2x = eta2[i]
        # The integral at the zone centre.
        him, qiwm = intmom(a, eta2x, qr0r, r0r, 0.0, 0.0)
        hi[i] = him * hs2[i]
        qiw[i] = qiwm * hs2[i]
        r0p = math.log(r0[i])
        rphip = c.cln * hr[i]
