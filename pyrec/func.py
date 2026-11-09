"""Python version of func.f (YREC 7be4133): the effective gravity and
surface element at one colatitude of a distorted equipotential."""

import math
import numpy as np

from pyrec import constants_mod


def func(x: float, g: float, s: float, r0: np.ndarray, hs: np.ndarray, aint: float, q: float, w2: float, a: float,
         i: int) -> tuple[float, float]:
    """At colatitude x on the equipotential of shell i, r = r0 (1 - a P2),
    compute the effective gravity g = |grad Phi| and the surface element
    s, the integrands of QGAUSS's angular averages.

    Translates: func.f#FUNC (7be4133)

    Arguments: x (theta), r0 (mean radii of the equipotentials), hs
    (log m), aint and q (the distortion integral and its integrand, from
    TRAPZD), w2 (Omega**2), a (the distortion amplitude), i (the shell);
    g and s are only set.

    Returns (g, s).

    """

    c = constants_mod

    # The angle's functions (Fortran lines 26–29).
    cs = math.cos(x) * math.sin(x)
    ss = math.sin(x) ** 2
    p2 = 1.0 - 1.5 * ss

    # The derivatives of r0 with respect to r and theta, and r on the
    # equipotential (Fortran lines 30–36).
    qr0r = 1.0 / (1.0 - 4.0 * a * p2)
    qr0t = 3.0 * cs * a * qr0r / (1.0 - a * p2)
    r = r0[i] * (1.0 - a * p2)
    r2 = r ** 2
    const = c.c4pi * c.cc13 / (r2 * r2)

    # The radial and angular derivatives of the potential, the gravity and
    # the surface element (Fortran lines 37–47).
    # Question (see the page): the comment assumes dr0/dr = 1.
    qphir = (math.exp(c.cln * (c.cgl + hs[i])) / r2
             - const * p2 * (3.0 * aint - r * q * qr0r) - w2 * r * ss)
    qphith = const * (p2 * q * qr0t - 3.0 * aint * cs) - w2 * r * cs
    g = math.sqrt(qphir ** 2 + qphith ** 2)
    s = (c.c4pi * r * math.sin(x) * r0[i]
         * math.sqrt((1.0 - a * p2) ** 2 + (3.0 * cs * a) ** 2))
    return g, s
