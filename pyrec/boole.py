"""Python version of boole.f (YREC development, 09c3b9b)."""

import numpy as np

from pyrec.splinj import splinj
from pyrec.splint import splint

SCALEX = 1.0e-11
SCALEY = 1.0e7


def boole(x, y, n, nn, s):
    """Integrate tabulated data with Boole's rule after spline resampling.

    The data are rescaled, fitted with a cubic spline, evaluated at ``nn``
    evenly spaced points, integrated with Boole's rule (5-point
    Newton-Cotes) and scaled back. Only use: the acoustic depth in CALCAD.

    Translates: boole.f#BOOLE (09c3b9b)

    Parameters
    ----------
    x, y : ndarray, shape (n + 1,)
        The tabulated abscissae and values, ``x[1..n]``, ``y[1..n]``.
    n : int
        Number of tabulated points.
    nn : int
        Number of points of the even grid; should be 1 + 4p.
    s : float
        Output only (Fortran declares it ``S(1)``, used as a scalar).

    Returns
    -------
    s : float
        The integral of ``y`` over ``x``.
    """
    xd = np.zeros(n + 1)
    yd = np.zeros(n + 1)
    d2x = np.zeros(n + 1)
    xi = np.zeros(nn + 1)
    yi = np.zeros(nn + 1)
    klo = khi = 0

    # Rescale radius and sound-speed vectors to values near 1 (label 22).
    for i in range(1, n + 1):
        xd[i] = x[i] * SCALEX
        yd[i] = y[i] * SCALEY

    # Interpolate onto an even grid (done for all data sets; label 20).
    splinj(xd, yd, d2x, n)
    for i in range(1, nn + 1):
        xi[i] = xd[1] + (i - 1) * (xd[n] - xd[1]) / (nn - 1)
        yi[i], klo, khi = splint(xd, yd, n, d2x, xi[i], yi[i], klo, khi)

    # How many sets of four intervals; Fortran's integer division.
    # Question (see the page): leftover points are silently dropped
    # when nn - 1 is not a multiple of 4.
    p = (nn - 1) // 4
    h = (xd[n] - xd[1]) / (nn - 1)
    s = 0.0
    for i in range(0, p):  # label 21
        s = s + 2.0 * h * (7 * yi[1 + 4 * i] + 32 * yi[2 + 4 * i]
                           + 12 * yi[3 + 4 * i] + 32 * yi[4 + 4 * i]
                           + 7 * yi[5 + 4 * i]) / 45.0

    # Rescale the result back to actual units.
    s = s / (SCALEX * SCALEY)
    return s
