"""Python version of quadeos06.f (YREC 7be4133): quadratic interpolation
through three points, for the 2006 OPAL equation of state."""

import numpy as np

# The cache of abscissa factors, slots 1..30 (SAVE).
xx12 = np.zeros(31)
xx13 = np.zeros(31)
xx23 = np.zeros(31)
xx1sq = np.zeros(31)
xx1pxx2 = np.zeros(31)


def quadeos06(ic: int, i: int, x: float, y1: float, y2: float, y3: float, x1: float, x2: float, x3: float) -> float:
    """The value at x of the quadratic through (x1, y1), (x2, y2),
    (x3, y3). With ic = 0 the abscissa factors are computed and cached in
    slot i; otherwise slot i is reused.

    Translates: quadeos06.f#QUADEOS06 (7be4133)

    Returns the interpolated value.

    """

    # The abscissa factors, computed or reused (Fortran lines 22–34).
    if ic == 0:
        xx12[i] = 1.0 / (x1 - x2)
        xx13[i] = 1.0 / (x1 - x3)
        xx23[i] = 1.0 / (x2 - x3)
        xx1sq[i] = x1 * x1
        xx1pxx2[i] = x1 + x2

    # The coefficients and the value (Fortran lines 35–40).
    c3 = (y1 - y2) * xx12[i]
    c3 = c3 - (y2 - y3) * xx23[i]
    c3 = c3 * xx13[i]
    c2 = (y1 - y2) * xx12[i] - xx1pxx2[i] * c3
    c1 = y1 - x1 * c2 - xx1sq[i] * c3
    return c1 + x * (c2 + x * c3)
