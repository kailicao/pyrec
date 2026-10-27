"""Python version of splinj.f (YREC 09c3b9b): the second derivatives of a
natural cubic spline, with arrays sized N (the algorithm of SPLINC)."""

import numpy as np


def splinj(x, y, y2, n):
    """The second derivatives y2 of the natural cubic spline through
    (x, y) (y2 = 0 at both ends), for interpolation with SPLINT.

    Translates: splinj.f#SPLINJ (09c3b9b)

    Arguments: x[1..n] (increasing), y[1..n]; y2[1..n], filled in place;
    n. None is returned.
    """
    # Question (see the page): U is an automatic array of size N.
    u = np.zeros(n + 1)
    y2[1] = 0.0
    u[1] = 0.0
    # Forward elimination.
    for i in range(2, n):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i])
                       - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
    qn = 0.0
    un = 0.0
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0)
    # Back substitution.
    for k in range(n - 1, 0, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
