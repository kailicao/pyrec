"""Python version of splinj.f (YREC 7be4133): the second derivatives of a
natural cubic spline, with arrays sized N (the algorithm of SPLINC)."""

import numpy as np


def splinj(x: np.ndarray, y: np.ndarray, y2: np.ndarray, n: int) -> None:
    """The second derivatives y2 of the natural cubic spline through
    (x, y) (y2 = 0 at both ends), for interpolation with SPLINT.

    Translates: splinj.f#SPLINJ (7be4133)

    Arguments: x[1..n] (increasing), y[1..n]; y2[1..n], filled in place;
    n. None is returned.

    """

    # U is an automatic array of size N (Fortran line 10).
    u = np.zeros(n + 1)

    # The natural spline's lower end (Fortran lines 13–14).
    y2[1] = 0.0
    u[1] = 0.0

    # Forward elimination (Fortran lines 15–21).
    for i in range(2, n):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i])
                       - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p

    # The natural spline's upper end (Fortran lines 22–24).
    qn = 0.0
    un = 0.0
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0)

    # Back substitution (Fortran lines 25–27).
    for k in range(n - 1, 0, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
