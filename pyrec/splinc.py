"""Python version of splinc.f (YREC 7be4133): the second derivatives of a
natural cubic spline, with arrays of the fixed size JSON."""

import numpy as np

from pyrec import dimensions_mod

JSON = dimensions_mod.json

# The work array U (SAVE).
u = np.zeros(JSON + 1)


def splinc(x: np.ndarray, y: np.ndarray, y2: np.ndarray, n: int) -> None:
    """The second derivatives y2 of the natural cubic spline through
    (x, y) (y2 = 0 at both ends), for interpolation with SPLINTD2.

    Translates: splinc.f#SPLINC (7be4133)

    Arguments: x[1..n] (increasing), y[1..n]; y2[1..n], filled in place;
    n (at most JSON). None is returned.

    """

    # The natural spline's lower end (Fortran lines 16–17).
    y2[1] = 0.0
    u[1] = 0.0

    # Forward elimination (Fortran lines 18–24).
    for i in range(2, n):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i])
                       - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p

    # The natural spline's upper end (Fortran lines 25–27).
    qn = 0.0
    un = 0.0
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0)

    # Back substitution (Fortran lines 28–30).
    for k in range(n - 1, 0, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
