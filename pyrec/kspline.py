"""Python version of kspline.f (YREC development, 7be4133)."""

import numpy as np

NM = 4


def kspline(x: np.ndarray, y: np.ndarray, y2: np.ndarray) -> None:
    """Second derivatives of a natural cubic spline through 4 points.

    *Numerical Recipes*' SPLINE with the size fixed at ``NM`` = 4 and
    natural end conditions (``y2[1] = y2[4] = 0``). KSPLINT evaluates it.

    Translates: kspline.f#KSPLINE (7be4133)

    Parameters
    ----------
    x, y : ndarray, shape (NM + 1,)
        The nodes and values, ``x[1..4]``, ``y[1..4]``.
    y2 : ndarray, shape (NM + 1,)
        Set in place: the second derivatives at the nodes.

    """

    # The work array U is SAVEd in the Fortran but fully reset each call.
    u = np.zeros(NM + 1)
    n = NM

    # The natural spline's lower end (Fortran lines 11–14).
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
