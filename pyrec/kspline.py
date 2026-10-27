"""Python version of kspline.f (YREC development, 09c3b9b)."""

import numpy as np

NM = 4


def kspline(x, y, y2):
    """Second derivatives of a natural cubic spline through 4 points.

    *Numerical Recipes*' SPLINE with the size fixed at ``NM`` = 4 and
    natural end conditions (``y2[1] = y2[4] = 0``). KSPLINT evaluates it.

    Translates: kspline.f#KSPLINE (09c3b9b)

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
    # Natural spline.
    y2[1] = 0.0
    u[1] = 0.0
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
    for k in range(n - 1, 0, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
