"""Python version of cspline.f (YREC 09c3b9b): the second derivatives of a
cubic spline (Numerical Recipes' SPLINE), for the opacity tables."""

import numpy as np

from pyrec import dimensions_mod

JSON = dimensions_mod.json

# The work array U (SAVE).
# Question (see the page): its fixed size JSON; N > JSON overflows it.
u = np.zeros(JSON + 1)


def cspline(x, y, n, yp1, ypn, y2):
    """The second derivatives y2 of the cubic spline through (x, y), with
    end slopes yp1, ypn; a slope of 0.99e30 or more gives a natural
    spline (zero second derivative) at that end.

    Translates: cspline.f#CSPLINE (09c3b9b)

    Arguments: x[1..n] (increasing), y[1..n]; n; yp1, ypn; y2[1..n],
    filled in place. None is returned.
    """
    if yp1 > 0.99e30:
        y2[1] = 0.0
        u[1] = 0.0
    else:
        y2[1] = -0.5
        u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1)
    # Forward elimination.
    for i in range(2, n):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i])
                       - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))
                / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p
    if ypn > 0.99e30:
        qn = 0.0
        un = 0.0
    else:
        qn = 0.5
        un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1])
                                          / (x[n] - x[n - 1]))
    y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0)
    # Back substitution.
    for k in range(n - 1, 0, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]
