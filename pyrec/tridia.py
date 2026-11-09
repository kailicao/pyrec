"""Python version of tridia.f (YREC 7be4133): the tridiagonal solver of the
implicit angular-momentum diffusion step, with each point's change of
angular momentum."""

import numpy as np

from pyrec import dimensions_mod
from pyrec import tridi_mod

# The Fortran's local array of starting angular velocities (SAVE).
u0 = np.zeros(dimensions_mod.json + 1)


def tridia(n: int, ei: np.ndarray, dj: np.ndarray, sumdj: float) -> float:
    """Solve the system a, b, c, d of tridi_mod (set by dcoeft) for the new
    angular velocities tridi_mod.u, and return each point's change of
    angular momentum dj (the surface's including the wind term, which
    dcoeft leaves in tridi_mod.gama[n]) and their sum.

    Translates: tridia.f#TRIDIA (7be4133)

    Arguments: n, the number of points; ei, the moments of inertia; dj, set
    in place; sumdj is passed in with any value.

    Returns sumdj.

    """

    a = tridi_mod.a
    b = tridi_mod.b
    c = tridi_mod.c
    d = tridi_mod.d
    u = tridi_mod.u
    gama = tridi_mod.gama

    # The wind term, read from gama[n] before the elimination reuses gama
    # as work space; the starting angular velocities (Fortran lines 22–25).
    dj[n] = gama[n]
    for i in range(1, n + 1):
        u0[i] = d[i]

    # Forward elimination (Fortran lines 26–33).
    bet = b[1]
    u[1] = d[1] / bet
    for j in range(2, n + 1):          # loop 1
        gama[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gama[j]
        if bet == 0:
            raise SystemExit('#TRIDIA:SINGULAR MATRIX')
        u[j] = (d[j] - a[j] * u[j - 1]) / bet

    # The surface's change of angular momentum, wind included (Fortran
    # lines 34–35).
    dj[n] = dj[n] + (u[n] - u0[n]) * ei[n]
    sumdj = dj[n]

    # Question (see the page): a no-op, u[n] unchanged up to rounding
    # (Fortran lines 36–37).
    fj = 1.0 + (u[n] - u0[n]) / u0[n]
    u[n] = u0[n] * fj

    # Back substitution, and each point's change of angular momentum
    # (Fortran lines 38–42).
    for j in range(n - 1, 0, -1):      # loop 2
        u[j] = u[j] - gama[j + 1] * u[j + 1]
        dj[j] = (u[j] - u0[j]) * ei[j]
        sumdj = sumdj + dj[j]

    return sumdj
