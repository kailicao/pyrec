"""Python version of ratext.f (YREC 7be4133): rational-function
extrapolation to zero step size, for the Bulirsch-Stoer integrator."""

import numpy as np

# Fortran line 11.
IMAX = 11
NMAX = 15
NCOL = 7

# The tableau and the abscissae, kept between calls (SAVE): the routine
# is stateful, and a call with iest = 1 starts afresh (Fortran lines 12–13).
x = np.zeros(IMAX + 1)
d = np.zeros((NMAX + 1, NCOL + 1))
fx = np.zeros(NCOL + 1)


def ratext(iest: int, xest: float, yest: np.ndarray, yz: np.ndarray, dy: np.ndarray, nv: int, nuse: int) -> None:
    """Add the iest-th estimate yest (at xest, the squared substep) and
    evaluate the next diagonal of the rational extrapolation tableau,
    with at most nuse earlier estimates (Numerical Recipes' RZEXTR).

    Translates: ratext.f#RATEXT (7be4133)

    Arguments: iest; xest; yest[1..nv]; yz, dy (set in place: the
    extrapolated values and their error estimates); nv; nuse.
    Returns nothing.

    """

    # Suspected bug (see the page of bsstep.f): no checks of nv, iest,
    # nuse against NMAX, IMAX, NCOL.
    # Save the independent variable; the first estimate starts the tableau
    # (Fortran lines 16–23).
    x[iest] = xest
    if iest == 1:
        for j in range(1, nv + 1):
            yz[j] = yest[j]
            d[j, 1] = yest[j]
            dy[j] = yest[j]
    else:
        # Use at most nuse earlier estimates (Fortran lines 25–29).
        m1 = min(iest, nuse)
        for k in range(1, m1):
            fx[k + 1] = x[iest - k] / xest

        # The next diagonal of the tableau (Fortran lines 30–53).
        for j in range(1, nv + 1):
            yy = yest[j]
            v = d[j, 1]
            c = yy
            d[j, 1] = yy
            for k in range(2, m1 + 1):
                b1 = fx[k] * v
                b = b1 - c
                if b != 0.0:
                    b = (c - v) / b
                    ddy = c * b
                    c = b1 * b
                else:
                    ddy = v
                v = d[j, k]
                d[j, k] = ddy
                yy = yy + ddy
            dy[j] = ddy
            yz[j] = yy
