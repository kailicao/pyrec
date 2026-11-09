"""Python version of tridiag_gs.f (YREC 7be4133): the general tridiagonal
solver (Thomas algorithm), used by gravitational settling and the Li run."""

import numpy as np

from pyrec import dimensions_mod

# The Fortran's local work array (SAVE; every element read is set first).
gama = np.zeros(dimensions_mod.json + 1)


def tridiag_gs(a: np.ndarray, b: np.ndarray, c: np.ndarray, ex_prime: np.ndarray, npt: int, ex: np.ndarray) -> None:
    """Solve a_j x_{j-1} + b_j x_j + c_j x_{j+1} = ex_prime_j, j = 1..npt,
    by forward elimination and back substitution, without pivoting.

    Translates: tridiag_gs.f#TRIDIAG_GS (7be4133)

    Arguments: the diagonals a, b, c and the right-hand side ex_prime
    (arrays, index 0 unused); npt, the number of points; ex, the solution,
    set in place. None is returned.

    """

    # Forward elimination (Fortran lines 17–24).
    bet = b[1]
    ex[1] = ex_prime[1] / bet
    for j in range(2, npt + 1):          # loop 1
        gama[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gama[j]
        if bet == 0:
            raise SystemExit('#TRIDIA:SINGULAR MATRIX')
        ex[j] = (ex_prime[j] - a[j] * ex[j - 1]) / bet

    # Back substitution (Fortran lines 25–27).
    for j in range(npt - 1, 0, -1):      # loop 2
        ex[j] = ex[j] - gama[j + 1] * ex[j + 1]
