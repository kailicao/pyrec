"""Python version of locate.f (YREC development, 7be4133)."""

import numpy as np


def locate(xx: np.ndarray, n: int, x: float, j: int) -> int:
    """Bisection search in a monotonic table, with 1 % end tolerances.

    Finds ``j`` such that ``x`` lies between ``xx[j]`` and ``xx[j+1]``
    (ascending or descending table); ``j = 0`` or ``n`` means out of
    range, except that values within 1 % of the end values are pulled
    back in (``j = 1`` or ``n - 1``).

    Translates: locate.f#LOCATE (7be4133)

    Parameters
    ----------
    xx : ndarray, shape (n + 1,)
        The table, ``xx[1..n]``.
    n : int
        Its length.
    x : float
        The value to locate.
    j : int
        Output only (passed in with any value).

    Returns
    -------
    j : int
        The lower index of the bracketing interval.

    """

    # Bisection (Fortran lines 18–29).
    jl = 0
    ju = n + 1
    while ju - jl > 1:  # label 10
        jm = (ju + jl) // 2
        if (xx[n] > xx[1]) == (x > xx[jm]):
            jl = jm
        else:
            ju = jm
    j = jl

    # The 1 % end tolerances (Fortran lines 30–35).
    # Question (see the page): the tolerances are relative to the end
    # values, so they act only for positive ends of an ascending table.
    if j == 0 and x > 0.99 * xx[1]:
        j = 1
    if j == n and x < 1.01 * xx[n]:
        j = n - 1
    return j
