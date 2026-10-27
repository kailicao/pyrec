"""Python version of ludcmp.f (YREC development, 09c3b9b)."""

import numpy as np

NMAX = 100
TINY = 1.0e-20


def ludcmp(a, n, np_, indx, d):
    """LU decomposition with partial pivoting (Crout, implicit scaling).

    *Numerical Recipes*' LUDCMP: ``a`` is replaced by L (below the
    diagonal, unit diagonal implied) and U; ``indx`` records the row
    swaps and ``d`` their parity. LUBKSB then solves systems with it.

    Translates: ludcmp.f#LUDCMP (09c3b9b)

    Parameters
    ----------
    a : ndarray, shape (np_ + 1, np_ + 1)
        In: the matrix, ``a[i, j]`` for ``A(I, J)``. Out, in place: its
        LU decomposition.
    n : int
        The order of the matrix.
    np_ : int
        The declared size of ``a`` (Fortran's ``NP``).
    indx : ndarray of int, shape (n + 1,)
        Set in place: the row permutation.
    d : float
        Output only (passed in with any value).

    Returns
    -------
    d : float
        +1 or -1, the parity of the row interchanges.
    """
    # The scale vector VV has a fixed size NMAX, and n is not checked
    # against it (question on the page).
    vv = np.zeros(NMAX + 1)
    imax = 0
    d = 1.0
    for i in range(1, n + 1):  # loop 12
        aamax = 0.0
        for j in range(1, n + 1):  # loop 11
            if abs(a[i, j]) > aamax:
                aamax = abs(a[i, j])
        if aamax == 0.0:
            raise SystemExit("Singular matrix.")
        vv[i] = 1.0 / aamax
    for j in range(1, n + 1):  # loop 19
        if j > 1:
            for i in range(1, j):  # loop 14
                sum_ = a[i, j]
                if i > 1:
                    for k in range(1, i):  # loop 13
                        sum_ = sum_ - a[i, k] * a[k, j]
                    a[i, j] = sum_
        aamax = 0.0
        for i in range(j, n + 1):  # loop 16
            sum_ = a[i, j]
            if j > 1:
                for k in range(1, j):  # loop 15
                    sum_ = sum_ - a[i, k] * a[k, j]
                a[i, j] = sum_
            dum = vv[i] * abs(sum_)
            if dum >= aamax:
                imax = i
                aamax = dum
        if j != imax:
            for k in range(1, n + 1):  # loop 17
                dum = a[imax, k]
                a[imax, k] = a[j, k]
                a[j, k] = dum
            d = -d
            vv[imax] = vv[j]
        indx[j] = imax
        if j != n:
            if a[j, j] == 0.0:
                a[j, j] = TINY
            dum = 1.0 / a[j, j]
            for i in range(j + 1, n + 1):  # loop 18
                a[i, j] = a[i, j] * dum
    if a[n, n] == 0.0:
        a[n, n] = TINY
    return d
