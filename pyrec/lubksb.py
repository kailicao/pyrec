"""Python version of lubksb.f (YREC development, 09c3b9b)."""


def lubksb(a, n, np_, indx, b):
    """Solve A x = b from the LU decomposition made by LUDCMP.

    *Numerical Recipes*' LUBKSB: forward substitution with L (unscrambling
    the row permutation ``indx`` and skipping the leading zeros of b), then
    back-substitution with U. The solution overwrites ``b``.

    Translates: lubksb.f#LUBKSB (09c3b9b)

    Parameters
    ----------
    a : ndarray, shape (np_ + 1, np_ + 1)
        The LU decomposition from LUDCMP, ``a[i, j]`` for ``A(I, J)``.
    n : int
        The order of the system.
    np_ : int
        The declared size of ``a`` (Fortran's ``NP``; renamed so as not
        to hide NumPy's usual name).
    indx : ndarray of int, shape (n + 1,)
        The row permutation from LUDCMP.
    b : ndarray, shape (n + 1,)
        In: the right-hand side. Out, in place: the solution.
    """
    ii = 0
    for i in range(1, n + 1):  # loop 12
        ll = indx[i]
        sum_ = b[ll]
        b[ll] = b[i]
        if ii != 0:
            for j in range(ii, i):  # loop 11
                sum_ = sum_ - a[i, j] * b[j]
        elif sum_ != 0.0:
            ii = i
        b[i] = sum_
    for i in range(n, 0, -1):  # loop 14
        sum_ = b[i]
        if i < n:
            for j in range(i + 1, n + 1):  # loop 13
                sum_ = sum_ - a[i, j] * b[j]
        b[i] = sum_ / a[i, i]
