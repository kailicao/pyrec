"""Python version of bandw.f (YREC 09c3b9b): a banded linear solver for
the angular-momentum transport in dadcoeft."""

import sys

NMAX = 8000
TINY = 1.0e-20


def _print_row(a, i):
    # FORMAT 911: I5, 1P10E12.3
    print("%5d" % i + "".join("%12.3E" % a[i, j] for j in range(1, 11)))


def bandw(a, nm, b):
    """Solve A x = b by Gaussian elimination without pivoting, for a band
    matrix with 4 subdiagonals and 5 superdiagonals: a[i, 5 + k] is the
    coefficient of x[i + k], k = -4..5 (diagonal in column 5).

    Translates: bandw.f#BANDW (09c3b9b)

    Arguments: a (NMAX+1 by 11, index 0 unused), overwritten; nm, the
    number of unknowns; b (NMAX+1), the right-hand side, replaced in place
    by the solution. A pivot below TINY stops the run (STOP 999).

    Returns nothing.
    """
    # Forward elimination, full band: rows 1..nm-5
    for i in range(1, nm - 5 + 1):
        if abs(a[i, 5]) < TINY:
            _print_row(a, i)
            sys.exit(999)
        b[i] = b[i] / a[i, 5]
        for j in range(6, 11):
            a[i, j] = a[i, j] / a[i, 5]
        for j in range(i + 1, i + 4 + 1):
            imj = i - j
            imj5 = 5 + imj          # column of x_i in row j
            b[j] = b[j] - a[j, imj5] * b[i]
            for k in range(imj5 + 1, 10 + imj + 1):
                a[j, k] = a[j, k] - a[j, imj5] * a[i, k + j - i]

    # The last rows, whose band is cut off by the end of the matrix
    for i in range(nm - 4, nm):
        ii = nm - i
        if abs(a[i, 5]) < TINY:
            _print_row(a, i)
            sys.exit(999)
        b[i] = b[i] / a[i, 5]
        for j in range(6, 5 + ii + 1):
            a[i, j] = a[i, j] / a[i, 5]
        for j in range(i + 1, nm + 1):
            imj = i - j
            imj5 = 5 + imj
            b[j] = b[j] - a[j, imj5] * b[i]
            for k in range(imj5 + 1, imj5 + ii + 1):
                a[j, k] = a[j, k] - a[j, imj5] * a[i, k + j - i]
    # The last row stops without printing
    if abs(a[nm, 5]) < TINY:
        sys.exit(999)
    b[nm] = b[nm] / a[nm, 5]

    # Back substitution (unit diagonal)
    for i in range(nm, 5, -1):
        for j in range(-1, -6, -1):
            b[i + j] = b[i + j] - b[i] * a[i + j, 5 - j]
    for i in range(5, 1, -1):
        jj = -(i - 1)
        for j in range(-1, jj - 1, -1):
            b[i + j] = b[i + j] - b[i] * a[i + j, 5 - j]
