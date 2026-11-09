"""Python version of search.f (YREC 7be4133): a binary search in an
increasing table, for the shape-preserving quadratic splines of McAllister
and Roulier."""

import numpy as np


def search(xtab: np.ndarray, num: int, s: float, lcn: int, fnd: int) -> tuple[int, int]:
    """Locate s in xtab[1..num], with xtab[1] <= s <= xtab[num]: lcn is the
    index of the last entry below s, or of the entry equal to s (then
    fnd = 1; otherwise fnd = 0).

    Translates: search.f#SEARCH (7be4133)

    Arguments: xtab[1..num] (increasing); num; s; lcn and fnd, set only.
    Returns (lcn, fnd).

    """

    # Fortran lines 51–53.
    first = 1
    last = num
    fnd = 0

    # The ends of the table (Fortran lines 55–63); equality of reals is
    # exact (see the page).
    if xtab[1] == s:
        return 1, 1
    if xtab[num] == s:  # label 10
        return num, 1

    # Bisection (labels 20–60, Fortran lines 67–91).
    while last - first != 1:  # label 20
        middle = (first + last) // 2
        if xtab[middle] < s:
            first = middle  # label 40
        elif xtab[middle] == s:
            return middle, 1  # label 50
        else:
            last = middle  # label 60
    lcn = first  # label 30
    return lcn, fnd
