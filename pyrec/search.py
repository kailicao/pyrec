"""Python version of search.f (YREC 09c3b9b): a binary search in an
increasing table, for the shape-preserving quadratic splines of McAllister
and Roulier."""


def search(xtab, num, s, lcn, fnd):
    """Locate s in xtab[1..num], with xtab[1] <= s <= xtab[num]: lcn is the
    index of the last entry below s, or of the entry equal to s (then
    fnd = 1; otherwise fnd = 0).

    Translates: search.f#SEARCH (09c3b9b)

    Arguments: xtab[1..num] (increasing); num; s; lcn and fnd, set only.
    Returns (lcn, fnd).
    """
    first = 1
    last = num
    fnd = 0

    # Question (see the page): exact equality of reals.
    if xtab[1] == s:
        return 1, 1
    if xtab[num] == s:  # label 10
        return num, 1

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
