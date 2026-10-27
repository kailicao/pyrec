"""Python version of inter3.f (YREC 09c3b9b): 3-point Lagrange weights
and derivative weights."""


def inter3(a, b, c, x):
    """Set b to the 3-point Lagrange interpolation weights at x for the
    nodes a[1:4], and c to the weights of the first derivative.

    Translates: inter3.f#INTER3 (09c3b9b)

    Arguments: a (the nodes), x.
    Arrays changed in place: b, c (index 0 unused).

    Returns nothing.
    """
    a32 = a[3] - a[2]
    a31 = a[3] - a[1]
    a21 = a[2] - a[1]
    d1 = a21 * a31
    d2 = -a21 * a32
    d3 = a31 * a32
    xa1 = x - a[1]
    xa2 = x - a[2]
    xa3 = x - a[3]
    b[1] = (xa2 * xa3) / d1
    b[2] = (xa1 * xa3) / d2
    b[3] = (xa1 * xa2) / d3
    c[1] = (xa2 + xa3) / d1
    c[2] = (xa1 + xa3) / d2
    c[3] = (xa1 + xa2) / d3
