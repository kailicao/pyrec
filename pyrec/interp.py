"""Python version of interp.f (YREC 09c3b9b): 4-point Lagrange weights
and derivative weights."""


def interp(a, b, c, x):
    """Set b to the 4-point Lagrange interpolation weights at x for the
    nodes a[1:5], and c to the weights of the first derivative.

    Translates: interp.f#INTERP (09c3b9b)

    Arguments: a (the nodes), x.
    Arrays changed in place: b, c (index 0 unused).

    Returns nothing.
    """
    a43 = a[4] - a[3]
    a42 = a[4] - a[2]
    a41 = a[4] - a[1]
    a32 = a[3] - a[2]
    a31 = a[3] - a[1]
    a21 = a[2] - a[1]
    d1 = -a21 * a31 * a41
    d2 = a21 * a32 * a42
    d3 = -a31 * a32 * a43
    d4 = a41 * a42 * a43
    xa1 = x - a[1]
    xa2 = x - a[2]
    xa3 = x - a[3]
    xa4 = x - a[4]
    b[1] = (xa2 * xa3 * xa4) / d1
    b[2] = (xa1 * xa3 * xa4) / d2
    b[3] = (xa1 * xa2 * xa4) / d3
    b[4] = (xa1 * xa2 * xa3) / d4
    c[1] = (xa3 * xa4 + xa2 * xa4 + xa2 * xa3) / d1
    c[2] = (xa3 * xa4 + xa1 * xa4 + xa1 * xa3) / d2
    c[3] = (xa2 * xa4 + xa1 * xa4 + xa1 * xa2) / d3
    c[4] = (xa2 * xa3 + xa1 * xa3 + xa1 * xa2) / d4
