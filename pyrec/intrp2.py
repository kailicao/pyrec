"""Python version of intrp2.f (YREC development, 09c3b9b)."""


def intrp2(a, b, x):
    """Weights of 4-point Lagrange interpolation at ``x``, no derivatives.

    The same as INTERP without its derivative weights. Used by the opacity
    tables, the diffusion regridding and PHYSIC.

    Translates: intrp2.f#INTRP2 (09c3b9b)

    Parameters
    ----------
    a : ndarray, shape (5,)
        The four nodes, ``a[1..4]``.
    b : ndarray, shape (5,)
        Set in place: the weights ``b[1..4]``, with
        ``f(x) ~ sum_i b[i] f(a[i])``.
    x : float
        The point of interpolation.
    """
    # SAVE in the Fortran is harmless: every local is set before use.
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
