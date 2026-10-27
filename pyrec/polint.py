"""Python version of polint.f (YREC 09c3b9b): polynomial interpolation by
Neville's algorithm (Numerical Recipes' POLINT)."""

import numpy as np

# Local work arrays (SAVE), fixed at 20 points.
c = np.zeros(21)
d = np.zeros(21)


def polint(xa, ya, n, x, y, dy):
    """Interpolate the polynomial through the n points (xa, ya) at x,
    starting from the nearest tabulated point.

    Translates: polint.f#POLINT (09c3b9b)

    Arguments: xa[1..n], ya[1..n] (n <= 20); x; y, dy (set). Returns
    (y, dy): the value at x and the last correction added, an error
    estimate.
    """
    ns = 1
    dif = abs(x - xa[1])
    for i in range(1, n + 1):
        dift = abs(x - xa[i])
        if dift < dif:
            ns = i
            dif = dift
        c[i] = ya[i]
        d[i] = ya[i]
    y = ya[ns]
    ns = ns - 1
    for j in range(1, n):
        for i in range(1, n - j + 1):
            ho = xa[i] - x
            hp = xa[i + j] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            if abs(den) < 1.0e-20:
                # Question (see the page): a stop without a message.
                raise SystemExit()
            den = w / den
            d[i] = hp * den
            c[i] = ho * den
        if 2 * ns < n - j:
            dy = c[ns + 1]
        else:
            dy = d[ns]
            ns = ns - 1
        y = y + dy
    return y, dy
