"""Python version of polint.f (YREC 7be4133): polynomial interpolation by
Neville's algorithm (Numerical Recipes' POLINT)."""

import numpy as np

from pyrec import iounits_mod

# Local work arrays (SAVE), fixed at 20 points.
c = np.zeros(21)
d = np.zeros(21)


def polint(xa: np.ndarray, ya: np.ndarray, n: int, x: float, y: float, dy: float) -> tuple[float, float]:
    """Interpolate the polynomial through the n points (xa, ya) at x,
    starting from the nearest tabulated point.

    Translates: polint.f#POLINT (7be4133)

    Arguments: xa[1..n], ya[1..n] (n <= 20); x; y, dy (set). Returns
    (y, dy): the value at x and the last correction added, an error
    estimate.

    """

    # The nearest tabulated point, and the tableau's first column
    # (Fortran lines 21–33).
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

    # Neville's tableau, one column per pass, adding one correction to y
    # each (Fortran lines 34–59).
    for j in range(1, n):
        for i in range(1, n - j + 1):
            ho = xa[i] - x
            hp = xa[i + j] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            if abs(den) < 1.0e-20:
                text = (f" POLINT: TWO ABSCISSAE COINCIDE{xa[i]:14.6E}"
                        f"{xa[i + j]:14.6E}{x:14.6E} N{n:3d}\n"
                        " RUN STOPPED")                       # format 900
                print(text)
                iounits_mod.write(iounits_mod.ishort, text)
                raise SystemExit(1)    # STOP 1, the run-time error code
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
