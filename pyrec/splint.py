"""Python version of splint.f (YREC 09c3b9b): cubic-spline interpolation at
one point (Numerical Recipes' SPLINT), with arrays sized N."""

from pyrec import iounits_mod


def splint(xa, ya, n, y2a, x, y, klo, khi):
    """The cubic-spline value y at x, from the table (xa, ya) and its second
    derivatives y2a (from SPLINJ or CSPLINE). Stops the run if two
    abscissas coincide; outside the table it extrapolates the end interval.

    Translates: splint.f#SPLINT (09c3b9b)

    Arguments: xa[1..n] (increasing), ya[1..n]; n; y2a[1..n]; x; y, klo and
    khi, set only (klo and khi are reset on entry). Returns (y, klo, khi),
    with klo, khi the bracketing indices.
    """
    klo = 1
    khi = n
    while khi - klo > 1:  # label 1
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if h == 0.0:
        iounits_mod.write(iounits_mod.ishort, " ERROR IN SPLINT ROUTINE.")
        raise SystemExit()
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi]
         + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * (h ** 2) / 6.0)
    return y, klo, khi
