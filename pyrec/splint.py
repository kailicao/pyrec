"""Python version of splint.f (YREC 7be4133): cubic-spline interpolation at
one point (Numerical Recipes' SPLINT), with arrays sized N."""

from pyrec import iounits_mod
import numpy as np


def splint(xa: np.ndarray, ya: np.ndarray, n: int, y2a: np.ndarray, x: float, y: float, klo: int,
           khi: int) -> tuple[float, int, int]:
    """The cubic-spline value y at x, from the table (xa, ya) and its second
    derivatives y2a (from SPLINJ or CSPLINE). Stops the run if two
    abscissas coincide; outside the table it extrapolates the end interval.

    Translates: splint.f#SPLINT (7be4133)

    Arguments: xa[1..n] (increasing), ya[1..n]; n; y2a[1..n]; x; y, klo and
    khi, set only (klo and khi are reset on entry). Returns (y, klo, khi),
    with klo, khi the bracketing indices.

    """

    # Bisection for the bracketing interval (Fortran lines 30–40).
    klo = 1
    khi = n
    while khi - klo > 1:  # label 1
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k

    # Stop if two abscissas coincide (Fortran lines 41–50).
    h = xa[khi] - xa[klo]
    if h == 0.0:
        text = (f" ERROR IN SPLINT ROUTINE: TWO ABSCISSAE COINCIDE AT"
                f"{xa[klo]:14.6E} X{x:14.6E}\n RUN STOPPED")   # format 900
        print(text)
        iounits_mod.write(iounits_mod.ishort, text)
        raise SystemExit(1)    # STOP 1, the run-time error code

    # The cubic spline's value (Fortran lines 51–54).
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi]
         + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * (h ** 2) / 6.0)
    return y, klo, khi
