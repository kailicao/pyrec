"""Python version of ksplint.f (YREC development, 7be4133)."""

from pyrec import iounits_mod
import numpy as np

NM = 4


def ksplint(xa: np.ndarray, ya: np.ndarray, y2a: np.ndarray, x: float, y: float) -> float:
    """Evaluate at ``x`` the 4-point cubic spline from KSPLINE.

    *Numerical Recipes*' SPLINT with ``NM`` = 4: bisection for the
    interval, then the cubic-spline formula. Outside the nodes it
    extrapolates with the end interval's cubic.

    Translates: ksplint.f#KSPLINT (7be4133)

    Parameters
    ----------
    xa, ya, y2a : ndarray, shape (NM + 1,)
        Nodes, values and second derivatives, indices 1..4.
    x : float
        The point of evaluation.
    y : float
        Output only (passed in with any value).

    Returns
    -------
    y : float
        The spline's value at ``x``.

    """

    # Bisection for the bracketing interval (Fortran lines 15–26).
    klo = 1
    khi = NM
    while khi - klo > 1:  # label 1
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k

    # Stop if two nodes coincide (Fortran lines 28–36).
    h = xa[khi] - xa[klo]
    if h == 0.0:
        text = (f" KSPLINT: TWO ABSCISSAE COINCIDE AT{xa[klo]:14.6E}"
                f" X{x:14.6E}\n RUN STOPPED")                  # format 900
        print(text)
        iounits_mod.write(iounits_mod.ishort, text)
        raise SystemExit(1)    # STOP 1, the run-time error code

    # The cubic spline's value (Fortran lines 38–41).
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi]
         + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * (h ** 2) / 6.0)
    return y
