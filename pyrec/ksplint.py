"""Python version of ksplint.f (YREC development, 09c3b9b)."""

NM = 4


def ksplint(xa, ya, y2a, x, y):
    """Evaluate at ``x`` the 4-point cubic spline from KSPLINE.

    *Numerical Recipes*' SPLINT with ``NM`` = 4: bisection for the
    interval, then the cubic-spline formula. Outside the nodes it
    extrapolates with the end interval's cubic.

    Translates: ksplint.f#KSPLINT (09c3b9b)

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
    klo = 1
    khi = NM
    while khi - klo > 1:  # label 1
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    if h == 0.0:
        print(" Ksplint failure")
        raise SystemExit
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    y = (a * ya[klo] + b * ya[khi]
         + ((a ** 3 - a) * y2a[klo] + (b ** 3 - b) * y2a[khi]) * (h ** 2) / 6.0)
    return y
