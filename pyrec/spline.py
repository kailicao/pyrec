"""Python version of spline.f (YREC 7be4133): the value of a
shape-preserving quadratic spline (McAllister and Roulier) at one point of
one interval."""


def _linear(xvals: float, xtabs: float, ytabs: float, xtabs1: float, ytabs1: float) -> float:
    """The fallback where a denominator would vanish (MHP 6/87): linear
    interpolation between the interval's ends. It runs only for xvals
    outside [xtabs, xtabs1] (see the page).

    """

    fs = (xvals - xtabs) / (xtabs1 - xtabs)
    return ytabs + fs * (ytabs1 - ytabs)


def spline(xvals: float, z1: float, z2: float, xtabs: float, ytabs: float, xtabs1: float, ytabs1: float, y1: float,
           y2: float, e2: float, w2: float, v2: float, ncase: int) -> float:
    """The value at xvals of the quadratic spline on [xtabs, xtabs1], made
    of quadratic pieces joined at one knot z1 (ncase 1-3) or at two knots
    y1 < z1 (ncase 4), with the parameters computed by MEVAL.

    Translates: spline.f#SPLINE (7be4133)

    Arguments: xvals; z1, z2 (the second knot and its value); xtabs, ytabs
    and xtabs1, ytabs1 (the interval's ends); y1, y2 (the first knot and
    its value, case 4); e2, w2, v2 (the middle control values); ncase.
    Returns the spline's value.

    """

    # Cases 1, 2 or 3: where xvals lies relative to the knot (Fortran
    # lines 51–87).
    if ncase != 4:
        if z1 < xvals:  # label 10
            if xtabs1 != z1:
                return (z2 * (xtabs1 - xvals) ** 2
                        + w2 * 2.0 * (xvals - z1) * (xtabs1 - xvals)
                        + ytabs1 * (xvals - z1) ** 2) / (xtabs1 - z1) ** 2
            return _linear(xvals, xtabs, ytabs, xtabs1, ytabs1)
        if z1 == xvals:  # label 20
            return z2
        if z1 != xtabs:  # label 30
            return (ytabs * (z1 - xvals) ** 2
                    + v2 * 2.0 * (xvals - xtabs) * (z1 - xvals)
                    + z2 * (xvals - xtabs) ** 2) / (z1 - xtabs) ** 2
        return _linear(xvals, xtabs, ytabs, xtabs1, ytabs1)

    # Case 4 (label 40): where xvals lies relative to the first knot
    # (Fortran lines 94–100); beyond it, label 70: where xvals lies
    # relative to the second knot (Fortran lines 118–147).
    if y1 < xvals:
        if z1 < xvals:  # label 100
            if xtabs1 != z1:
                return (z2 * (xtabs1 - xvals) ** 2
                        + w2 * 2.0 * (xvals - z1) * (xtabs1 - xvals)
                        + ytabs1 * (xvals - z1) ** 2) / (xtabs1 - z1) ** 2
            return _linear(xvals, xtabs, ytabs, xtabs1, ytabs1)
        if z1 == xvals:  # label 90
            return z2
        if z1 != y1:  # label 80
            return (y2 * (z1 - xvals) ** 2
                    + e2 * 2.0 * (xvals - y1) * (z1 - xvals)
                    + z2 * (xvals - y1) ** 2) / (z1 - y1) ** 2
        return _linear(xvals, xtabs, ytabs, xtabs1, ytabs1)

    # Case 4, up to the first knot (Fortran lines 102–113).
    if y1 == xvals:  # label 60
        return y2
    if y1 != xtabs:  # label 50
        return (ytabs * (y1 - xvals) ** 2
                + v2 * 2.0 * (xvals - xtabs) * (y1 - xvals)
                + y2 * (xvals - xtabs) ** 2) / (y1 - xtabs) ** 2
    return _linear(xvals, xtabs, ytabs, xtabs1, ytabs1)
