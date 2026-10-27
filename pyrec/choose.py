"""Python version of choose.f (YREC 09c3b9b): the case of the shape-
preserving quadratic spline (McAllister and Roulier) for one interval."""


def choose(p1, p2, m1, m2, q1, q2, eps, ncase):
    """Choose which of the four cases of McAllister and Roulier's
    quadratic spline the interval needs, from the end slopes and the
    slope of the chord.

    Translates: choose.f#CHOOSE (09c3b9b)

    Arguments: (p1, p2), (q1, q2), the points of interpolation; m1, m2,
    the slopes there; eps, the tolerance for slopes close to the chord's
    slope or twice it; ncase, passed in (with any value).

    Returns ncase, 1 to 4.
    """
    # The slope of the chord.
    # Question (see the page): no check that q1 != p1.
    spq = (q2 - p2) / (q1 - p1)
    # Modification by Marc Pinsonneault, to avoid a division by zero in
    # CASES: a zero end slope gives case 2 at once.
    if m1 == 0.0 or m2 == 0.0:
        return 2  # labels 9, 10
    if spq == 0.0:
        if m1 * m2 >= 0.0:
            return 2  # label 10
        return 1
    # label 20
    prod1 = spq * m1
    prod2 = spq * m2
    mref = abs(spq)
    mref1 = abs(m1)
    mref2 = abs(m2)
    if abs(spq - m1) <= eps * mref or abs(spq - m2) <= eps * mref:
        # label 30: an end slope close to the chord's: case 2 or 3.
        if prod1 < 0.0 or prod2 < 0.0:
            return _opposite_signs(prod1, prod2, mref, mref1, mref2, eps)
    else:
        if prod1 < 0.0 or prod2 < 0.0:
            return _opposite_signs(prod1, prod2, mref, mref1, mref2, eps)
        prod = (mref - mref1) * (mref - mref2)
        if prod < 0.0:
            # The two tangent lines meet between p1 and q1.
            return 1
    # label 40: all three slopes have the same sign.
    if mref1 > 2.0 * mref:
        # label 50
        if mref2 > (2.0 - eps) * mref:
            return 4  # label 70
        return 3
    if mref2 > 2.0 * mref:
        # label 60
        if mref1 > (2.0 - eps) * mref:
            return 4  # label 70
        return 3
    # Both tangents cross the midline: a knot at (p1 + q1)/2.
    return 2


def _opposite_signs(prod1, prod2, mref, mref1, mref2, eps):
    """Label 80 of CHOOSE: at least one end slope has the sign opposite to
    the chord's. Returns ncase."""
    if prod1 < 0.0 and prod2 < 0.0:
        return 2  # label 130
    if prod1 < 0.0:
        # label 90
        if mref2 > (1.0 + eps) * mref:
            return 1  # label 100
        return 2
    # label 110
    if mref1 > (1.0 + eps) * mref:
        return 1  # label 120
    return 2
