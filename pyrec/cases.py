"""Python version of cases.f (YREC 09c3b9b): the knots of a shape-preserving
quadratic spline on one interval (McAllister and Roulier)."""


def cases(p1, p2, m1, m2, q1, q2, e1, e2, v1, v2, w1, w2, z1, z2, y1, y2,
          ncase):
    """The knots of the quadratic spline on (p1, q1), for case ncase (1 to
    4) of McAllister and Roulier's method.

    Translates: cases.f#CASES (09c3b9b)

    Arguments: (p1, p2) and (q1, q2), the points of interpolation; m1, m2,
    the slopes there; ncase, which fixes the number and place of the
    knots. The outputs e1 ... y2 are passed in (with any value).

    Returns (e1, e2, v1, v2, w1, w2, z1, z2, y1, y2): the knots (v1, v2),
    (w1, w2), (z1, z2), and for case 4 also (e1, e2) and (y1, y2); those a
    case does not set are returned as passed in.
    """
    # Question (see the page): no check of the divisions by m2 - m1, m1, m2.
    if ncase not in (3, 4):
        if ncase == 2:
            # label 10: case 2.
            z1 = (p1 + q1) / 2.0
            v1 = (p1 + z1) / 2.0
            v2 = p2 + m1 * (v1 - p1)
            w1 = (z1 + q1) / 2.0
            w2 = q2 + m2 * (w1 - q1)
            z2 = (v2 + w2) / 2.0
            return e1, e2, v1, v2, w1, w2, z1, z2, y1, y2
        # Case 1 (also any NCASE other than 2, 3, 4).
        z1 = (p2 - q2 + m2 * q1 - m1 * p1) / (m2 - m1)
        ztwo = p2 + m1 * (z1 - p1)
        v1 = (p1 + z1) / 2.0
        v2 = (p2 + ztwo) / 2.0
        w1 = (z1 + q1) / 2.0
        w2 = (ztwo + q2) / 2.0
        z2 = v2 + ((w2 - v2) / (w1 - v1)) * (z1 - v1)
        return e1, e2, v1, v2, w1, w2, z1, z2, y1, y2
    # label 20: the parameters used in both cases 3 and 4.
    c1 = p1 + (q2 - p2) / m1
    d1 = q1 + (p2 - q2) / m2
    h1 = 2.0 * c1 - p1
    j1 = 2.0 * d1 - q1
    mbar1 = (q2 - p2) / (h1 - p1)
    mbar2 = (p2 - q2) / (j1 - q1)
    if ncase == 3:
        k1 = (p2 - q2 + q1 * mbar2 - p1 * mbar1) / (mbar2 - mbar1)
        if abs(m1) > abs(m2):
            z1 = (k1 + p1) / 2.0  # label 30
        else:
            z1 = (k1 + q1) / 2.0
        # label 40
        v1 = (p1 + z1) / 2.0
        v2 = p2 + m1 * (v1 - p1)
        w1 = (q1 + z1) / 2.0
        w2 = q2 + m2 * (w1 - q1)
        z2 = v2 + ((w2 - v2) / (w1 - v1)) * (z1 - v1)
        return e1, e2, v1, v2, w1, w2, z1, z2, y1, y2
    # label 50: case 4.
    y1 = (p1 + c1) / 2.0
    v1 = (p1 + y1) / 2.0
    v2 = m1 * (v1 - p1) + p2
    z1 = (d1 + q1) / 2.0
    w1 = (q1 + z1) / 2.0
    w2 = m2 * (w1 - q1) + q2
    mbar3 = (w2 - v2) / (w1 - v1)
    y2 = mbar3 * (y1 - v1) + v2
    z2 = mbar3 * (z1 - v1) + v2
    e1 = (y1 + z1) / 2.0
    e2 = mbar3 * (e1 - v1) + v2
    return e1, e2, v1, v2, w1, w2, z1, z2, y1, y2
