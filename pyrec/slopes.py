"""Python version of slopes.f (YREC 09c3b9b): the derivatives at the nodes
for the shape-preserving quadratic splines of McAllister and Roulier, and
KC's SAFEDIVIDE."""


def slopes(xtab, ytab, mtab, num):
    """The derivative mtab[i] at each node, chosen so that the osculatory
    quadratic spline has one extra knot per interval and keeps convexity
    and monotonicity where the data allow.

    Translates: slopes.f#SLOPES (09c3b9b)

    Arguments: xtab[1..num], ytab[1..num] (num >= 3); mtab[1..num], set in
    place; num. None is returned.
    """
    num1 = num - 1
    im1 = 1
    i = 2
    i1 = 3
    # Question (see the page): unprotected divisions remain.
    ydif1 = ytab[2] - ytab[1]
    ydif2 = ytab[3] - ytab[2]
    m1 = ydif1 / (xtab[2] - xtab[1])
    m1s = m1
    m2 = ydif2 / (xtab[3] - xtab[2])
    m2s = m2
    # The interior nodes.
    while True:
        if m1 == 0.0 or m2 == 0.0 or (m1 * m2) <= 0.0:  # label 10
            mtab[i] = 0.0  # label 20
        elif abs(m1) > abs(m2):
            xbar = (ydif2 / m1) + xtab[i]  # label 30
            xhat = (xbar + xtab[i1]) / 2.0
            mtab[i] = ydif2 / (xhat - xtab[i])
        else:
            xbar = (-ydif1 / m2) + xtab[i]  # label 40
            xhat = (xtab[im1] + xbar) / 2.0
            mtab[i] = ydif1 / (xtab[i] - xhat)
        im1 = i  # label 50
        i = i1
        i1 = i1 + 1
        if i > num1:
            break
        ydif1 = ytab[i] - ytab[im1]
        ydif2 = ytab[i1] - ytab[i]
        m1 = ydif1 / (xtab[i] - xtab[im1])
        m2 = safedivide(ydif2, (xtab[i1] - xtab[i]), m2)

    # The last node (label 60).
    if (m1 * m2) < 0.0:
        mtab[num] = 2.0 * m2  # label 80
    else:
        xmid = (xtab[num1] + xtab[num]) / 2.0
        yxmid = mtab[num1] * (xmid - xtab[num1]) + ytab[num1]
        mtab[num] = safedivide((ytab[num] - yxmid), (xtab[num] - xmid),
                               mtab[num])
        if (mtab[num] * m2) < 0.0:
            mtab[num] = 0.0  # label 70

    # The first node (label 90).
    if (m1s * m2s) < 0.0:
        mtab[1] = 2.0 * m1s  # label 110
        return
    xmid = (xtab[1] + xtab[2]) / 2.0
    yxmid = mtab[2] * (xmid - xtab[2]) + ytab[2]
    mtab[1] = (yxmid - ytab[1]) / (xmid - xtab[1])
    if (mtab[1] * m1s) < 0.0:
        mtab[1] = 0.0  # label 100


def safedivide(numerator, denominator, result):
    """numerator / denominator, or 0 when either is zero (KC 2025-05-31).

    Translates: slopes.f#SAFEDIVIDE (09c3b9b)

    Arguments: numerator, denominator; result, set only. Returns result.
    """
    result = 0.0
    if numerator != 0.0 and denominator != 0.0:
        result = numerator / denominator
    return result
