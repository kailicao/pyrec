"""Python version of meval.f (YREC development, 09c3b9b)."""

from pyrec import dimensions_mod  # noqa: F401  (JSON, the declared size)
from pyrec.cases import cases
from pyrec.choose import choose
from pyrec.search import search
from pyrec.spline import spline


def meval(xval, yval, xtab, ytab, mtab, num, nume, eps, err):
    """Evaluate McAllister and Roulier's shape-preserving quadratic spline.

    Walks through the evaluation points ``xval`` and the data intervals
    together: points equal to a data abscissa get the data value; for each
    interval holding points, CHOOSE and CASES set the spline's parameters
    once and SPLINE evaluates it; points beyond the ends are extrapolated
    with the end intervals' splines. Called only by OSPLIN.

    Translates: meval.f#MEVAL (09c3b9b)

    Parameters
    ----------
    xval : ndarray
        The points of evaluation, ``xval[1..nume]``, nondecreasing.
    yval : ndarray
        Set in place: the spline's values at ``xval``.
    xtab, ytab, mtab : ndarray
        The data abscissae (increasing), ordinates and slopes,
        indices 1..num.
    num, nume : int
        The number of data points and of points of evaluation.
    eps : float
        Relative tolerance passed to CHOOSE.
    err : int
        Output only (passed in with any value).

    Returns
    -------
    err : int
        0 normal; 1 extrapolation occurred; 2 ``xval`` not sorted (then
        nothing is computed).
    """
    # SAVE in the Fortran: every local is set before it is read within a
    # call, so all are kept local here. The spline's parameters start at
    # zero; CASES returns those a case does not set as passed in, and
    # SPLINE reads only those the case set.
    e1 = e2 = v1 = v2 = w1 = w2 = z1 = z2 = y1 = y2 = 0.0
    start = 1
    end = nume
    err = 0
    if nume != 1:
        # Determine if xval is nondecreasing (loop 10).
        for i in range(1, nume):
            if xval[i + 1] < xval[i]:
                # Question (see the page): OSPLIN does not check err, and
                # yval is left as it was.
                err = 2
                return err  # label 230

    # Question (see the page): the Fortran's comment here says points
    # beyond the data get the end values; the code extrapolates instead.
    # Points below the first data abscissa (label 20, loop 30).
    for i in range(1, nume + 1):
        if xval[i] >= xtab[1]:
            break
        start = i + 1
    # Points above the last data abscissa (label 40, loop 50).
    for i in range(1, nume + 1):
        ind = nume + 1 - i
        if xval[ind] <= xtab[num]:
            break
        end = ind - 1

    # Extrapolate below the first data point (label 60).
    if start != 1:
        err = 1
        ncase = choose(xtab[1], ytab[1], mtab[1], mtab[2], xtab[2],
                       ytab[2], eps, 0)
        e1, e2, v1, v2, w1, w2, z1, z2, y1, y2 = cases(
            xtab[1], ytab[1], mtab[1], mtab[2], xtab[2], ytab[2],
            e1, e2, v1, v2, w1, w2, z1, z2, y1, y2, ncase)
        for i in range(1, start):  # loop 70
            yval[i] = spline(xval[i], z1, z2, xtab[1], ytab[1], xtab[2],
                             ytab[2], y1, y2, e2, w2, v2, ncase)
        if nume == 1:
            return err

    # label 80
    goto_200 = (nume == 1) and (end != nume)
    if not goto_200:
        # The interval of the first in-range point. Possible bug (see
        # update-findings.md): if every point lies below xtab[1], start is
        # nume + 1 and the Fortran reads XVAL(NUME+1) here.
        lcn, fnd = search(xtab, num, xval[start], 0, 0)
        lcn1 = lcn + 1
        if fnd != 0:
            # Points equal to a data point take its value (label 90).
            while True:
                yval[start] = ytab[lcn]
                start1 = start
                start = start + 1
                if start > nume:
                    return err
                if xval[start1] != xval[start]:
                    break
            # Advance until a point lies strictly inside an interval
            # (labels 100, 110, 120).
            while xval[start] >= xtab[lcn1]:  # label 100
                if xval[start] == xtab[lcn1]:
                    while True:  # label 110
                        yval[start] = ytab[lcn1]
                        start1 = start
                        start = start + 1
                        if start > nume:
                            return err
                        if xval[start] != xval[start1]:
                            break
                lcn = lcn1  # label 120
                lcn1 = lcn1 + 1

        # The points within range of the data (label 130). With
        # extrapolation below, the parameters of interval 1 are known.
        if not (lcn == 1 and err == 1):
            ncase = choose(xtab[lcn], ytab[lcn], mtab[lcn], mtab[lcn1],
                           xtab[lcn1], ytab[lcn1], eps, 0)
            e1, e2, v1, v2, w1, w2, z1, z2, y1, y2 = cases(
                xtab[lcn], ytab[lcn], mtab[lcn], mtab[lcn1],
                xtab[lcn1], ytab[lcn1],
                e1, e2, v1, v2, w1, w2, z1, z2, y1, y2, ncase)
        for i in range(start, end + 1):  # label 140, loop 190
            if xval[i] < xtab[lcn1]:
                # label 150: the parameters of this interval are known.
                yval[i] = spline(xval[i], z1, z2, xtab[lcn], ytab[lcn],
                                 xtab[lcn1], ytab[lcn1], y1, y2, e2, w2,
                                 v2, ncase)
            elif xval[i] == xtab[lcn1]:
                yval[i] = ytab[lcn1]  # label 160
            else:
                while True:
                    lcn = lcn1  # label 170
                    lcn1 = lcn + 1
                    if xval[i] < xtab[lcn1]:
                        # label 180, then 150
                        ncase = choose(xtab[lcn], ytab[lcn], mtab[lcn],
                                       mtab[lcn1], xtab[lcn1], ytab[lcn1],
                                       eps, 0)
                        e1, e2, v1, v2, w1, w2, z1, z2, y1, y2 = cases(
                            xtab[lcn], ytab[lcn], mtab[lcn], mtab[lcn1],
                            xtab[lcn1], ytab[lcn1],
                            e1, e2, v1, v2, w1, w2, z1, z2, y1,
                            y2, ncase)
                        yval[i] = spline(xval[i], z1, z2, xtab[lcn],
                                         ytab[lcn], xtab[lcn1], ytab[lcn1],
                                         y1, y2, e2, w2, v2, ncase)
                        break
                    if xval[i] == xtab[lcn1]:
                        yval[i] = ytab[lcn1]  # label 160
                        break

        # The points beyond the last data point.
        if end == nume:
            return err
        goto_200 = not (lcn1 == num and xval[end] != xtab[num])

    if goto_200:
        # label 200: extrapolation above the last data point.
        err = 1
        num1 = max(num - 1, 1)
        ncase = choose(xtab[num1], ytab[num1], mtab[num1], mtab[num],
                       xtab[num], ytab[num], eps, 0)
        e1, e2, v1, v2, w1, w2, z1, z2, y1, y2 = cases(
            xtab[num1], ytab[num1], mtab[num1], mtab[num], xtab[num],
            ytab[num], e1, e2, v1, v2, w1, w2, z1, z2, y1, y2,
            ncase)
    # label 210
    end1 = end + 1
    num1 = max(num - 1, 1)
    for i in range(end1, nume + 1):  # loop 220
        yval[i] = spline(xval[i], z1, z2, xtab[num1], ytab[num1],
                         xtab[num], ytab[num], y1, y2, e2, w2, v2, ncase)
    return err
