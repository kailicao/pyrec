"""Python version of osplin.f (YREC 09c3b9b): interpolation by a
shape-preserving osculatory quadratic spline (ACM Algorithm 574)."""

import numpy as np

from pyrec import dimensions_mod
from pyrec.meval import meval
from pyrec.slopes import slopes

# The local slope array M (SAVE).
m = np.zeros(dimensions_mod.json + 1)


def osplin(xval, yval, xtab, ytab, n, k):
    """Interpolate the table (xtab, ytab) of n points to the k abscissas
    xval, with a spline that preserves the data's monotonicity and
    convexity: SLOPES gives the slopes at the data points, MEVAL evaluates.

    Translates: osplin.f#OSPLIN (09c3b9b)

    Arguments: xval[1..k], the abscissas; yval, set in place to the values;
    xtab[1..n], ytab[1..n], the table; n, k. Returns nothing.
    """
    slopes(xtab, ytab, m, n)
    eps = 1.0e-04
    # Question (see the page): MEVAL's error code is ignored.
    err = meval(xval, yval, xtab, ytab, m, n, k, eps, 0)
