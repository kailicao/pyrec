"""Python version of osplin.f (YREC 7be4133): interpolation by a
shape-preserving osculatory quadratic spline (ACM Algorithm 574)."""

import numpy as np

from pyrec import dimensions_mod
from pyrec.meval import meval
from pyrec.slopes import slopes

# The local slope array M (SAVE).
m = np.zeros(dimensions_mod.json + 1)


def osplin(xval: np.ndarray, yval: np.ndarray, xtab: np.ndarray, ytab: np.ndarray, n: int, k: int) -> None:
    """Interpolate the table (xtab, ytab) of n points to the k abscissas
    xval, with a spline that preserves the data's monotonicity and
    convexity: SLOPES gives the slopes at the data points, MEVAL evaluates.

    Translates: osplin.f#OSPLIN (7be4133)

    Arguments: xval[1..k], the abscissas; yval, set in place to the values;
    xtab[1..n], ytab[1..n], the table; n, k. Returns nothing.

    """

    # The slopes at the data points (Fortran line 42).
    slopes(xtab, ytab, m, n)

    # The tolerance, and the evaluation (Fortran lines 45–47).
    # Suspected bug (see the page): MEVAL's error code is ignored.
    eps = 1.0e-04
    err = meval(xval, yval, xtab, ytab, m, n, k, eps, 0)  # noqa: F841  (unused in the Fortran too)
