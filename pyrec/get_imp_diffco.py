"""Python version of get_imp_diffco.f (YREC 7be4133): the implicit
diffusion step's tridiagonal matrix, after correcting the coefficients."""

import numpy as np


def get_imp_diffco(alpha: np.ndarray, ecod2_h: np.ndarray, ex_h: np.ndarray, eqcod2x_h: np.ndarray, a: np.ndarray,
                   b: np.ndarray, c: np.ndarray, npt: int) -> None:
    """Correct the diffusion coefficients at the half points for the
    last change in X, then build the backward-Euler tridiagonal matrix
    (a, b, c) with zero-flux ends.

    Translates: get_imp_diffco.f#GET_IMP_DIFFCO (7be4133)

    Arguments: alpha (the step's coefficients per point), ex_h (the last
    change in X), eqcod2x_h (dD/dX), npt (number of points).
    Arrays changed in place: ecod2_h (D at the half points), a, b, c.

    Returns nothing.

    """

    # Correct the coefficients for the last change in X (Fortran lines
    # 19–21). Suspected bug (see GRSETT's page): GRSETT's ex_h is one
    # midpoint off.
    for i in range(1, npt):  # label 10
        ecod2_h[i] = ecod2_h[i] + eqcod2x_h[i] * ex_h[i]

    # The tridiagonal matrix (Fortran lines 23–34).
    a[1] = 0.0
    b[1] = 1.0 + alpha[1] * ecod2_h[1]
    c[1] = -alpha[1] * ecod2_h[1]
    for i in range(2, npt):  # label 20
        a[i] = -alpha[i] * ecod2_h[i - 1]
        b[i] = 1.0 + alpha[i] * (ecod2_h[i - 1] + ecod2_h[i])
        c[i] = -alpha[i] * ecod2_h[i]
    a[npt] = -alpha[npt] * ecod2_h[npt - 1]
    b[npt] = 1.0 + alpha[npt] * ecod2_h[npt - 1]
    c[npt] = 0.0
