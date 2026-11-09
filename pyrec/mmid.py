"""Python version of mmid.f (YREC 7be4133): the modified midpoint step of
Numerical Recipes, the inner step of BSSTEP."""

from collections.abc import Callable

import numpy as np

# The Fortran's local arrays, kept by SAVE (Fortran line 19).
ym = np.zeros(4)
yn = np.zeros(4)


def _call_deriv(deriv, x0, y, dydx, b, fpl, ftl, gl, latmo, lderiv, locond,
                lsave, rl, teffl, x, z, kount, ksaha):
    # DERIV is QENV, which returns (rl, kenv, ksaha), or QATM, which returns
    # (katm, ksaha); the Fortran passes RL, KOUNT and KSAHA by reference.
    result = deriv(x0, y, dydx, b, fpl, ftl, gl, latmo, lderiv, locond,
                   lsave, rl, teffl, x, z, kount, ksaha)
    if len(result) == 3:
        rl, kount, ksaha = result
    else:
        kount, ksaha = result
    return rl, kount, ksaha


def mmid(y: np.ndarray, dydx: np.ndarray, nvar: int, xs: float, htot: float, nstep: int, yout: np.ndarray,
         deriv: Callable, b: float, fpl: float, ftl: float, gl: float, latmo: bool, lderiv: bool, locond: bool,
         lsave: bool, rl: float, teffl: float, x: float, z: float, kount: int, ksaha: int) -> tuple[float, int, int]:
    """Advance y from xs to xs + htot in nstep substeps by the modified
    midpoint method; the result goes to yout.

    Translates: mmid.f#MMID (7be4133)

    Arguments: y and dydx at xs (three variables, index 0 unused; nvar of
    them used); the step htot and the number of substeps nstep; yout, set
    in place; deriv, the derivative routine (qenv or qatm); the arguments
    from b to ksaha are only passed on to deriv.

    Returns (rl, kount, ksaha), as deriv leaves them.

    """

    # The substep (Fortran line 28).
    h = htot / float(nstep)

    # First step (Fortran lines 29–34).
    for i in range(1, nvar + 1):                         # loop 10
        ym[i] = y[i]
        yn[i] = y[i] + dydx[i] * h
    x0 = xs + h

    # yout holds the derivatives in between (Fortran lines 35–38).
    rl, kount, ksaha = _call_deriv(deriv, x0, yn, yout, b, fpl, ftl, gl,
                                   latmo, lderiv, locond, lsave, rl, teffl,
                                   x, z, kount, ksaha)
    h2 = 2.0 * h

    # General step (Fortran lines 39–49).
    for n in range(2, nstep + 1):                        # loop 30
        for i in range(1, nvar + 1):
            swap = ym[i] + h2 * yout[i]
            ym[i] = yn[i]
            yn[i] = swap
        x0 = x0 + h
        rl, kount, ksaha = _call_deriv(deriv, x0, yn, yout, b, fpl, ftl, gl,
                                       latmo, lderiv, locond, lsave, rl,
                                       teffl, x, z, kount, ksaha)

    # Last step (Fortran lines 50–54).
    for i in range(1, nvar + 1):                         # loop 40
        yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i])
    return rl, kount, ksaha
