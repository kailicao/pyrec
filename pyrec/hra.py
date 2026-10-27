"""Python version of hra.f (YREC 09c3b9b): the Harvard Reference
Atmosphere's T(tau) relation."""

import math


def hra(tau):
    """Return log10 T at optical depth tau from a ninth-degree polynomial
    in log10 tau (a fit to the Harvard Reference Atmosphere).

    Translates: hra.f#HRA (09c3b9b)

    Arguments: tau.

    Returns log10 T.
    """
    # Question (see the page): the range of tau the fit is valid for.
    tl = math.log10(tau)
    x = 3.81152046471
    x = x + 0.146133736471 * tl
    tl2 = tl * tl
    x = x + 0.0267719174279 * tl2
    tl3 = tl2 * tl
    x = x - 0.029280655317 * tl3
    tl4 = tl3 * tl
    x = x - 0.0123814456666 * tl4
    tl5 = tl4 * tl
    x = x + 0.00285734990893 * tl5
    tl6 = tl5 * tl
    x = x + 0.0024575213331 * tl6
    tl7 = tl6 * tl
    x = x + 0.000521560431455 * tl7
    tl8 = tl7 * tl
    x = x + 4.71770176883e-5 * tl8
    tl9 = tl8 * tl
    x = x + 1.58685112637e-6 * tl9
    return x
