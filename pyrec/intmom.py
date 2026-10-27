"""Python version of intmom.f (YREC 09c3b9b): the moment of inertia per
unit mass of a rotationally distorted shell."""

import numpy as np

from pyrec import constants_mod

# DATA: the coefficients of the series.
c = np.zeros(5 + 1)
c[1:6] = [0.6, 0.72857142857143, 0.50884115884115, 0.49728653699242,
          0.42852079651326]


def intmom(a, eta2x, qr0r, r0r, him, qiwm):
    """Return the moment of inertia per unit mass of a distorted shell,
    (2/3) qr0r r0r (1 + 0.15 sum_j c_j (j eta2 + 5) a**j), and the sum's
    a d/da counterpart, qiwm.

    Translates: intmom.f#INTMOM (09c3b9b)

    Arguments: a (the distortion amplitude), eta2x (eta_2), qr0r, r0r
    (from MOMI: 1 and r0**5/r**3); him, qiwm are only set.

    Returns (him, qiwm).
    """
    # Question (see the page): the source of the series.
    h1 = constants_mod.cc23 * qr0r * r0r
    sum1 = 0.0
    sum2 = 0.0
    for j in range(1, 5 + 1):  # label 10
        sum1 = sum1 + c[j] * (j * eta2x + 5.0) * a ** j
        sum2 = sum2 + j * c[j] * (j * eta2x + 5.0) * a ** j
    him = h1 * (1.0 + 1.5e-1 * sum1)
    # Question (see the page): a dI/da, not dI/dOmega.
    qiwm = h1 * (1.5e-1 * sum2)
    return him, qiwm
