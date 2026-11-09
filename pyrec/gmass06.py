"""Python version of gmass06.f (YREC 7be4133): the mean mass per nucleus
and number fractions for the 2006 OPAL equation of state."""

import numpy as np

# DATA: the ground-state energies of the ions and the nuclear charges, in
# the order Ne, O, N, C, He, H.
eion = np.zeros(7 + 1)
eion[1:7] = [-3394.873554, -1974.86545, -1433.92718, -993.326315,
             -76.1959403, -15.29409]
anum = np.zeros(6 + 1)
anum[1:7] = [10., 8., 7., 6., 2., 1.]
amas = np.zeros(7 + 1)


def gmass06(x: float, z: float, amoles: float, eground: float, fracz: float, frac: np.ndarray) -> tuple[float, float,
            float, float]:
    """Split Z in fixed number proportions among C, N, O and Ne, and
    compute the number fractions frac (Ne, O, N, C, He, H), the
    particles per nucleus amoles, the ground-state energy eground and the
    mean mass per nucleus (without electrons).

    Translates: gmass06.f#GMASS06 (7be4133)

    Arguments: x, z; amoles, eground, fracz are only set.
    Array changed in place: frac (indices 1-6).

    Returns (gmass06, amoles, eground, fracz).

    """

    # The metals' number fractions and the atomic masses (Fortran lines
    # 22–32).
    xc = 0.247137766
    xn = .0620782
    xo = .52837118
    xne = .1624188
    amas[7] = 1.0079
    amas[6] = 4.0026
    amas[5] = 12.011
    amas[4] = 14.0067
    amas[3] = 15.9994
    amas[2] = 20.179
    amas[1] = 0.00054858

    # The number fractions (Fortran lines 34–47).
    fracz = z / (xc * amas[5] + xn * amas[4] + xo * amas[3] + xne * amas[2])
    xc2 = fracz * xc
    xn2 = fracz * xn
    xo2 = fracz * xo
    xne2 = fracz * xne
    xh = x / amas[7]
    xhe = (1.0 - x - z) / amas[6]
    xtot = xh + xhe + xc2 + xn2 + xo2 + xne2
    frac[6] = xh / xtot
    frac[5] = xhe / xtot
    frac[4] = xc2 / xtot
    frac[3] = xn2 / xtot
    frac[2] = xo2 / xtot
    frac[1] = xne2 / xtot

    # Ground-state energy and particles per nucleus (Fortran lines 48–54).
    eground = 0.0
    amoles = 0.0
    for i in range(1, 6 + 1):
        eground = eground + eion[i] * frac[i]
        amoles = amoles + (1.0 + anum[i]) * frac[i]
    anume = amoles - 1.0  # noqa: F841  (unused in the Fortran too)

    # The mean mass per nucleus (Fortran lines 55–58).
    gmass = 0.0
    for i in range(2, 7 + 1):
        gmass = gmass + amas[i] * frac[i - 1]
    return gmass, amoles, eground, fracz
