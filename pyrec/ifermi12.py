"""Python version of ifermi12.f (YREC 09c3b9b): the inverse Fermi-Dirac
integral of order 1/2 (Antia 1993)."""

import math

# DATA: Antia's coefficients (index 0 unused; the rest of the Fortran's
# 12-element arrays is never set or used).
AN = 0.5
M1, K1, M2, K2 = 4, 3, 6, 5
A1 = [0.0, 1.999266880833e4, 5.702479099336e3, 6.610132843877e2,
      3.818838129486e1, 1.0]
B1 = [0.0, 1.771804140488e4, -2.014785161019e3, 9.130355392717e1,
      -1.670718177489e0]
A2 = [0.0, -1.277060388085e-2, 7.187946804945e-2, -4.262314235106e-1,
      4.997559426872e-1, -1.285579118012e0, -3.930805454272e-1, 1.0]
B2 = [0.0, -9.745794806288e-3, 5.485432756838e-2, -3.299466243260e-1,
      4.077841975923e-1, -1.145531476975e0, -6.067091689181e-2]


def ifermi12(f):
    """The degeneracy parameter eta for which F_{1/2}(eta) = f, by Antia's
    two rational approximations, split at f = 4.

    Translates: ifermi12.f#IFERMI12 (09c3b9b)

    Arguments: f. Returns the inverse.
    """
    # Question (see the page): the normalization of F_{1/2}.
    if f < 4.0:
        rn = f + A1[M1]
        for i in range(M1 - 1, 0, -1):
            rn = rn * f + A1[i]
        den = B1[K1 + 1]
        for i in range(K1, 0, -1):
            den = den * f + B1[i]
        return math.log(f * rn / den)
    ff = 1.0 / f**(1.0 / (1.0 + AN))
    rn = ff + A2[M2]
    for i in range(M2 - 1, 0, -1):
        rn = rn * ff + A2[i]
    den = B2[K2 + 1]
    for i in range(K2, 0, -1):
        den = den * ff + B2[i]
    return rn / (den * ff)
