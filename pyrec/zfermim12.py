"""Python version of zfermim12.f (YREC 09c3b9b): the Fermi-Dirac integral
of order -1/2 (Antia 1993)."""

import math

# DATA: Antia's coefficients (index 0 unused).
M1, K1, M2, K2 = 7, 7, 11, 11
A1 = [0.0, 1.71446374704454e7, 3.88148302324068e7, 3.16743385304962e7,
      1.14587609192151e7, 1.83696370756153e6, 1.14980998186874e5,
      1.98276889924768e3, 1.0]
B1 = [0.0, 9.67282587452899e6, 2.87386436731785e7, 3.26070130734158e7,
      1.77657027846367e7, 4.81648022267831e6, 6.13709569333207e5,
      3.13595854332114e4, 4.35061725080755e2]
A2 = [0.0, -4.46620341924942e-15, -1.58654991146236e-12,
      -4.44467627042232e-10, -6.84738791621745e-8, -6.64932238528105e-6,
      -3.69976170193942e-4, -1.12295393687006e-2, -1.60926102124442e-1,
      -8.52408612877447e-1, -7.45519953763928e-1, 2.98435207466372e0, 1.0]
B2 = [0.0, -2.23310170962369e-15, -7.94193282071464e-13,
      -2.22564376956228e-10, -3.43299431079845e-8, -3.33919612678907e-6,
      -1.86432212187088e-4, -5.69764436880529e-3, -8.34904593067194e-2,
      -4.78770844009440e-1, -4.99759250374148e-1, 1.86795964993052e0,
      4.16485970495288e-1]


def zfermim12(x):
    """F_{-1/2}(x) by Antia's two rational approximations, in exp(x) below
    x = 2 and in 1/x^2 above.

    Translates: zfermim12.f#ZFERMIM12 (09c3b9b)

    Arguments: x, the degeneracy parameter. Returns the integral.
    """
    if x < 2.0:
        xx = math.exp(x)
        rn = xx + A1[M1]
        for i in range(M1 - 1, 0, -1):
            rn = rn * xx + A1[i]
        den = B1[K1 + 1]
        for i in range(K1, 0, -1):
            den = den * xx + B1[i]
        return xx * rn / den
    xx = 1.0 / (x * x)
    rn = xx + A2[M2]
    for i in range(M2 - 1, 0, -1):
        rn = rn * xx + A2[i]
    den = B2[K2 + 1]
    for i in range(K2, 0, -1):
        den = den * xx + B2[i]
    return math.sqrt(x) * rn / den
