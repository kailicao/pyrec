"""Python version of qgauss.f (YREC 09c3b9b): ten-point Gauss-Legendre
quadrature over the polar angle on one equipotential shell."""

import numpy as np

from pyrec.func import func

# DATA: the five positive nodes and their weights of the 10-point rule.
x = np.array([0.0, .14887433898163, .43339539412925, .67940956829902,
              .86506336668899, .97390652851717])
w = np.array([0.0, .29552422471475, .26926671931, .21908636251598,
              .14945134915058, .06667134430869])


def qgauss(g0g, ginvg, sphig, b, r0, hs, aint, q, w2, a, i):
    """Integrate g dS, dS/g and dS over theta in [0, b] on shell i, with
    FUNC giving the effective gravity g and the surface element S.

    Translates: qgauss.f#QGAUSS (09c3b9b)

    Arguments: g0g, ginvg, sphig (set); b, the upper limit; r0, hs, aint,
    q, w2, a, i, passed through to FUNC. Returns (g0g, ginvg, sphig).
    """
    xm = 0.5 * b
    xr = xm
    g0g = 0.0
    ginvg = 0.0
    sphig = 0.0
    g = s = g2 = s2 = 0.0
    for j in range(1, 6):
        dx = xr * x[j]
        g, s = func(xm + dx, g, s, r0, hs, aint, q, w2, a, i)
        g2, s2 = func(xm - dx, g2, s2, r0, hs, aint, q, w2, a, i)
        g0g = g0g + w[j] * (g * s + g2 * s2)
        ginvg = ginvg + w[j] * (s / g + s2 / g2)
        sphig = sphig + w[j] * (s + s2)
    g0g = g0g * xr
    ginvg = ginvg * xr
    sphig = sphig * xr
    return g0g, ginvg, sphig
