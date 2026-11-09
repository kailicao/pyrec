"""Python version of fpft.f (YREC 7be4133): the rotational form factors
f_P, f_T and the mean effective gravity on each equipotential."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import quadd_mod
from pyrec.polint import polint
from pyrec.qgauss import qgauss
from pyrec.shape import shape
from pyrec.trapzd import trapzd


def fpft(hd: np.ndarray, hr: np.ndarray, hs: np.ndarray, m: int, omega: np.ndarray, eta2: np.ndarray, fp: np.ndarray,
         ft: np.ndarray, hg: np.ndarray, r0: np.ndarray) -> None:
    """For every shell of a rotating model, outward: the distortion
    amplitude, the distortion integral (trapezoid rule with one Romberg
    step), the angular averages by Gauss-Legendre quadrature, and from
    them the form factors f_P, f_T and <g>; also the spherical,
    rotational and distortion potentials in quadd_mod.

    Translates: fpft.f#FPFT (7be4133)

    Arguments: hd, hr, hs (log rho, log r, log m, cgs), m (number of
    shells), omega (angular velocity per shell).
    Arrays changed in place: eta2, r0 (from SHAPE), fp, ft, hg, and
    quadd_mod.phisp, phirot, phidis.

    Returns nothing.

    """

    c = constants_mod

    # Work arrays (Fortran lines 31–33).
    h = np.zeros(20 + 1)
    xa = np.zeros(20 + 1)
    ya = np.zeros(20 + 1)
    aint0 = np.zeros(10 + 1)

    # The run of r0 and eta2 for the model, and the constants (Fortran lines 37–44).
    shape(hd, hr, hs, 1, m, omega, eta2, r0)
    b = 0.125 * c.c4pi
    cg = math.exp(c.cln * c.cgl)
    eps = 1.0e-6
    jmax = 2
    k = 2
    km = k - 1

    # The values at the previous shell (at the centre for the first) (Fortran lines 45–55).
    mp = 0.0
    eta22p = 2.0
    w2p = omega[1] ** 2
    rhop = math.exp(c.cln * hd[1])
    r0p = 0.0
    aintp = 0.0
    aintt = 0.0  # noqa: F841  (unused in the Fortran too)
    qp = 0.0
    q = 0.0
    aint1 = 0.0

    for i in range(1, m + 1):  # label 100

        # The shell's distortion amplitude and the values at it (Fortran lines 58–65).
        rho = math.exp(c.cln * hd[i])
        gm = cg * math.exp(c.cln * hs[i])
        fact = 5.0 * c.cc13 * omega[i] ** 2 / (gm * (2.0 + eta2[i]))
        r03 = r0[i] ** 3
        a = fact * r03
        sm = math.exp(c.cln * hs[i])
        eta22 = eta2[i] + 2.0
        w2 = omega[i] ** 2

        # The distortion integral over this step: trapezoid rule, refined
        # and extrapolated to zero step (Romberg) (Fortran lines 66–85).
        h[1] = 1.0
        for j in range(1, jmax + 1):  # label 40
            aint0[j], q = trapzd(r0p, r0[i], aint0[j], j, rho, rhop, sm, mp,
                                 w2, w2p, eta22, eta22p, q, qp)
            if j >= k:
                n1 = j - km
                for j1 in range(1, k + 1):  # label 35
                    xa[j1] = h[n1]
                    ya[j1] = aint0[n1]
                    n1 = n1 + 1
                aint1, dint = polint(xa, ya, k, 0.0, 0.0, 0.0)
                # Question (see the page): with jmax = 2 the result is
                # used whether or not the tolerance is met.
                if abs(dint) < eps * abs(aint1):
                    break  # GOTO 50
            aint0[j + 1] = aint0[j]
            h[j + 1] = 0.25 * h[j]
        aint = aintp + aint1  # label 50

        # <g> and <g^-1> by Gaussian quadrature, and the form factors (Fortran lines 86–92).
        g0, ginv0, sphi = qgauss(0.0, 0.0, 0.0, b, r0, hs, aint, q, w2, a, i)
        hg[i] = g0 / sphi
        ginv = math.log10(ginv0)
        f0 = c.c4pil + 4.0 * hr[i] - ginv
        fp[i] = math.exp(c.cln * (f0 - c.cgl - hs[i]))
        ft[i] = math.exp(c.cln * (f0 + c.c4pil)) / g0

        # Output data: the potentials (Fortran lines 93–98).
        rphi = math.exp(c.cln * hr[i])
        rphi3 = rphi ** 3
        quadd_mod.phisp[i] = gm / rphi
        quadd_mod.phirot[i] = w2 * rphi ** 2
        quadd_mod.phidis[i] = c.c4pi * c.cc13 * aint / rphi3

        # The values for the next shell (Fortran lines 99–105).
        aintp = aint
        qp = q
        mp = sm
        eta22p = eta22
        w2p = w2
        rhop = rho
        r0p = r0[i]
