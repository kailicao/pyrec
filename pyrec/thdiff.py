"""Python version of thdiff.f (YREC 09c3b9b): Burgers' equations for the
diffusion of the species in a fully ionised plasma (Thoul et al. 1994)."""

import math

import numpy as np

from pyrec.lubksb import lubksb
from pyrec.ludcmp import ludcmp

MMAX = 20
NMAX = 42


def thdiff(m, a, z, x, cl, ap, at, ax):
    """Solve Burgers' equations (Thoul, Bahcall & Loeb 1994) for the
    coefficients of the pressure, temperature and concentration gradients
    in each species' diffusion velocity.

    Translates: thdiff.f#THDIFF (09c3b9b)

    Arguments: m, the number of species (1 hydrogen, 2 helium, 3..m-1 heavy
    elements, m the electrons); a, z, x, their atomic masses, charges and
    mass fractions (index 0 unused); cl[s, t], the Coulomb logarithms; ap,
    at and ax[s, t], set in place. x[m] is overwritten with the electrons'
    mass fraction. None is returned.
    """
    # Question (see the page): no check that m <= MMAX.
    indx = np.zeros(NMAX + 1, dtype=int)
    c = np.zeros(MMAX + 1)
    xx = np.zeros((MMAX + 1, MMAX + 1))
    y = np.zeros((MMAX + 1, MMAX + 1))
    yy = np.zeros((MMAX + 1, MMAX + 1))
    k = np.zeros((MMAX + 1, MMAX + 1))
    nu = np.zeros(NMAX + 1)
    alpha = np.zeros(NMAX + 1)
    gamma = np.zeros((NMAX + 1, NMAX + 1))
    delta = np.zeros((NMAX + 1, NMAX + 1))
    ga = np.zeros(NMAX + 1)

    ko = 2.0
    n = 2 * m + 2
    # Concentrations, normalised to the electrons.
    for i in range(1, m + 1):
        c[i] = 0.0
    cc = 0.0
    ac = 0.0
    temp = 0.0
    for i in range(1, m):
        temp = temp + z[i] * x[i] / a[i]
    for i in range(1, m):
        c[i] = x[i] / a[i] / temp
    c[m] = 1.0
    for i in range(1, m + 1):
        cc = cc + c[i]
        ac = ac + a[i] * c[i]
    # Question (see the page): the caller's x[m] is overwritten.
    x[m] = a[m] / ac
    # Mass ratios and collision coefficients.
    for i in range(1, m + 1):
        for j in range(1, m + 1):
            xx[i, j] = a[j] / (a[i] + a[j])
            y[i, j] = a[i] / (a[i] + a[j])
            yy[i, j] = 3.0 * y[i, j] + 1.3 * xx[i, j] * a[j] / a[i]
            k[i, j] = (1.0 * cl[i, j]
                       * math.sqrt(a[i] * a[j] / (a[i] + a[j])) * c[i] * c[j]
                       * z[i]**2 * z[j]**2)
    # Right-hand sides: alpha (pressure), nu (temperature), gamma
    # (concentrations).
    for i in range(1, m + 1):
        alpha[i] = c[i] / cc
        nu[i] = 0.0
        for j in range(1, m + 1):
            gamma[i, j] = 0.0
        for j in range(1, m + 1):
            if j != 2 and j != m:
                gamma[i, j] = -c[j] / cc + c[2] / cc * z[j] * c[j] / z[2] / c[2]
                if j == i:
                    gamma[i, j] = gamma[i, j] + 1.0
                if i == 2:
                    gamma[i, j] = gamma[i, j] - z[j] * c[j] / z[2] / c[2]
                gamma[i, j] = gamma[i, j] * c[i] / cc
        for j in range(m + 1, n + 1):
            gamma[i, j] = 0.0
    for i in range(m + 1, n - 1):
        alpha[i] = 0.0
        nu[i] = 2.5 * c[i - m] / cc
        for j in range(1, n + 1):
            gamma[i, j] = 0.0
    alpha[n - 1] = 0.0
    nu[n - 1] = 0.0
    for j in range(1, n + 1):
        gamma[n - 1, j] = 0.0
    alpha[n] = 0.0
    nu[n] = 0.0
    for j in range(1, n + 1):
        gamma[n, j] = 0.0
    # The matrix delta.
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            delta[i, j] = 0.0
    for i in range(1, m + 1):            # momentum equations
        for j in range(1, m + 1):
            if j == i:
                for l in range(1, m + 1):
                    if l != i:
                        delta[i, j] = delta[i, j] - k[i, l]
            else:
                delta[i, j] = k[i, j]
        for j in range(m + 1, n - 1):
            if j - m == i:
                for l in range(1, m + 1):
                    if l != i:
                        delta[i, j] = delta[i, j] + 0.6 * xx[i, l] * k[i, l]
            else:
                delta[i, j] = -0.6 * y[i, j - m] * k[i, j - m]
        delta[i, n - 1] = c[i] * z[i]
        delta[i, n] = -c[i] * a[i]
    for i in range(m + 1, n - 1):        # heat-flow equations
        for j in range(1, m + 1):
            if j == i - m:
                for l in range(1, m + 1):
                    if l != i - m:
                        delta[i, j] = delta[i, j] + 1.5 * xx[i - m, l] * k[i - m, l]
            else:
                delta[i, j] = -1.5 * xx[i - m, j] * k[i - m, j]
        for j in range(m + 1, n - 1):
            if j - m == i - m:
                for l in range(1, m + 1):
                    if l != i - m:
                        delta[i, j] = (delta[i, j] - y[i - m, l] * k[i - m, l]
                                       * (1.6 * xx[i - m, l] + yy[i - m, l]))
                delta[i, j] = delta[i, j] - 0.8 * k[i - m, i - m]
            else:
                delta[i, j] = (2.7 * k[i - m, j - m] * xx[i - m, j - m]
                               * y[i - m, j - m])
        delta[i, n - 1] = 0.0
        delta[i, n] = 0.0
    for j in range(1, m + 1):            # current neutrality
        delta[n - 1, j] = c[j] * z[j]
    for j in range(m + 1, n + 1):
        delta[n - 1, j] = 0.0
    for j in range(1, m + 1):            # zero mass flow
        delta[n, j] = c[j] * a[j]
    for j in range(m + 1, n + 1):
        delta[n, j] = 0.0
    # Factorise once, solve for each right-hand side.
    d = 0.0
    d = ludcmp(delta, n, NMAX, indx, d)
    lubksb(delta, n, NMAX, indx, alpha)
    lubksb(delta, n, NMAX, indx, nu)
    for j in range(1, n + 1):
        for i in range(1, n + 1):
            ga[i] = gamma[i, j]
        lubksb(delta, n, NMAX, indx, ga)
        for i in range(1, n + 1):
            gamma[i, j] = ga[i]
    # Question (see the page): the factor ko*ac*cc against the comment's p/K0.
    for i in range(1, m + 1):
        alpha[i] = alpha[i] * ko * ac * cc
        nu[i] = nu[i] * ko * ac * cc
        for j in range(1, m + 1):
            gamma[i, j] = gamma[i, j] * ko * ac * cc
    for i in range(1, m + 1):
        ap[i] = alpha[i]
        at[i] = nu[i]
        for j in range(1, m + 1):
            ax[i, j] = gamma[i, j]
