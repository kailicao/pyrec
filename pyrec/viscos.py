"""Python version of viscos.f (YREC 09c3b9b): the kinematic viscosity and
radiative thermal diffusivity of every shell."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import scrtch_mod
from pyrec import temp_mod

# DATA: the atomic mass unit (g), and the mass and charge of HCOMP's first
# 11 species (index 0 unused; species 3, the total Z, has charge 0).
# Question (see the page): written here, not taken from the constants;
# some masses differ from standard values.
AMU = 1.6605655e-24
WEIGHT = np.array([0.0, 1.007825, 4.0026, 1.0, 3.01603, 12.0,
                   13.00335, 14.00307, 15.0011, 15.99491, 16.99491,
                   17.99491])
Z = np.array([0.0, 1.0, 2.0, 0.0, 2.0, 6.0, 6.0, 7.0, 7.0, 8.0,
              8.0, 8.0])

# The Fortran's local arrays (SAVE).
viscx = np.zeros(12)
viscy = np.zeros(12)
compn = np.zeros(12)


def viscos(hcomp, hd, ht, m):
    """For every shell, the kinematic viscosity (radiative, Ledoux 1958,
    plus molecular, Spitzer 1962) into temp_mod.visc, and the radiative
    thermal diffusivity into temp_mod.thdif; the Endal & Sofia values are
    computed for comparison and discarded.

    Translates: viscos.f#VISCOS (09c3b9b)

    Arguments: hcomp[species, shell], log rho hd and log T ht (index 0
    unused); m, the number of shells. Reads the opacity scrtch_mod.so and
    temp_mod.cp (from physic). None is returned.
    """
    cln = constants_mod.cln
    so = scrtch_mod.so
    visc = temp_mod.visc
    thdif = temp_mod.thdif
    cp = temp_mod.cp
    # Question (see the page): the comment says convective zones are
    # skipped, but every shell is computed.
    for im in range(1, m + 1):              # loop 100
        o = so[im]
        zf = 0.0
        comsum = 0.0
        for i in range(1, 12):              # loop 10
            compn[i] = hcomp[i, im] / WEIGHT[i]
            # Question (see the page): species 3 (total Z) enters comsum.
            comsum = comsum + compn[i]
            zf = zf + compn[i] * Z[i]
        zf = zf / comsum
        t = math.exp(cln * ht[im])
        t2 = t**2
        rho = math.exp(cln * hd[im])
        # Question (see the page): mean charge per ion, not per unit mass.
        elecn = zf * rho / AMU
        viscr = 6.7282653e-26 * t2 * t2 / (o * rho**2)
        viscr2 = 3.36e-25 * t2 * t2 / (hcomp[1, im] + 1.0) / rho**2
        factor = min(1.0, 4.2e5 / t)
        lambd1 = 9.424536845 + 0.5 * math.log(t2 * t * factor / (elecn * zf**2))
        c = 3.125e-15 * math.sqrt(t) * t2
        vismol = 0.0
        for i in range(1, 12):              # loop 30
            if i == 3:
                continue
            con = (c * compn[i] * math.sqrt(WEIGHT[i])
                   / ((lambd1 - math.log(Z[i])) * Z[i]**2))
            total = 0.0                     # SUM
            for j in range(1, 12):          # loop 20
                if j == 3:
                    continue
                total = total + compn[j] * Z[j]**2 * math.sqrt(
                    (WEIGHT[i] + WEIGHT[j]) / WEIGHT[j])
            if abs(total) < 1.0e-38:
                continue
            viscx[i] = con / total
            if viscx[i] > 0.0:
                vismol = vismol + viscx[i]
        vismol = vismol / rho
        visc[im] = viscr + vismol
        # Endal & Sofia, for comparison only (question: never used).
        vismo2 = 0.0
        con = 2.21e-15 * math.sqrt(t) * t2
        for i in range(1, 12):              # loop 40
            if i == 3:
                continue
            viscy[i] = con * hcomp[i, im] / (math.sqrt(WEIGHT[i]) * Z[i]**4
                                             * (lambd1 - math.log(Z[i])))
            if viscy[i] > 0.0:
                vismo2 = vismo2 + viscy[i]
        vismo2 = vismo2 / rho / comsum
        visc2 = viscr2 + vismo2             # noqa: F841 (unused, as in Fortran)
        thdif[im] = (1.6e1 * constants_mod.cc13 * 5.669e-5 * t * t2
                     / (o * rho**2 * cp[im]))
