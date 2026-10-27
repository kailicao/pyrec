"""Python version of equal_to_model.f (YREC 09c3b9b): the settling
changes from the equally spaced grid back to the model points."""

import numpy as np

from pyrec import confac_mod
from pyrec import gravez_mod
from pyrec import nmlpar_mod
from pyrec.intrp2 import intrp2


def _scale_metals(hcomp, i, zz):
    """Set Z of shell i to zz, scale the metal species (rows 5-11) by
    Z_new/Z_old, and give helium the rest."""
    # Question (see the page): a shell with Z = 0 divides by zero.
    zz2 = zz / hcomp[3, i]
    hcomp[3, i] = zz
    for j in range(5, 11 + 1):
        hcomp[j, i] = zz2 * hcomp[j, i]
    hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]


def equal_to_model(dt, er, ex, ibegin, iend, npt, hcomp, hqpr, hru, hs1,
                   htu, m, stot):
    """Apply the changes in X (and, with nmlpar_mod.ldifz, in Z) found on
    the equally spaced grid to the model points, by four-point Lagrange
    interpolation inside the grid and the end values outside it, within
    the composition limits; then convert the model back from the
    solver's units.

    Translates: equal_to_model.f#EQUAL_TO_MODEL (09c3b9b)

    Arguments: er, ex (grid radii and changes in X, npt points); ibegin,
    iend (the model points bounding the grid); m (number of shells).
    Arrays changed in place: hcomp, hqpr, hru, hs1, htu.

    Returns (dt, stot), back in physical units.
    """
    tabler = np.zeros(4 + 1)
    facinterp = np.zeros(4 + 1)
    ez = gravez_mod.ez

    # Below the grid: the change at the first grid point.
    xminim = 0.0
    for i in range(ibegin, 0, -1):  # label 10
        hcomp[1, i] = max(hcomp[1, i] + ex[1], xminim)
    # Metals sink and hydrogen rises, so the failsafes are opposite.
    if nmlpar_mod.ldifz:
        for i in range(ibegin, 0, -1):
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            zz = min(hcomp[3, i] + ez[1], zmax)
            _scale_metals(hcomp, i, zz)
    else:
        for i in range(ibegin, 0, -1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # Inside the grid: four-point Lagrange interpolation.
    jmin = 2
    for i in range(ibegin + 1, iend):  # label 20
        for j in range(jmin, npt + 1):  # label 30
            if j == 0:
                # Question (see the page): never true.
                print(' line 47 etm')
            if er[j] >= hru[i]:
                k0 = max(j - 2, 1)
                k0 = min(k0, npt - 3)
                if k0 == 0:
                    k0 = 1
                jmin = j
                break  # GOTO 40
        else:
            k0 = npt - 3
            jmin = npt
        # label 40
        for k in range(1, 4 + 1):  # label 50
            tabler[k] = er[k0 + k - 1]
        radmod = hru[i]
        intrp2(tabler, facinterp, radmod)
        dxmod = (facinterp[1] * ex[k0] + facinterp[2] * ex[k0 + 1]
                 + facinterp[3] * ex[k0 + 2] + facinterp[4] * ex[k0 + 3])
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + dxmod, xmax)
        if nmlpar_mod.ldifz:
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            dzmod = (facinterp[1] * ez[k0] + facinterp[2] * ez[k0 + 1]
                     + facinterp[3] * ez[k0 + 2]
                     + facinterp[4] * ez[k0 + 3])
            zz = min(hcomp[3, i] + dzmod, zmax)
            _scale_metals(hcomp, i, zz)
        else:
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # Above the grid: the change at the last grid point.
    for i in range(iend, m + 1):  # label 60
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + ex[npt], xmax)
    if nmlpar_mod.ldifz:
        zmin = 0.0
        for i in range(iend, m + 1):
            zz = max(hcomp[3, i] + ez[npt], zmin)
            _scale_metals(hcomp, i, zz)
    else:
        for i in range(iend, m + 1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # Back to physical units.
    for i in range(1, m + 1):  # label 70
        hru[i] = hru[i] / confac_mod.con_rad
        htu[i] = htu[i] / confac_mod.con_temp
        hs1[i] = hs1[i] / confac_mod.con_mass
        hqpr[i] = hqpr[i] * confac_mod.con_rad
    dt = dt * confac_mod.con_time
    stot = stot / confac_mod.con_mass
    return dt, stot
