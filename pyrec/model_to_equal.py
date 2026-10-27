"""Python version of model_to_equal.f (YREC 09c3b9b): the diffusion
coefficients and abundances of the older settling routine (GRSETT) on an
equally spaced grid in radius."""

import math

import numpy as np

from pyrec import dimensions_mod
from pyrec import gravez_mod
from pyrec import nmlpar_mod
from pyrec.interp import interp

JSON = dimensions_mod.json

# The Fortran's local arrays (SAVE).
facinterp = np.zeros(5)
facderiv = np.zeros(5)
tabler = np.zeros(5)


def _lagrange(f, k0):
    # FACINTERP(1)*F(K0) + ... + FACINTERP(4)*F(K0+3), in the Fortran's order.
    return (facinterp[1] * f[k0] + facinterp[2] * f[k0 + 1]
            + facinterp[3] * f[k0 + 2] + facinterp[4] * f[k0 + 3])


def _find_k0(r, jmin, iend, m, hru):
    # Labels 20/30 (and 80/90): the first model point J >= JMIN at or above
    # radius r gives the four interpolation points K0..K0+3, kept inside
    # 1..M; returns (k0, jmin).
    for j in range(jmin, iend + 1):
        if hru[j] >= r:
            k0 = max(j - 2, 1)
            k0 = min(k0, m - 3)
            return k0, j
    return m - 3, m


def model_to_equal(cod1, cod2, hcomp, hru, hs1, ibegin, iend, m, qcod1x,
                   qcod2x, dr, ecod1, ecod1_h, ecod2_h, em, em_h, eqcod1x_h,
                   eqcod2x_h, er, ex, ex_h, npt, cod1z, cod2z, qcod1z,
                   qcod2z):
    """Interpolate the settling coefficients and abundances of GRSETT from
    the model points (ibegin..iend) to an equally spaced grid in radius, at
    the grid points and the midpoints (_h).

    Translates: model_to_equal.f#MODEL_TO_EQUAL (09c3b9b)

    Arguments: the coefficients cod1, cod2, their derivatives in X qcod1x,
    qcod2x, and the metal ones cod1z, cod2z, qcod1z, qcod2z (from
    SETUP_GRSETT); hcomp[k, j], the radii hru and masses hs1 (Bahcall-Loeb
    units); m shells. dr and npt are passed in with any value. Arrays have
    index 0 unused.

    Sets in place: em, ex, ecod1 on the grid points; em_h, ecod1_h, ecod2_h,
    eqcod1x_h, eqcod2x_h, ex_h on the midpoints; er (the grid radii, after
    holding the midpoints' radii); with LDIFZ, the metal arrays of
    gravez_mod.

    Returns (dr, npt): the grid spacing and the number of grid points.
    """
    gz = gravez_mod
    ldifz = nmlpar_mod.ldifz
    # The minimum model point spacing in radius, over the whole region.
    drtot = hru[iend] - hru[ibegin]
    drmin = drtot
    for i in range(ibegin + 1, iend + 1):
        drmin = min(drmin, hru[i] - hru[i - 1])
    # The minimum number of equally spaced points such that DR <= DRMIN.
    npt = int(drtot / drmin)
    if math.fmod(drtot, drmin) != 0.0:
        npt = npt + 1
    # JvS: trap for models with NPT = 1.
    if npt == 1:
        npt = npt + 1
    npt = min(npt, JSON)
    dr = drtot / float(npt - 1)

    # Equally spaced radii at the zone midpoints (held in ER for now).
    # Question (see the page): ER is overwritten below.
    er[1] = hru[ibegin] + 0.5 * dr
    if npt == 2:
        er[2] = er[1] + dr
    else:
        for i in range(2, npt):  # DO 10 (its debug PRINT can never fire)
            er[i] = er[i - 1] + dr

    # First midpoint: linear interpolation.
    # Question (see the page): the search runs over model indices 2..NPT.
    iu = npt
    for k in range(2, npt + 1):  # DO 15
        if hru[k] >= er[1]:
            iu = k
            break
    fx = (er[1] - hru[iu - 1]) / (hru[iu] - hru[iu - 1])

    def lin(f):
        return f[iu - 1] + fx * (f[iu] - f[iu - 1])

    em_h[1] = lin(hs1)
    ecod1_h[1] = lin(cod1)
    ecod2_h[1] = lin(cod2)
    eqcod1x_h[1] = lin(qcod1x)
    eqcod2x_h[1] = lin(qcod2x)
    ex_h[1] = hcomp[1, iu - 1] + fx * (hcomp[1, iu] - hcomp[1, iu - 1])
    if ldifz:
        gz.ecod1z_h[1] = lin(cod1z)
        gz.ecod2z_h[1] = lin(cod2z)
        gz.eqcod1z_h[1] = lin(qcod1z)
        gz.eqcod2z_h[1] = lin(qcod2z)
        # Question (see the page): row 8, not row 3 (Z).
        gz.ez_h[1] = hcomp[8, iu - 1] + fx * (hcomp[8, iu] - hcomp[8, iu - 1])
    # "Center derivative": FAC1, FAC2, DELR are never used (question 2).
    dr1 = er[1] - hru[iu - 1]
    dr2 = hru[iu] - er[1]
    if dr2 > dr1:
        fac1, fac2, delr = 1.0, dr1 / dr2, 2.0 * dr1
    elif dr1 > dr2:
        fac1, fac2, delr = dr2 / dr1, 1.0, 2.0 * dr2
    else:
        fac1, fac2, delr = 1.0, 1.0, dr1 + dr2

    # Other midpoints: 4-point Lagrange interpolation.
    jmin = ibegin + 1
    for i in range(2, npt):  # DO 50
        k0, jmin = _find_k0(er[i], jmin, iend, m, hru)
        for k in range(1, 5):
            tabler[k] = hru[k0 + k - 1]
        interp(tabler, facinterp, facderiv, er[i])
        em_h[i] = _lagrange(hs1, k0)
        ecod1_h[i] = _lagrange(cod1, k0)
        ecod2_h[i] = _lagrange(cod2, k0)
        eqcod1x_h[i] = _lagrange(qcod1x, k0)
        eqcod2x_h[i] = _lagrange(qcod2x, k0)
        ex_h[i] = _lagrange(hcomp[1, :], k0)
        if ldifz:
            # Metal diffusion, treated as fully ionized iron.
            gz.ecod1z_h[i] = _lagrange(cod1z, k0)
            gz.ecod2z_h[i] = _lagrange(cod2z, k0)
            gz.eqcod1z_h[i] = _lagrange(qcod1z, k0)
            gz.eqcod2z_h[i] = _lagrange(qcod2z, k0)
            gz.ez_h[i] = _lagrange(hcomp[8, :], k0)

    # Equally spaced radii at the zone centers.
    er[1] = hru[ibegin]
    for i in range(2, npt + 1):
        er[i] = er[i - 1] + dr
    # First point: by definition, at the starting point.
    em[1] = hs1[ibegin]
    ex[1] = hcomp[1, ibegin]
    ecod1[1] = cod1[ibegin]
    if ldifz:
        gz.ez[1] = hcomp[8, ibegin]
        gz.ecod1z[1] = cod1z[ibegin]
    # Other points: 4-point Lagrange interpolation.
    jmin = ibegin + 1
    for i in range(2, npt):  # DO 70
        k0, jmin = _find_k0(er[i], jmin, iend, m, hru)
        for k in range(1, 5):
            tabler[k] = hru[k0 + k - 1]
        interp(tabler, facinterp, facderiv, er[i])
        em[i] = _lagrange(hs1, k0)
        ecod1[i] = _lagrange(cod1, k0)
        ex[i] = _lagrange(hcomp[1, :], k0)
        if ldifz:
            gz.ecod1z[i] = _lagrange(cod1z, k0)
            gz.ez[i] = _lagrange(hcomp[8, :], k0)
    # Last point: by definition, at the ending point.
    em[npt] = hs1[iend]
    ex[npt] = hcomp[1, iend]
    ecod1[npt] = cod1[iend]
    if ldifz:
        gz.ez[npt] = hcomp[8, iend]
        gz.ecod1z[npt] = cod1z[iend]
    return dr, npt
