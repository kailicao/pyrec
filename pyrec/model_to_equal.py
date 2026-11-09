"""Python version of model_to_equal.f (YREC 7be4133): the diffusion
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


def _lagrange(f: np.ndarray, k0: int) -> float:
    """FACINTERP(1)*F(K0) + ... + FACINTERP(4)*F(K0+3), in the Fortran's
    order.

    """

    return (facinterp[1] * f[k0] + facinterp[2] * f[k0 + 1]
            + facinterp[3] * f[k0 + 2] + facinterp[4] * f[k0 + 3])


def _find_k0(r: float, jmin: int, iend: int, m: int, hru: np.ndarray) -> tuple[int, int]:
    """Labels 20/30 (and 80/90): the first model point J >= JMIN at or
    above radius r gives the four interpolation points K0..K0+3, kept
    inside 1..M; returns (k0, jmin) (Fortran lines 133–146, 219–232).

    """

    for j in range(jmin, iend + 1):
        if hru[j] >= r:
            k0 = max(j - 2, 1)
            k0 = min(k0, m - 3)
            return k0, j
    return m - 3, m


def model_to_equal(cod1: np.ndarray, cod2: np.ndarray, hcomp: np.ndarray, hru: np.ndarray, hs1: np.ndarray, ibegin: int,
                   iend: int, m: int, qcod1x: np.ndarray, qcod2x: np.ndarray, dr: float, ecod1: np.ndarray,
                   ecod1_h: np.ndarray, ecod2_h: np.ndarray, em: np.ndarray, em_h: np.ndarray, eqcod1x_h: np.ndarray,
                   eqcod2x_h: np.ndarray, er: np.ndarray, ex: np.ndarray, ex_h: np.ndarray, npt: int, cod1z: np.ndarray,
                   cod2z: np.ndarray, qcod1z: np.ndarray, qcod2z: np.ndarray) -> tuple[float, int]:
    """Interpolate the settling coefficients and abundances of GRSETT from
    the model points (ibegin..iend) to an equally spaced grid in radius, at
    the grid points and the midpoints (_h).

    Translates: model_to_equal.f#MODEL_TO_EQUAL (7be4133)

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

    # The minimum model point spacing in radius, over the whole region
    # (Fortran lines 55–59).
    drtot = hru[iend] - hru[ibegin]
    drmin = drtot
    for i in range(ibegin + 1, iend + 1):
        drmin = min(drmin, hru[i] - hru[i - 1])

    # The minimum number of equally spaced points such that DR <= DRMIN, at
    # least 4 for the 4-point Lagrange stencils (Fortran lines 62–73).
    npt = int(drtot / drmin)
    if math.fmod(drtot, drmin) != 0.0:
        npt = npt + 1
    npt = max(npt, 4)
    npt = min(npt, JSON)
    dr = drtot / float(npt - 1)

    # Equally spaced radii at the zone midpoints, held in ER for now
    # (Fortran lines 76–85). The NPT == 2 branch cannot be taken.
    er[1] = hru[ibegin] + 0.5 * dr
    if npt == 2:
        er[2] = er[1] + dr
    else:
        for i in range(2, npt):  # DO 10 (its debug PRINT can never fire)
            er[i] = er[i - 1] + dr

    # First midpoint: linear interpolation (Fortran lines 91–109).
    # Suspected bug (see the page): the search runs over model indices
    # 2..NPT.
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
        # Row 8, not row 3 (Z); its values cancel in GRSETT's change.
        gz.ez_h[1] = hcomp[8, iu - 1] + fx * (hcomp[8, iu] - hcomp[8, iu - 1])

    # "Center derivative": FAC1, FAC2, DELR are never used (Fortran lines
    # 111–125).
    dr1 = er[1] - hru[iu - 1]
    dr2 = hru[iu] - er[1]
    if dr2 > dr1:
        fac1, fac2, delr = 1.0, dr1 / dr2, 2.0 * dr1  # noqa: F841  (unused in the Fortran too)
    elif dr1 > dr2:
        fac1, fac2, delr = dr2 / dr1, 1.0, 2.0 * dr2  # noqa: F841  (unused in the Fortran too)
    else:
        fac1, fac2, delr = 1.0, 1.0, dr1 + dr2  # noqa: F841  (unused in the Fortran too)

    # Other midpoints: 4-point Lagrange interpolation (Fortran lines
    # 131–192).
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

    # Equally spaced radii at the zone centers (Fortran lines 195–198).
    er[1] = hru[ibegin]
    for i in range(2, npt + 1):
        er[i] = er[i - 1] + dr

    # First point: by definition, at the starting point (Fortran lines
    # 205–211).
    em[1] = hs1[ibegin]
    ex[1] = hcomp[1, ibegin]
    ecod1[1] = cod1[ibegin]
    if ldifz:
        gz.ez[1] = hcomp[8, ibegin]
        gz.ecod1z[1] = cod1z[ibegin]

    # Other points: 4-point Lagrange interpolation (Fortran lines 217–263).
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

    # Last point: by definition, at the ending point (Fortran lines
    # 266–273).
    em[npt] = hs1[iend]
    ex[npt] = hcomp[1, iend]
    ecod1[npt] = cod1[iend]
    if ldifz:
        gz.ez[npt] = hcomp[8, iend]
        gz.ecod1z[npt] = cod1z[iend]

    return dr, npt
