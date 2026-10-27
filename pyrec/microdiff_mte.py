"""Python version of microdiff_mte.f (YREC 09c3b9b): "model to equal",
the model interpolated onto an equally spaced grid in radius."""

import math

import numpy as np

from pyrec import dimensions_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod
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


def microdiff_mte(nlight, cslight, hcomp, hqpr, hru, hs1, ibegin,
                  iend, m, dr, npt, oden, otem, em, er, eden, etem, ehqpr,
                  edelr, ex, ey, ez, eli, em_h, er_h, eden_h, etem_h,
                  ehqpr_h, edelr_h, ex_h, ey_h, ez_h, eli_h):
    """Interpolate the model between shells ibegin and iend onto an equally
    spaced grid in radius, at the grid points (e...) and the midpoints
    between them (e..._h), for MICRODIFF.

    Translates: microdiff_mte.f#MICRODIFF_MTE (09c3b9b)

    Arguments: nlight light elements, at the hcomp rows cslight[1..nlight];
    the model's hcomp[k, j], hqpr, radius hru (Rsun), mass hs1, density
    oden and temperature otem, with m shells; the gradient comes from
    scrtch_mod.sdel[2, j]. dr and npt are passed in with any value.
    Arrays have index 0 unused.

    Sets in place: em, er, eden, etem, ehqpr, edelr, ex, ey, ez, eli[kk, i]
    on the npt grid points, and the same with _h on the npt - 1 midpoints.

    Returns (dr, npt): the grid spacing and the number of grid points.
    """
    sdel = scrtch_mod.sdel
    ldifli = nmlpar_mod.ldifli
    # Find the minimum model point spacing in radius, in the 20 layers below
    # the surface convection zone.
    # Question (see the page): only these 20 layers.
    drtot = hru[iend] - hru[ibegin]
    drmin = drtot
    istart = max(iend - 20, ibegin + 1)
    for i in range(istart, iend + 1):
        drmin = min(drmin, hru[i] - hru[i - 1])
    # The minimum number of equally spaced points such that DR <= DRMIN.
    # Question (see the page): NPT points make NPT - 1 intervals.
    npt = int(drtot / drmin)
    if math.fmod(drtot, drmin) != 0.0:
        npt = npt + 1
    # JvS: trap for models with NPT = 1.
    if npt == 1:
        npt = npt + 1
    # Ensure that the number of points does not exceed JSON.
    halfjson = JSON
    npt = min(npt, halfjson)
    dr = drtot / float(npt - 1)

    # Equally spaced radii at the zone midpoints.
    er_h[1] = hru[ibegin] + 0.5 * dr
    if npt == 2:
        er_h[2] = er_h[1] + dr
    else:
        for i in range(2, npt):  # DO 10 (its debug PRINT can never fire)
            er_h[i] = er_h[i - 1] + dr

    # First midpoint: linear interpolation.
    # Question (see the page): the search runs over model indices 2..NPT.
    iu = npt
    for k in range(2, npt + 1):  # DO 15
        if hru[k] >= er_h[1]:
            iu = k
            break
    fx = (er_h[1] - hru[iu - 1]) / (hru[iu] - hru[iu - 1])

    def lin(f):
        return f[iu - 1] + fx * (f[iu] - f[iu - 1])

    em_h[1] = lin(hs1)
    eden_h[1] = lin(oden)
    etem_h[1] = lin(otem)
    ehqpr_h[1] = lin(hqpr)
    edelr_h[1] = sdel[2, iu - 1] + fx * (sdel[2, iu] - sdel[2, iu - 1])
    ex_h[1] = hcomp[1, iu - 1] + fx * (hcomp[1, iu] - hcomp[1, iu - 1])
    ey_h[1] = hcomp[2, iu - 1] + fx * (hcomp[2, iu] - hcomp[2, iu - 1])
    ez_h[1] = hcomp[3, iu - 1] + fx * (hcomp[3, iu] - hcomp[3, iu - 1])
    if ldifli:
        for kk in range(1, nlight + 1):
            ii = cslight[kk]
            eli_h[kk, 1] = (hcomp[ii, iu - 1]
                            + fx * (hcomp[ii, iu] - hcomp[ii, iu - 1]))

    # Other midpoints: 4-point Lagrange interpolation. JMIN is the uppermost
    # model point above the previous grid point.
    jmin = ibegin + 1
    for i in range(2, npt):  # DO 50
        k0, jmin = _find_k0(er_h[i], jmin, iend, m, hru)
        for k in range(1, 5):
            tabler[k] = hru[k0 + k - 1]
        interp(tabler, facinterp, facderiv, er_h[i])
        em_h[i] = _lagrange(hs1, k0)
        eden_h[i] = _lagrange(oden, k0)
        etem_h[i] = _lagrange(otem, k0)
        ehqpr_h[i] = _lagrange(hqpr, k0)
        edelr_h[i] = _lagrange(sdel[2, :], k0)
        ex_h[i] = _lagrange(hcomp[1, :], k0)
        ey_h[i] = _lagrange(hcomp[2, :], k0)
        ez_h[i] = _lagrange(hcomp[3, :], k0)
        if ldifli:
            for kk in range(1, nlight + 1):
                eli_h[kk, i] = _lagrange(hcomp[cslight[kk], :], k0)

    # Equally spaced radii at the zone centers.
    er[1] = hru[ibegin]
    for i in range(2, npt + 1):
        er[i] = er[i - 1] + dr

    # First point: by definition, at the starting point.
    em[1] = hs1[ibegin]
    eden[1] = oden[ibegin]
    etem[1] = otem[ibegin]
    ehqpr[1] = hqpr[ibegin]
    edelr[1] = sdel[2, ibegin]
    ex[1] = hcomp[1, ibegin]
    ey[1] = hcomp[2, ibegin]
    ez[1] = hcomp[3, ibegin]
    if ldifli:
        for kk in range(1, nlight + 1):
            eli[kk, 1] = hcomp[cslight[kk], ibegin]
    # Other points: 4-point Lagrange interpolation.
    jmin = ibegin + 1
    for i in range(2, npt):  # DO 70
        k0, jmin = _find_k0(er[i], jmin, iend, m, hru)
        for k in range(1, 5):
            tabler[k] = hru[k0 + k - 1]
        interp(tabler, facinterp, facderiv, er[i])
        em[i] = _lagrange(hs1, k0)
        eden[i] = _lagrange(oden, k0)
        etem[i] = _lagrange(otem, k0)
        ehqpr[i] = _lagrange(hqpr, k0)
        edelr[i] = _lagrange(sdel[2, :], k0)
        ex[i] = _lagrange(hcomp[1, :], k0)
        ey[i] = _lagrange(hcomp[2, :], k0)
        ez[i] = _lagrange(hcomp[3, :], k0)
        if ldifli:
            for kk in range(1, nlight + 1):
                eli[kk, i] = _lagrange(hcomp[cslight[kk], :], k0)
    # Last point: by definition, at the ending point.
    em[npt] = hs1[iend]
    eden[npt] = oden[iend]
    etem[npt] = otem[iend]
    ehqpr[npt] = hqpr[iend]
    edelr[npt] = sdel[2, iend]
    ex[npt] = hcomp[1, iend]
    ey[npt] = hcomp[2, iend]
    ez[npt] = hcomp[3, iend]
    if ldifli:
        for kk in range(1, nlight + 1):
            eli[kk, npt] = hcomp[cslight[kk], iend]
    return dr, npt
