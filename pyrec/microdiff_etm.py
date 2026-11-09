"""Python version of microdiff_etm.f (YREC 7be4133): "equal to model",
the changes in abundance mapped back from the equal grid to the model."""

import numpy as np

from pyrec import confac_mod
from pyrec import nmlpar_mod
from pyrec.intrp2 import intrp2

# The Fortran's local arrays (SAVE).
tabler = np.zeros(5)
facinterp = np.zeros(5)


def _rescale_metals(hcomp, i, zz):
    # HCOMP(3,I) set to ZZ, the metal rows 5..11 scaled by ZZ/Z, and Y by
    # closure: the same block at the centre, the interior and the envelope
    # (Fortran lines 53–58, 112–117, 143–148).
    # Suspected bug (see the page): divides by zero if hcomp[3, i] = 0.
    zz2 = zz / hcomp[3, i]
    hcomp[3, i] = zz
    for j in range(5, 12):
        hcomp[j, i] = zz2 * hcomp[j, i]
    hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]


def microdiff_etm(dt: float, er: np.ndarray, edelx: np.ndarray, edelz: np.ndarray, edelli: np.ndarray, ibegin: int,
                  iend: int, npt: int, hcomp: np.ndarray, hqpr: np.ndarray, hru: np.ndarray, hs1: np.ndarray,
                  htu: np.ndarray, m: int, stot: float, nlight: int, cslight: np.ndarray) -> tuple[float, float]:
    """Interpolate the changes in X, Z and the light elements from the
    equal grid back to the model points, apply them to hcomp (helium by
    closure), and convert the variables MICRODIFF_SETUP changed back to cgs.

    Translates: microdiff_etm.f#MICRODIFF_ETM (7be4133)

    Arguments: dt (Bahcall and Loeb units); the grid radii er and the
    changes edelx, edelz, edelli[kk, i] on npt points; the region ibegin,
    iend; the model hcomp[k, j], hqpr, hru, hs1, htu with m shells; stot;
    nlight light elements at the hcomp rows cslight[1..nlight]. Arrays
    have index 0 unused.

    Sets in place: hcomp; hru, htu, hs1, hqpr converted back.

    Returns (dt, stot), converted back to seconds and grams.

    """

    ldifz = nmlpar_mod.ldifz
    ldifli = nmlpar_mod.ldifli

    # The centre (points ibegin..1): the change at the first grid point
    # (Fortran lines 40–45).
    xminim = 0.0
    for i in range(ibegin, 0, -1):                       # loop 10
        hcomp[1, i] = max(hcomp[1, i] + edelx[1], xminim)

    # Metals and helium at the centre (Fortran lines 46–64).
    # Question (see the page): which guards are intended where.
    if ldifz:
        for i in range(ibegin, 0, -1):
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            zz = min(hcomp[3, i] + edelz[1], zmax)
            _rescale_metals(hcomp, i, zz)
    else:
        for i in range(ibegin, 0, -1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # Light elements at the centre, floored at zero (Fortran lines 65–74).
    if ldifli:
        for kk in range(1, nlight + 1):
            xminim = 0.0
            for i in range(ibegin, 0, -1):               # loop 15
                ii = cslight[kk]
                hcomp[ii, i] = max(hcomp[ii, i] + edelli[kk, 1], xminim)

    # The radiative interior: four-point Lagrange interpolation (Fortran
    # lines 75–132).
    jmin = 2
    for i in range(ibegin + 1, iend):                    # loop 20
        # The grid points around the model point (Fortran lines 76–93).
        for j in range(jmin, npt + 1):                   # loop 30
            if er[j] >= hru[i]:
                k0 = max(j - 2, 1)
                k0 = min(k0, npt - 3)
                # Handles the case NPT = 3 (JVS).
                if k0 == 0:
                    k0 = 1
                jmin = j
                break                                    # GOTO 40
        else:
            k0 = npt - 3
            jmin = npt
        # label 40

        # The interpolation factors, and X (Fortran lines 95–105).
        for k in range(1, 5):                            # loop 50
            tabler[k] = er[k0 + k - 1]
        radmod = hru[i]
        intrp2(tabler, facinterp, radmod)
        dxmod = (facinterp[1] * edelx[k0] + facinterp[2] * edelx[k0 + 1]
                 + facinterp[3] * edelx[k0 + 2] + facinterp[4] * edelx[k0 + 3])
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + dxmod, xmax)

        # Metals and helium (Fortran lines 106–120).
        if ldifz:
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            dzmod = (facinterp[1] * edelz[k0] + facinterp[2] * edelz[k0 + 1]
                     + facinterp[3] * edelz[k0 + 2]
                     + facinterp[4] * edelz[k0 + 3])
            zz = min(hcomp[3, i] + dzmod, zmax)
            _rescale_metals(hcomp, i, zz)
        else:
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

        # Light elements (Fortran lines 121–131).
        if ldifli:
            for kk in range(1, nlight + 1):
                ii = cslight[kk]
                dxmod = (facinterp[1] * edelli[kk, k0]
                         + facinterp[2] * edelli[kk, k0 + 1]
                         + facinterp[3] * edelli[kk, k0 + 2]
                         + facinterp[4] * edelli[kk, k0 + 3])
                # Suspected bug (see the page): no floor at zero here.
                hcomp[ii, i] = hcomp[ii, i] + dxmod

    # The envelope (points iend..m): the change at the last grid point
    # (Fortran lines 134–137).
    for i in range(iend, m + 1):                         # loop 60
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + edelx[npt], xmax)

    # Metals and helium in the envelope (Fortran lines 138–154).
    if ldifz:
        zmin = 0.0
        for i in range(iend, m + 1):
            zz = max(hcomp[3, i] + edelz[npt], zmin)
            _rescale_metals(hcomp, i, zz)
    else:
        for i in range(iend, m + 1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # Light elements in the envelope (Fortran lines 155–163).
    # Suspected bug (see the page): no floor at zero here.
    if ldifli:
        for kk in range(1, nlight + 1):
            for i in range(iend, m + 1):                 # loop 65
                ii = cslight[kk]
                hcomp[ii, i] = hcomp[ii, i] + edelli[kk, npt]

    # Back to cgs units; MICRODIFF does not use these (Fortran lines
    # 165–172).
    for i in range(1, m + 1):                            # loop 70
        hru[i] = hru[i] / confac_mod.con_rad
        htu[i] = htu[i] / confac_mod.con_temp
        hs1[i] = hs1[i] / confac_mod.con_mass
        hqpr[i] = hqpr[i] * confac_mod.con_rad
    dt = dt * confac_mod.con_time
    stot = stot / confac_mod.con_mass

    return dt, stot
