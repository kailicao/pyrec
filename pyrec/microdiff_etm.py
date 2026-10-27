"""Python version of microdiff_etm.f (YREC 09c3b9b): "equal to model",
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
    # closure (the same block at the centre, the interior and the envelope).
    # Question (see the page): fails if hcomp[3, i] = 0.
    zz2 = zz / hcomp[3, i]
    hcomp[3, i] = zz
    for j in range(5, 12):
        hcomp[j, i] = zz2 * hcomp[j, i]
    hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]


def microdiff_etm(dt, er, edelx, edelz, edelli, ibegin, iend, npt, hcomp,
                  hqpr, hru, hs1, htu, m, stot, nlight, cslight):
    """Interpolate the changes in X, Z and the light elements from the
    equal grid back to the model points, apply them to hcomp (helium by
    closure), and convert the variables MICRODIFF_SETUP changed back to cgs.

    Translates: microdiff_etm.f#MICRODIFF_ETM (09c3b9b)

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
    # The centre (points ibegin..1): the change at the first grid point.
    xminim = 0.0
    for i in range(ibegin, 0, -1):                       # loop 10
        hcomp[1, i] = max(hcomp[1, i] + edelx[1], xminim)
    # Question (see the page): which guards are intended where.
    if ldifz:
        for i in range(ibegin, 0, -1):
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            zz = min(hcomp[3, i] + edelz[1], zmax)
            _rescale_metals(hcomp, i, zz)
    else:
        for i in range(ibegin, 0, -1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
    if ldifli:
        for kk in range(1, nlight + 1):
            xminim = 0.0
            for i in range(ibegin, 0, -1):               # loop 15
                ii = cslight[kk]
                hcomp[ii, i] = max(hcomp[ii, i] + edelli[kk, 1], xminim)
    # The radiative interior: four-point Lagrange interpolation.
    jmin = 2
    for i in range(ibegin + 1, iend):                    # loop 20
        for j in range(jmin, npt + 1):                   # loop 30
            if er[j] >= hru[i]:
                k0 = max(j - 2, 1)
                k0 = min(k0, npt - 3)
                # JVS fix for NPT = 3?
                # Question (see the page): k0 = -1 for npt = 2 is not caught.
                if k0 == 0:
                    k0 = 1
                jmin = j
                break                                    # GOTO 40
        else:
            k0 = npt - 3
            jmin = npt
        # label 40
        for k in range(1, 5):                            # loop 50
            tabler[k] = er[k0 + k - 1]
        radmod = hru[i]
        intrp2(tabler, facinterp, radmod)
        dxmod = (facinterp[1] * edelx[k0] + facinterp[2] * edelx[k0 + 1]
                 + facinterp[3] * edelx[k0 + 2] + facinterp[4] * edelx[k0 + 3])
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + dxmod, xmax)
        if ldifz:
            zmax = 1.0 - hcomp[1, i] - hcomp[4, i]
            dzmod = (facinterp[1] * edelz[k0] + facinterp[2] * edelz[k0 + 1]
                     + facinterp[3] * edelz[k0 + 2]
                     + facinterp[4] * edelz[k0 + 3])
            zz = min(hcomp[3, i] + dzmod, zmax)
            _rescale_metals(hcomp, i, zz)
        else:
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
        if ldifli:
            for kk in range(1, nlight + 1):
                ii = cslight[kk]
                dxmod = (facinterp[1] * edelli[kk, k0]
                         + facinterp[2] * edelli[kk, k0 + 1]
                         + facinterp[3] * edelli[kk, k0 + 2]
                         + facinterp[4] * edelli[kk, k0 + 3])
                # Question (see the page): no floor at zero here.
                hcomp[ii, i] = hcomp[ii, i] + dxmod
    # The envelope (points iend..m): the change at the last grid point.
    for i in range(iend, m + 1):                         # loop 60
        xmax = 1.0 - hcomp[3, i] - hcomp[4, i]
        hcomp[1, i] = min(hcomp[1, i] + edelx[npt], xmax)
    if ldifz:
        zmin = 0.0
        for i in range(iend, m + 1):
            zz = max(hcomp[3, i] + edelz[npt], zmin)
            _rescale_metals(hcomp, i, zz)
    else:
        for i in range(iend, m + 1):
            hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
    if ldifli:
        for kk in range(1, nlight + 1):
            for i in range(iend, m + 1):                 # loop 65
                ii = cslight[kk]
                hcomp[ii, i] = hcomp[ii, i] + edelli[kk, npt]
    # Back to cgs units.
    for i in range(1, m + 1):                            # loop 70
        hru[i] = hru[i] / confac_mod.con_rad
        htu[i] = htu[i] / confac_mod.con_temp
        hs1[i] = hs1[i] / confac_mod.con_mass
        hqpr[i] = hqpr[i] * confac_mod.con_rad
    dt = dt * confac_mod.con_time
    stot = stot / confac_mod.con_mass
    return dt, stot
