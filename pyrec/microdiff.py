"""Python version of microdiff.f (YREC 09c3b9b): the driver of the newer
gravitational settling (Bahcall and Loeb 1989) of helium, metals and light
elements."""

import numpy as np

from pyrec import dimensions_mod
from pyrec import nmlpar_mod
from pyrec.microdiff_etm import microdiff_etm
from pyrec.microdiff_mte import microdiff_mte
from pyrec.microdiff_run import microdiff_run
from pyrec.microdiff_setup import microdiff_setup

JSON = dimensions_mod.json
NLIGHT = 3

# DATA: the light elements (6Li, 7Li, 9Be), their HCOMP rows, atomic
# weights and charges (index 0 unused).
CSLIGHT = np.array([0, 13, 14, 15])
ASLIGHT = np.array([0.0, 6.015, 7.016, 9.012])
ZSLIGHT = np.array([0.0, 3.0, 3.0, 4.0])

# The Fortran's local arrays and scalars (SAVE).
hru = np.zeros(JSON + 1)
htu = np.zeros(JSON + 1)
oden = np.zeros(JSON + 1)
otem = np.zeros(JSON + 1)
em = np.zeros(JSON + 1)
er = np.zeros(JSON + 1)
eden = np.zeros(JSON + 1)
etem = np.zeros(JSON + 1)
ehqpr = np.zeros(JSON + 1)
edelr = np.zeros(JSON + 1)
ex = np.zeros(JSON + 1)
ey = np.zeros(JSON + 1)
ez = np.zeros(JSON + 1)
eli = np.zeros((NLIGHT + 1, JSON + 1))
espec = np.zeros((4, JSON + 1))
edxdr = np.zeros(JSON + 1)
em_h = np.zeros(JSON + 1)
er_h = np.zeros(JSON + 1)
eden_h = np.zeros(JSON + 1)
etem_h = np.zeros(JSON + 1)
ehqpr_h = np.zeros(JSON + 1)
edelr_h = np.zeros(JSON + 1)
ex_h = np.zeros(JSON + 1)
ey_h = np.zeros(JSON + 1)
ez_h = np.zeros(JSON + 1)
eli_h = np.zeros((NLIGHT + 1, JSON + 1))
espec_h = np.zeros((4, JSON + 1))
edxdr_h = np.zeros(JSON + 1)
# Question (see the page): kept from an earlier call when not set (edelx
# when LDIFY is off).
edelx = np.zeros(JSON + 1)
edelz = np.zeros(JSON + 1)
edelli = np.zeros((NLIGHT + 1, JSON + 1))
ibegin = 0
iend = 0
lallcz = False
dr = 0.0
npt = 0


def microdiff(dt, hcomp, hqpr, hr, hd, hs1, ht, lc, m, stot):
    """Settle hydrogen (LDIFY), the metals as iron (LDIFZ) and the light
    elements (LDIFLI) over one step: move the model onto an equally spaced
    grid in radius, diffuse each species there, and map the changes back.

    Translates: microdiff.f#MICRODIFF (09c3b9b)

    Arguments: the step dt (s); the model's hcomp[k, j], hqpr (dlnP/dr),
    hr (log r), hd (log rho), hs1 (mass, g), ht (log T) and convective
    flags lc, with m shells; the total mass stot (g). Arrays have index 0
    unused.

    Sets in place: hcomp. hs1 and hqpr are converted to Bahcall and Loeb's
    units and back.

    Returns (dt, stot): converted and converted back by the callees, so
    equal to the values passed in up to rounding.
    """
    global ibegin, iend, lallcz, dr, npt
    dt, stot, ibegin, iend, lallcz = microdiff_setup(
        dt, hqpr, hr, hd, hs1, ht, lc, m, stot, hcomp, hru, htu,
        ibegin, iend, lallcz, oden, otem)
    # Skip settling for fully convective models.
    if lallcz:
        return dt, stot
    # Transform to an equally spaced grid in radius (E: at the grid points,
    # E..._H: at the midpoints).
    dr, npt = microdiff_mte(NLIGHT, CSLIGHT, hcomp, hqpr, hru, hs1, ibegin,
                            iend, m, dr, npt, oden, otem, em, er, eden, etem,
                            ehqpr, edelr, ex, ey, ez, eli, em_h, er_h, eden_h,
                            etem_h, ehqpr_h, edelr_h, ex_h, ey_h, ez_h, eli_h)
    # The three species of the Thoul calculation: 1 = X, 2 = Y, 3 = metal.
    for i in range(1, npt):
        espec[1, i] = ex[i]
        espec[2, i] = ey[i]
        espec_h[1, i] = ex_h[i]
        espec_h[2, i] = ey_h[i]
    espec[1, npt] = ex[npt]
    espec[2, npt] = ey[npt]
    # Diffuse hydrogen, with the metals as iron for the third species.
    if nmlpar_mod.ldify:
        j = 1
        thea = 55.86
        thez = 26.0
        for i in range(1, npt):
            espec[3, i] = ez[i]
            espec_h[3, i] = ez_h[i]
        espec[3, npt] = ez[npt]
        microdiff_run(dr, dt, stot, npt, em, er, eden, etem, ehqpr, edelr,
                      espec, edxdr, em_h, er_h, eden_h, etem_h, ehqpr_h,
                      edelr_h, espec_h, edxdr_h, thea, thez, j)
        for i in range(1, npt + 1):
            edelx[i] = espec[1, i]
        # Restore the original hydrogen vectors.
        for i in range(1, npt):
            espec[1, i] = ex[i]
            espec_h[1, i] = ex_h[i]
        espec[1, npt] = ex[npt]
    # Diffuse the metals, as iron.
    # Question (see the page): iron for all the metals.
    if nmlpar_mod.ldifz:
        j = 3
        thea = 55.86
        thez = 26.0
        for i in range(1, npt):
            espec[3, i] = ez[i]
            espec_h[3, i] = ez_h[i]
        espec[3, npt] = ez[npt]
        microdiff_run(dr, dt, stot, npt, em, er, eden, etem, ehqpr, edelr,
                      espec, edxdr, em_h, er_h, eden_h, etem_h, ehqpr_h,
                      edelr_h, espec_h, edxdr_h, thea, thez, j)
        for i in range(1, npt + 1):
            edelz[i] = espec[3, i]
    # Diffuse the light elements, each in turn as the third species.
    if nmlpar_mod.ldifli:
        j = 3
        for ii in range(1, NLIGHT + 1):
            thea = ASLIGHT[ii]
            thez = ZSLIGHT[ii]
            for i in range(1, npt):
                espec[3, i] = eli[ii, i]
                espec_h[3, i] = eli_h[ii, i]
            # Question (see the page): the Fortran's I after DO I=1,NPT-1,
            # that is max(npt, 1), stands for NPT.
            i = max(npt, 1)
            espec[3, i] = eli[ii, npt]
            microdiff_run(dr, dt, stot, npt, em, er, eden, etem, ehqpr,
                          edelr, espec, edxdr, em_h, er_h, eden_h, etem_h,
                          ehqpr_h, edelr_h, espec_h, edxdr_h, thea, thez, j)
            for i in range(1, npt + 1):
                edelli[ii, i] = espec[3, i]
    # Back to the model grid; helium by closure.
    dt, stot = microdiff_etm(dt, er, edelx, edelz, edelli, ibegin, iend, npt,
                             hcomp, hqpr, hru, hs1, htu, m, stot, NLIGHT,
                             CSLIGHT)
    return dt, stot
