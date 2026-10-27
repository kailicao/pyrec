"""Python version of codiff.f (YREC 09c3b9b): diffusion coefficients for
rotational mixing."""

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import nmlpar_mod
from pyrec import temp2_mod

JSON = dimensions_mod.json

# Formerly COMMON/ADVEC/; declared, kept by SAVE, and never used.
fadv = np.zeros(JSON + 1)
fadv0 = np.zeros(JSON + 1)


def codiff(hru, m, rmid, cod, cod2, vfc):
    """Diffusion coefficients at the interfaces between shells: cod for
    angular momentum, from the sum of the rotational-instability velocities
    times the radius, and cod2 for composition (with the extra factor f_c
    and per-mechanism weights, or the variable f_c in vfc); options for a
    constant background (LCODM) and for no angular-momentum transport
    (LNOJ).

    Translates: codiff.f#CODIFF (09c3b9b)

    Arguments: hru, the radius of each shell (not logged); m, the number of
    shells; rmid, cod, cod2 (set here, element i at the interface between
    shells i-1 and i; elements 1 of cod, cod2 are not set); vfc, the
    variable f_c. Arrays have index 0 unused and are changed in place.

    Returns nothing.
    """
    nml = nmlpar_mod
    t2 = temp2_mod
    # Question (see the page): the length scale; HLE is set to zero.
    t2.hle[1] = 0.0
    rmid[1] = 0.0
    for i in range(2, m + 1):
        t2.hle[i] = 0.0
        rmid[i] = 0.5 * (hru[i] + hru[i - 1])
    con1 = constants_mod.c4pi * nml.fw
    # Question (see the page): CODM is added to cod only, never to cod2.
    if not nml.lcodm:
        if not nml.lvfc:
            con2 = con1 * nml.fc
            for i in range(2, m + 1):
                cod[i] = con1 * (t2.ves[i] + t2.vgsf[i] + t2.vss[i]) * rmid[i]
                cod2[i] = con2 * (nml.fesc * t2.ves[i] + nml.fgsfc * t2.vgsf[i]
                                  + nml.fssc * t2.vss[i]) * rmid[i]
        else:
            for i in range(2, m + 1):
                cod[i] = con1 * (t2.ves[i] + t2.vgsf[i] + t2.vss[i]) * rmid[i]
                cod2[i] = cod[i] * vfc[i]
    else:
        if not nml.lvfc:
            con2 = con1 * nml.fc
            for i in range(2, m + 1):
                cod[i] = con1 * (nml.codm + (t2.ves[i] + t2.vgsf[i]
                                             + t2.vss[i]) * rmid[i])
                cod2[i] = con2 * (nml.fesc * t2.ves[i] + nml.fgsfc * t2.vgsf[i]
                                  + nml.fssc * t2.vss[i]) * rmid[i]
        else:
            for i in range(2, m + 1):
                cod[i] = con1 * (t2.ves[i] + t2.vgsf[i] + t2.vss[i]) * rmid[i]
                cod2[i] = cod[i] * vfc[i]
                cod[i] = cod[i] + con1 * nml.codm
    # No angular-momentum transport (MHP 9/93).
    if nml.lnoj:
        for i in range(2, m + 1):
            cod[i] = 0.0
