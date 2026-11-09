"""Python version of dcoeft.f (YREC 7be4133): the tridiagonal system for
one implicit step of angular-momentum transport on the even grid."""

import numpy as np

from pyrec import difad_mod
from pyrec import nmlpar_mod
from pyrec import tridi_mod


def dcoeft(ecod: np.ndarray, dr: float, dt: float, ei: np.ndarray, ej: np.ndarray, ew: np.ndarray, ntot: int,
           wind1: float, wind2: float, lfix: bool) -> None:
    """Set up the implicit tridiagonal system for Omega (diffusion, and with
    ldifad a centred advective term), with no flux at the centre, no flux
    or a fixed Omega at the surface, and the wind's loss of angular
    momentum in the surface row.

    Translates: dcoeft.f#DCOEFT (7be4133)

    Arguments: ecod, the diffusion coefficients; dr, the grid step in chi;
    dt, the step; ei, ej, ew, the moments of inertia, angular momenta and
    Omega; ntot, the number of points; wind1, wind2, the wind's change of
    Omega (explicit, implicit); lfix, a fixed surface Omega.

    Returns nothing: sets tridi_mod.a, b, c, d and gama[ntot] (the wind's
    Delta J, read by TRIDIA).

    """

    a = tridi_mod.a
    b = tridi_mod.b
    c = tridi_mod.c
    d = tridi_mod.d

    # The starting Omega, and the wind's Delta J in the surface row
    # (passed to TRIDIA in gama[ntot]) (Fortran lines 78–83).
    fact0 = dt / dr
    for i in range(1, ntot + 1):  # label 10
        d[i] = ew[i]
    tridi_mod.gama[ntot] = -0.5 * (wind1 + wind2) * ei[ntot]
    d[ntot] = d[ntot] * (1.0 + tridi_mod.gama[ntot] / ej[ntot])

    if not nmlpar_mod.ldifad:
        # Diffusion only; no flux below the first shell (Fortran lines
        # 84–91).
        fact = fact0 / ei[1]
        a[1] = 0.0
        b[1] = 1.0 + fact * ecod[2]
        c[1] = -fact * ecod[2]

        # The interior rows (Fortran lines 92–97).
        for i in range(2, ntot):  # label 20
            fact = fact0 / ei[i]
            a[i] = -fact * ecod[i]
            b[i] = 1.0 + fact * (ecod[i] + ecod[i + 1])
            c[i] = -fact * ecod[i + 1]

        # The surface row: no flux, or a fixed Omega (Fortran lines
        # 98–107).
        fact = fact0 / ei[ntot]
        if not lfix:
            a[ntot] = -fact * ecod[ntot]
            b[ntot] = 1.0 + fact * ecod[ntot]
        else:
            a[ntot] = 0.0
            b[ntot] = 1.0
        c[ntot] = 0.0

    else:
        # Diffusion (ecod4) and advection (ecod3); unused? (SECULR calls
        # DCOEFT only without ldifad). The first row (Fortran lines
        # 108–119).
        ecod3 = difad_mod.ecod3
        ecod4 = difad_mod.ecod4
        fact = fact0 / ei[1]
        facta = 0.5 * dt / ei[1]
        a[1] = 0.0
        b[1] = 1.0 + fact * ecod4[2] - facta * ecod3[2]
        c[1] = -fact * ecod4[2] - facta * ecod3[2]

        # The interior rows (Fortran lines 120–127).
        for i in range(2, ntot):
            fact = fact0 / ei[i]
            facta = 0.5 * dt / ei[i]
            a[i] = -fact * ecod4[i] + facta * ecod3[i]
            b[i] = (1.0 + fact * (ecod4[i] + ecod4[i + 1])
                    + facta * (ecod3[i] - ecod3[i + 1]))
            c[i] = -fact * ecod4[i + 1] - facta * ecod3[i + 1]

        # The surface row (Fortran lines 128–138).
        fact = fact0 / ei[ntot]
        facta = 0.5 * dt / ei[ntot]
        if not lfix:
            a[ntot] = -fact * ecod4[ntot] + facta * ecod3[ntot]
            b[ntot] = 1.0 + fact * ecod4[ntot] + facta * ecod3[ntot]
        else:
            a[ntot] = 0.0
            b[ntot] = 1.0
        c[ntot] = 0.0
