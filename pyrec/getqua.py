"""Python version of getqua.f (YREC 7be4133): the quadrupole moment of a
rotating star (Sweet 1950, as given in Zahn 1992)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import quadru_mod
from pyrec import tridi_mod

# SAVE: the work array rho * Omega**2 per shell.
dw2 = np.zeros(dimensions_mod.json + 1)


def getqua(hd: np.ndarray, hg: np.ndarray, hru: np.ndarray, omega: np.ndarray, m: int) -> None:
    """Discretize the second-order equation for the quadrupole term phi
    (centred on the smaller neighbouring interval), with phi -> 0 at the
    centre and phi ~ r**-3 outside the star, and solve the tridiagonal
    system by the Thomas algorithm; phi goes to quadru_mod.quad, g to
    quadru_mod.gg.

    Translates: getqua.f#GETQUA (7be4133)

    Arguments: hd (log rho), hg (gravity), hru (radius, unlogged), omega,
    m (number of shells).
    Changed in place: quadru_mod.quad, gg, and the work arrays of
    tridi_mod.

    Returns nothing.

    """

    c = constants_mod
    t = tridi_mod

    # rho * Omega**2 per shell, and 4 pi G (Fortran lines 48–51).
    for i in range(1, m + 1):
        dw2[i] = math.exp(c.cln * hd[i]) * omega[i] ** 2
    c4pig = math.exp(c.cln * (c.c4pil + c.cgl))

    # The central row (Fortran lines 52–65).
    # Question (see the page): drho is not divided by r2 - r1, and the
    # right-hand side uses Omega_1**2 drho/dr.
    dr = min(hru[1], hru[2] - hru[1])
    fpl = dr / (hru[2] - hru[1])
    fmi = dr / hru[1]
    fac = 1.0 / dr ** 2
    dri = 1.0 / dr / hru[1]
    t.a[1] = 0.0
    drho = (math.exp(c.cln * hd[2]) - math.exp(c.cln * hd[1])) / dr
    t.b[1] = (-fac * (fpl + fmi) + dri * (fmi - fpl) - 6.0 / hru[1] ** 2
              - c4pig * drho / hg[1])
    t.c[1] = fpl * (fac + dri)
    t.d[1] = -c4pig * drho * (hru[1] * omega[1]) ** 2 / 3.0 / hg[1]

    # The interior rows (Fortran lines 66–86).
    for i in range(2, m):
        drpl = hru[i + 1] - hru[i]
        drmi = hru[i] - hru[i - 1]
        dr = min(drmi, drpl)
        fpl = dr / drpl
        fmi = dr / drmi
        fac = 1.0 / dr ** 2
        dri = 1.0 / dr / hru[i]
        # Suspected bug (see the page): the last term should use hd[i - 1].
        # Question (see the page): the sign of the drho term in b.
        drho = 0.5 * (fpl * math.exp(c.cln * hd[i + 1])
                      + (fmi - fpl) * math.exp(c.cln * hd[i])
                      - fmi * math.exp(c.cln * hd[i])) / dr
        t.a[i] = fmi * (fac - dri)
        t.b[i] = (-fac * (fpl + fmi) + dri * (fmi - fpl) - 6.0 / hru[i] ** 2
                  - c4pig * drho / hg[i])
        t.c[i] = fpl * (fac + dri)
        drhow2 = 0.5 * (fpl * dw2[i + 1] + (fmi - fpl) * dw2[i]
                        - fmi * dw2[i - 1]) / dr
        t.d[i] = -drhow2 * c4pig * hru[i] ** 2 / 3.0 / hg[i]

    # The surface row: dphi/dr + 3 phi/r = 0 (Fortran lines 87–97).
    dr = hru[m] - hru[m - 1]
    rplu = hru[m] + dr
    fx = (hru[m] / rplu) ** 3
    fac = 1.0 / dr ** 2
    t.a[m] = fac
    t.b[m] = (fx - 2.0) * fac - 1.2e1 / hru[m] ** 2
    t.c[m] = 0.0
    t.d[m] = 0.0

    # The Thomas algorithm (from Numerical Recipes) (Fortran lines 98–110).
    if t.b[1] == 0.0:
        _zero_pivot(1, m)
    bet = t.b[1]
    t.u[1] = t.d[1] / bet
    for j in range(2, m + 1):
        t.gama[j] = t.c[j - 1] / bet
        bet = t.b[j] - t.a[j] * t.gama[j]
        if bet == 0.0:
            _zero_pivot(j, m)
        t.u[j] = (t.d[j] - t.a[j] * t.u[j - 1]) / bet

    # Back-substitution (Fortran lines 111–114).
    for j in range(m - 1, 0, -1):
        t.u[j] = t.u[j] - t.gama[j + 1] * t.u[j + 1]

    # The results to QUADRU_MOD (Fortran lines 115–118).
    for i in range(1, m + 1):
        quadru_mod.quad[i] = t.u[i]
        quadru_mod.gg[i] = hg[i]


def _zero_pivot(j, m):
    """Label 90: report the zero pivot at zone j to the terminal and the
    .short file, and stop (STOP 1, the run-time error code) (Fortran
    lines 122–128)."""

    text = (f" GETQUA: ZERO PIVOT IN THE TRIDIAGONAL SOLVE AT ZONE{j:5d}"
            f" OF{m:5d}\n RUN STOPPED")  # format 900
    print(text)
    iounits_mod.write(iounits_mod.ishort, text)
    raise SystemExit(1)
