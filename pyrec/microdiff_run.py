"""Python version of microdiff_run.f (YREC 09c3b9b): one species diffused
over one step on the equally spaced grid, for MICRODIFF."""

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.get_imp_diffco import get_imp_diffco
from pyrec.lax_wendrof1 import lax_wendrof1
from pyrec.lax_wendrof2 import lax_wendrof2
from pyrec.microdiff_cod import microdiff_cod
from pyrec.tridiag_gs import tridiag_gs

JSON = dimensions_mod.json

# The Fortran's local arrays (SAVE).
ed = np.zeros(JSON + 1)
ed_h = np.zeros(JSON + 1)
ed_orig = np.zeros(JSON + 1)
ed_orig_h = np.zeros(JSON + 1)
ed_p = np.zeros(JSON + 1)
ed_prime = np.zeros(JSON + 1)
ecod1 = np.zeros(JSON + 1)
ecod1_h = np.zeros(JSON + 1)
ecod2 = np.zeros(JSON + 1)
ecod2_h = np.zeros(JSON + 1)
eqcod2x_h = np.zeros(JSON + 1)
alpha = np.zeros(JSON + 1)
a = np.zeros(JSON + 1)
b = np.zeros(JSON + 1)
c = np.zeros(JSON + 1)


def microdiff_run(dr, dt, stot, npt, em, er, eden, etem, ehqpr, edelr,
                  espec, edxdr, em_h, er_h, eden_h, etem_h, ehqpr_h,
                  edelr_h, espec_h, edxdr_h, thea, thez, j):
    """Diffuse species column j over one step on the equal grid: the
    advective term by the two-step Lax-Wendroff scheme, then the diffusive
    term implicitly (a tridiagonal solve, iterated up to NITER_GS times to
    the tolerance GRTOL). Replaces espec[j, 1..npt] by the change in the
    abundance.

    Translates: microdiff_run.f#MICRODIFF_RUN (09c3b9b)

    Arguments: the spacing dr; dt and stot in Bahcall and Loeb's units;
    npt grid points; the grid variables em, er, eden, etem, ehqpr, edelr
    and their midpoint versions (_h, npt - 1 of them); the abundances
    espec[k, i] and espec_h[k, i] (k = 1..3); the hydrogen gradient edxdr,
    edxdr_h (set by MICRODIFF_COD on the hydrogen call); the species'
    weight thea, charge thez and column j. Arrays have index 0 unused.

    Sets in place: espec[j, i] (the change), espec_h[j, i] (the half-step
    abundances), edxdr, edxdr_h. Returns nothing.
    """
    # Dummy vectors for the diffusing element, and its original run.
    for i in range(1, npt):
        ed[i] = espec[j, i]
        ed_h[i] = espec_h[j, i]
    ed[npt] = espec[j, npt]
    for i in range(1, npt):
        ed_orig[i] = ed[i]
        ed_orig_h[i] = ed_h[i]
    ed_orig[npt] = ed[npt]
    # LDOLI: the Lax-Wendroff routines do not use the old metal vectors.
    ldoli = True
    # Coefficients at the grid points.
    microdiff_cod(npt, espec, er, eden, etem, ehqpr, edelr, ecod1, ecod2,
                  edxdr, thea, thez, j)
    # Lax-Wendroff, first half step: the change at the midpoints.
    lax_wendrof1(dt, ecod1, em, npt, stot, ed_h, ldoli)
    for i in range(1, npt):
        espec_h[j, i] = ed_orig_h[i] + ed_h[i]
    # Coefficients at the npt - 1 midpoints, from the provisional solution.
    npt2 = npt - 1
    microdiff_cod(npt2, espec_h, er_h, eden_h, etem_h, ehqpr_h, edelr_h,
                  ecod1_h, ecod2_h, edxdr_h, thea, thez, j)
    # Lax-Wendroff, full step at the points.
    lax_wendrof2(dt, ecod1_h, em_h, ed, npt, stot, ldoli)
    # The implicit term: ed_p is the abundance at the start of the step,
    # ed_prime the abundance without the second term.
    for i in range(1, npt + 1):
        ed_p[i] = ed_orig[i]
        ed_prime[i] = ed[i]
    fac = constants_mod.c4pi * dt / dr
    alpha[1] = fac / em_h[1]
    for i in range(2, npt):
        alpha[i] = fac / (em_h[i] - em_h[i - 1])
    alpha[npt] = fac / (stot - em_h[npt - 1])
    ishort = iounits_mod.ishort
    converged = False
    dxmax = 0.0
    imax = 0
    for it in range(1, nmlpar_mod.niter_gs + 1):        # loop 100
        # Question (see the page): midpoint i is taken between points i-1
        # and i here, between i and i+1 in MICRODIFF_MTE.
        for i in range(2, npt + 1):
            ed_h[i] = 0.5 * (ed[i] + ed[i - 1] - ed_p[i] - ed_p[i - 1])
        # Question (see the page): eqcod2x_h is zeroed, so every iteration
        # solves the same system.
        for i in range(1, npt):
            ed_p[i] = ed[i]
            eqcod2x_h[i] = 0.0
        ed_p[npt] = ed[npt]
        get_imp_diffco(alpha, ecod2_h, ed_h, eqcod2x_h, a, b, c, npt)
        tridiag_gs(a, b, c, ed_prime, npt, ed)
        # Question (see the page): the largest signed change, not |change|.
        dxmax = 0.0
        imax = 0
        for i in range(1, npt + 1):
            dx = ed[i] - ed_p[i]
            if dx > dxmax:
                dxmax = dx
                imax = i
        iounits_mod.write(ishort, f" ITERATION {it:3d} DXMAX "
                          f"{dxmax:10.2E} IMAX {imax:4d}")     # format 90
        if dxmax < nmlpar_mod.grtol:
            converged = True
            break                                           # GOTO 120
    if not converged:
        text = (f" MICRODIFF FAILED TO CONVERGE TO WITHIN "
                f"{nmlpar_mod.grtol:9.3E} IN {nmlpar_mod.niter_gs:3d}"
                f"ITERATIONS\n LAST ITERATION CHANGE IN D {dxmax:9.3E}"
                f" IN EQUALLY SPACED SHELL {imax:5d}")     # format 110
        iounits_mod.write(iounits_mod.iowr, text)
        iounits_mod.write(ishort, text)
    # label 120: the change in the diffused element.
    for i in range(1, npt + 1):
        espec[j, i] = ed[i] - ed_orig[i]
