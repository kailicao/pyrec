"""Python version of grsett.f (YREC development, 09c3b9b): gravitational
settling of helium (and metals) over one time step."""

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import gravez_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.equal_to_model import equal_to_model
from pyrec.get_imp_diffco import get_imp_diffco
from pyrec.lax_wendrof1 import lax_wendrof1
from pyrec.lax_wendrof2 import lax_wendrof2
from pyrec.model_to_equal import model_to_equal
from pyrec.setup_grsett import setup_grsett
from pyrec.tridiag_gs import tridiag_gs

JSON = dimensions_mod.json

# GRSETT has SAVE: its local arrays are module-level. COD1Z..QCOD2Z were the
# COMMON block /GRAVSZ/ until rp_1; GRSETT now owns them, SETUP_GRSETT fills
# them and MODEL_TO_EQUAL reads them.
cod1z = np.zeros(JSON + 1)
cod2z = np.zeros(JSON + 1)
qcod1z = np.zeros(JSON + 1)
qcod2z = np.zeros(JSON + 1)
em = np.zeros(JSON + 1)
er = np.zeros(JSON + 1)
a = np.zeros(JSON + 1)
b = np.zeros(JSON + 1)
c = np.zeros(JSON + 1)
ex = np.zeros(JSON + 1)
hru = np.zeros(JSON + 1)
htu = np.zeros(JSON + 1)
em_h = np.zeros(JSON + 1)
ex_p = np.zeros(JSON + 1)
ex_orig = np.zeros(JSON + 1)
alpha = np.zeros(JSON + 1)
ex_h = np.zeros(JSON + 1)
cod1 = np.zeros(JSON + 1)
cod2 = np.zeros(JSON + 1)
qcod1x = np.zeros(JSON + 1)
qcod2x = np.zeros(JSON + 1)
ecod1 = np.zeros(JSON + 1)
ecod1_h = np.zeros(JSON + 1)
ex_prime = np.zeros(JSON + 1)
ecod2_h = np.zeros(JSON + 1)
eqcod1x_h = np.zeros(JSON + 1)
eqcod2x_h = np.zeros(JSON + 1)
ez_orig = np.zeros(JSON + 1)
ez_p = np.zeros(JSON + 1)


def grsett(dt, hcomp, hqpr, hr, hd, hs1, ht, lc, m, stot):
    """Settling of helium (and, with LDIFZ, of the metals) in one step.

    Regrids to equally spaced radii, treats the settling term explicitly
    (two-step Lax-Wendroff) and the diffusion term implicitly (an iterated
    tridiagonal solve, for X only), and regrids the change back.

    Translates: grsett.f#GRSETT (09c3b9b)

    Parameters
    ----------
    dt : float
        Time step (s). SETUP_GRSETT rescales it to the Bahcall-Loeb units
        and EQUAL_TO_MODEL restores it, so it is passed through both.
    hcomp : array
        Mass fractions; X and Y updated.
    hqpr, hr, hd, hs1, ht, lc : arrays
        d ln P / dr, log r, log rho, m (g), log T, convective flags. HQPR
        and HS1 are rescaled in place by SETUP_GRSETT.
    m : int
    stot : float
        Total mass; rescaled and restored like ``dt``.

    Returns
    -------
    dt, stot
    """
    gz = gravez_mod
    nml = nmlpar_mod
    # Coefficients, convection-zone edges, and the fully convective case.
    dt, stot, ibegin, iend, lallcz = setup_grsett(
        dt, hqpr, hr, hd, hs1, ht, lc, m, stot, cod1, cod2,
        hcomp, hru, htu, 0, 0, False, qcod1x, qcod2x, cod1z, cod2z,
        qcod1z, qcod2z)
    if lallcz:
        return dt, stot
    # To an equally spaced grid in radius.
    dr, npt = model_to_equal(cod1, cod2, hcomp, hru, hs1,
                             ibegin, iend, m, qcod1x, qcod2x, 0.0,
                             ecod1, ecod1_h, ecod2_h, em, em_h,
                             eqcod1x_h, eqcod2x_h, er, ex, ex_h, 0,
                             cod1z, cod2z, qcod1z, qcod2z)
    for i in range(1, npt + 1):                          # label 5
        ex_orig[i] = ex[i]
    if nml.ldifz:
        for i in range(1, npt + 1):
            ez_orig[i] = gz.ez[i]

    # Explicit term: half step, coefficients linearized, full step.
    ldoli = False
    lax_wendrof1(dt, ecod1, em, npt, stot, ex_h, ldoli)
    for i in range(1, npt):                              # label 10
        ecod1_h[i] = ecod1_h[i] + ex_h[i] * eqcod1x_h[i]
        ex_p[i] = ex_orig[i]
    ex_p[npt] = ex_orig[npt]
    if nml.ldifz:
        for i in range(1, npt):
            gz.ecod1z_h[i] = gz.ecod1z_h[i] + gz.ez_h[i] * gz.eqcod1z_h[i]
            ez_p[i] = ez_orig[i]
        ez_p[npt] = ez_orig[npt]
    lax_wendrof2(dt, ecod1_h, em_h, ex, npt, stot, ldoli)
    for i in range(1, npt + 1):                          # label 15
        ex_prime[i] = ex[i]

    # Implicit term.
    fac = constants_mod.c4pi * dt / dr
    alpha[1] = fac / em_h[1]
    for i in range(2, npt):                              # label 20
        alpha[i] = fac / (em_h[i] - em_h[i - 1])
    alpha[npt] = fac / (stot - em_h[npt - 1])
    converged = False
    for iter_ in range(1, nml.niter_gs + 1):             # label 100
        # Question (see the page): EX_H(I) is the change between points
        # I-1 and I, but GET_IMP_DIFFCO applies it to ECOD2_H(I), between
        # I and I+1; EX_H(1) is left from LAX_WENDROF1.
        for i in range(2, npt + 1):                      # label 30
            ex_h[i] = 0.5 * (ex[i] + ex[i - 1] - ex_p[i] - ex_p[i - 1])
        for i in range(1, npt + 1):                      # label 40
            ex_p[i] = ex[i]
        get_imp_diffco(alpha, ecod2_h, ex_h, eqcod2x_h, a, b, c, npt)
        tridiag_gs(a, b, c, ex_prime, npt, ex)
        # Question (see the page): the largest signed change, not ABS(DX).
        dxmax = 0.0
        imax = 0
        for i in range(1, npt + 1):                      # label 50
            dx = ex[i] - ex_p[i]
            if dx > dxmax:
                dxmax = dx
                imax = i
        iounits_mod.write(iounits_mod.ishort,
                          f" ITERATION {iter_:3d} DXMAX {dxmax:10.2E} "
                          f"IMAX {imax:4d}")             # format 90
        if dxmax < nml.grtol:
            converged = True
            break                                        # GOTO 120
    if not converged:
        # Question (see the page): failure is only reported.
        text = (f" GRSETT FAILED TO CONVERGE TO WITHIN {nml.grtol:9.3E} IN "
                f"{nml.niter_gs:3d}ITERATIONS\n LAST ITERATION CHANGE IN X "
                f"{dxmax:9.3E} IN EQUALLY SPACED SHELL {imax:5d}")  # format 110
        iounits_mod.write(iounits_mod.iowr, text)
        iounits_mod.write(iounits_mod.ishort, text)

    # label 120: return the change, not the abundance.
    for i in range(1, npt + 1):                          # label 130
        ex[i] = ex[i] - ex_orig[i]
    if nml.ldifz:
        # Question (see the page): ZZ1 and ZZ2 are not used.
        zz1 = gz.ez[1] - ez_orig[1]
        zz2 = gz.ez[npt] - ez_orig[npt]
        for i in range(1, npt + 1):
            gz.ez[i] = gz.ez[i] - ez_orig[i]
    # Back to the model's points (and the units restored).
    dt, stot = equal_to_model(dt, er, ex, ibegin, iend, npt, hcomp, hqpr,
                              hru, hs1, htu, m, stot)
    return dt, stot
