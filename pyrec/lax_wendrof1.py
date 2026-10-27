"""Python version of lax_wendrof1.f (YREC development, 09c3b9b)."""

from pyrec import constants_mod
from pyrec import dimensions_mod  # noqa: F401  (JSON, the declared size)
from pyrec import gravez_mod
from pyrec import nmlpar_mod


def lax_wendrof1(dt, ecod1, em, npt, stot, ex_h, ldoli):
    """First (half) step of the Lax-Wendroff scheme for settling.

    The change in X at the midpoints of the equal grid over half a time
    step, from the difference of the settling fluxes, with zero flux at
    the centre and the surface; the same for Z (``gravez_mod.ez_h`` from
    ``gravez_mod.ecod1z``) when metal diffusion is on and ``ldoli`` is
    false. LAX_WENDROF2 follows.

    Translates: lax_wendrof1.f#LAX_WENDROF1 (09c3b9b)

    Parameters
    ----------
    dt : float
        Time step (s).
    ecod1 : ndarray
        The settling coefficients (fluxes) on the equal grid, 1..npt.
    em : ndarray
        The mass at each grid point, 1..npt.
    npt : int
        Number of grid points.
    stot : float
        Total mass.
    ex_h : ndarray
        Set in place: the half-step changes in X, ``ex_h[1..npt-1]``.
    ldoli : bool
        True skips the metals.

    Sets ``gravez_mod.ez_h[1..npt-1]`` with metal diffusion. Nothing is
    returned.
    """
    ecod1z = gravez_mod.ecod1z
    ez_h = gravez_mod.ez_h
    dt2 = 0.5 * dt * constants_mod.c4pi
    # Central boundary condition (F_1 = 0, m_1 = 0).
    emass = em[2]
    dex_h = dt2 * ecod1[2] / emass
    ex_h[1] = dex_h
    # General case (loop 10).
    for i in range(2, npt - 1):
        emass = em[i + 1] - em[i]
        dex_h = dt2 * (ecod1[i + 1] - ecod1[i]) / emass
        ex_h[i] = dex_h
    # Surface boundary condition (F_N = 0, m_N = M).
    # Question (see the page): this mass spans M - m_{N-1}.
    emass = stot - em[npt - 1]
    dex_h = -dt2 * ecod1[npt - 1] / emass
    ex_h[npt - 1] = dex_h
    # Metal diffusion.
    if nmlpar_mod.ldifz and not ldoli:
        emass = em[2]
        dez_h = dt2 * ecod1z[2] / emass
        ez_h[1] = dez_h
        for i in range(2, npt - 1):
            emass = em[i + 1] - em[i]
            dez_h = dt2 * (ecod1z[i + 1] - ecod1z[i]) / emass
            ez_h[i] = dez_h
        emass = stot - em[npt - 1]
        dez_h = -dt2 * ecod1z[npt - 1] / emass
        ez_h[npt - 1] = dez_h
