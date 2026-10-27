"""Python version of lax_wendrof2.f (YREC development, 09c3b9b)."""

from pyrec import constants_mod
from pyrec import dimensions_mod  # noqa: F401  (JSON, the declared size)
from pyrec import gravez_mod
from pyrec import nmlpar_mod


def lax_wendrof2(dt, ecod1_h, em_h, ex, npt, stot, ldoli):
    """Second (full) step of the Lax-Wendroff scheme for settling.

    Adds to X at the grid points the change over a full time step, from
    the difference of the settling fluxes at the midpoints, with zero flux
    at the centre and the surface; the same for Z (``gravez_mod.ez`` from
    ``gravez_mod.ecod1z_h``) when metal diffusion is on and ``ldoli`` is
    false. Follows LAX_WENDROF1.

    Translates: lax_wendrof2.f#LAX_WENDROF2 (09c3b9b)

    Parameters
    ----------
    dt : float
        Time step (s).
    ecod1_h : ndarray
        The settling coefficients at the midpoints, 1..npt-1.
    em_h : ndarray
        The mass at the midpoints, 1..npt-1.
    ex : ndarray
        Updated in place: X at the grid points, 1..npt.
    npt : int
        Number of grid points.
    stot : float
        Total mass.
    ldoli : bool
        True skips the metals.

    Updates ``gravez_mod.ez[1..npt]`` with metal diffusion. Nothing is
    returned.
    """
    ecod1z_h = gravez_mod.ecod1z_h
    ez = gravez_mod.ez
    dt2 = dt * constants_mod.c4pi
    # Central boundary condition.
    emass = em_h[1]
    dex = dt2 * ecod1_h[1] / emass
    ex[1] = ex[1] + dex
    # General case (loop 10).
    for i in range(2, npt):
        emass = em_h[i] - em_h[i - 1]
        dex = dt2 * (ecod1_h[i] - ecod1_h[i - 1]) / emass
        ex[i] = ex[i] + dex
    # Surface boundary condition.
    emass = stot - em_h[npt - 1]
    dex = -dt2 * ecod1_h[npt - 1] / emass
    ex[npt] = ex[npt] + dex
    # Metal diffusion (MHP 3/94).
    if nmlpar_mod.ldifz and not ldoli:
        emass = em_h[1]
        dez = dt2 * ecod1z_h[1] / emass
        ez[1] = ez[1] + dez
        for i in range(2, npt):
            emass = em_h[i] - em_h[i - 1]
            dez = dt2 * (ecod1z_h[i] - ecod1z_h[i - 1]) / emass
            ez[i] = ez[i] + dez
        emass = stot - em_h[npt - 1]
        dez = -dt2 * ecod1z_h[npt - 1] / emass
        ez[npt] = ez[npt] + dez
