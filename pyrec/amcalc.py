"""Python version of amcalc.f (YREC 7be4133): the structure factor of the
magnetic-braking torque (Matt et al. 2012)."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod


def amcalc(smass: float, bl: float, teffl: float, structfactor: float) -> float:
    """Compute the structure factor of the wind torque: the mass, R/Rsun,
    L/Lsun, the photospheric pressure over the solar one and, if LROSSBY,
    the convective overturn time over the solar one, each raised to its
    exponent from nmlpar_mod.

    Translates: amcalc.f#AMCALC (7be4133)

    Arguments: smass (mass, solar units), bl (log L/Lsun), teffl
    (log Teff); structfactor is output only (the value passed is unused).

    Returns structfactor.

    """

    # Mass (Fortran line 70).
    rmass = smass

    # Radius, from L and Teff (Stefan-Boltzmann), in solar units (Fortran
    # lines 72–74).
    rl = 0.5 * (bl + constants_mod.clsunl - constants_mod.c4pil
                - constants_mod.csigl - 4.0 * teffl)
    rtot = math.exp(constants_mod.cln * rl)
    rrad = rtot / constants_mod.crsun

    # Luminosity (Fortran line 76).
    rlum = 10.0 ** bl

    # Photospheric pressure, interpolated within the step (Fortran line 78).
    rphot = (10.0 ** (ovrtrn_mod.pphot0 + ovrtrn_mod.fracstep
                      * (ovrtrn_mod.pphot - ovrtrn_mod.pphot0))
             / 10.0 ** nmlpar_mod.pmmsolp)

    # Convective overturn timescale (Fortran lines 80–84).
    if nmlpar_mod.lrossby:
        rtau = (ovrtrn_mod.taucz0 + ovrtrn_mod.fracstep
                * (ovrtrn_mod.taucz - ovrtrn_mod.taucz0)) / nmlpar_mod.pmmsoltau
    else:
        rtau = 1.0

    # Combine them all (Fortran lines 86–87).
    structfactor = (rmass ** nmlpar_mod.exm * rrad ** nmlpar_mod.exr
                    * rlum ** nmlpar_mod.exl * rphot ** nmlpar_mod.expr
                    * rtau ** nmlpar_mod.extau)
    return structfactor
