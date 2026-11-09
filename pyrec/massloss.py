"""Python version of massloss.f (YREC 7be4133): the driver for mass loss or
accretion in a step."""

import math
import numpy as np

from pyrec import atmprt_mod
from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import masschg2_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.alsurfp import alsurfp
from pyrec.eqstat import eqstat
from pyrec.mdot import mdot

# Efficiency factor for the thermal energy content of accreted matter (DATA
# FACC/1.0D0/) (MHP 5/02).
FACC = 1.0

# SAVEd locals passed to MDOT that are set only for accretion: for a loss
# MDOT gets the values of an earlier call (undefined on the first; 0.0 here).
# MDOT does not use them (see the page).
EACC = 0.0
ETHAV = 0.0
RMU = 0.0


def massloss(bl: float, dage: float, delts: float, hcomp: np.ndarray, hd: np.ndarray, hjm: np.ndarray, hp: np.ndarray,
             hr: np.ndarray, hs: np.ndarray, hs1: np.ndarray, hs2: np.ndarray, hstot: float, ht: np.ndarray, jenv: int,
             lnew: bool, m: int, omega: np.ndarray, smass: float, teffl: float, senvold: float,
             lnewfit: bool) -> tuple[float, float, bool, float, float, bool]:
    """Set the mass-change rate (constant, or a Reimers formula), compute
    the convective envelope's mass and, for accretion, the energetics and
    entropies of the envelope and of accreted matter, then call MDOT to
    change the model's mass. dmdt > 0 is accretion, < 0 a loss.
    Experimental, and valid for Allard atmospheres only (llp 06/15/2009).

    Translates: massloss.f#MASSLOSS (7be4133)

    Arguments: bl (log L/Lsun), dage (Gyr), delts (s); the model hcomp,
    hd, hjm, hp, hr, hs, hs1, hs2, ht, omega (arrays, index 0 unused,
    changed in place by MDOT); hstot (log total mass, g); jenv (base of
    the surface convection zone); lnew; m (number of shells); smass
    (Msun); teffl (log Teff); senvold; lnewfit is output only.

    Returns (delts, hstot, lnew, smass, senvold, lnewfit): the first four
    and senvold as changed by MDOT (unchanged if LMDOT is false), and
    lnewfit, True only for accretion (the caller then restores the
    envelope fitting point).

    """

    global EACC, ETHAV, RMU
    c = constants_mod
    nml = nmlpar_mod

    # Initialize mass loss at the default rate (Fortran lines 80–87).
    dmdt = nml.dmdt0
    if nml.lmdot:
        ldomdot = True
    else:
        ldomdot = False
        lnewfit = False
        return delts, hstot, lnew, smass, senvold, lnewfit  # label 9999

    # Global quantities in cgs: radius, mass, age (Fortran lines 92–97).
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    rtot = 10.0**rl
    stot = smass * c.cmsun
    agetot = dage * 1.0e9 * c.csecyr  # noqa: F841  (unused in the Fortran too)

    # Reimers formula, MDOT = K*L/g/R (overwrites the constant rate)
    # (Fortran lines 100–103).
    if ldomdot and nml.lreimer:
        gsurf = 10.0**c.cgl * stot / rtot**2
        # Suspected bug (see the page): cgs, but MDOT reads Msun/yr.
        dmdt = nml.creim * 10.0**(bl + c.clsunl) / gsurf / rtot

    # Convective overturn time (computed earlier, in MIXCZ's call) (Fortran
    # lines 106–122).
    if jenv < m:
        # Question (see the page): a debug line.
        print(f' MASSLOSS: TAUCZ (YR), R/RSUN {ovrtrn_mod.taucz / c.csecyr:12.4E}'
              f'{rtot / c.crsun:12.4E}')
        deuter_mod.jcz = jenv
    else:
        ovrtrn_mod.taucz = 0.0
        deuter_mod.jcz = m

    # The mass of the convective envelope (MHP 8/10) (Fortran lines
    # 126–129).
    sumdm = 0.0
    for j in range(jenv, m + 1):
        sumdm = sumdm + hs2[j]

    # Accretion (Fortran lines 131–228).
    if dmdt > 0.0:
        lnewfit = True
        t = 10.0**ht[jenv]
        p = 10.0**hp[jenv]
        d = 10.0**hd[jenv]
        beta = 1.0 - (c.ca3 * t**4 / p)
        RMU = p * beta / (d * t)

        # Kinetic energy per gram of accreted matter (virial theorem)
        # (Fortran lines 145–146).
        EACC = FACC * stot * (10.0**c.cgl) / rtot
        eacc0 = EACC

        # Mass-weighted means over the convective envelope (Fortran lines
        # 149–177).
        sumeth = 0.0
        ethaccbar = 0.0
        masschg2_mod.scen = 0.0
        for j in range(jenv, m + 1):
            t = 10.0**ht[j]
            p = 10.0**hp[j]
            d = 10.0**hd[j]
            r = 10.0**hr[j]
            delacc = 0.5 * EACC * (rtot / r - 1.0)
            ethacc = EACC + delacc
            ethaccbar = ethaccbar + hs2[j] * ethacc
            beta = 1.0 - (c.ca3 * t**4 / p)
            RMU = p * beta / (d * t)
            slocal = RMU * (1.5 * math.log(t) - math.log(d))
            masschg2_mod.scen = masschg2_mod.scen + slocal * hs2[j]
            eth = p * beta / d
            sumeth = sumeth + eth * hs2[j]
        ETHAV = sumeth / sumdm
        EACC = ethaccbar / sumdm
        masschg2_mod.scen = masschg2_mod.scen / sumdm

        # The photosphere from the Allard tables, and its entropy (Fortran
        # lines 178–205).
        lprt = False
        gl = c.cgl + hstot - 2.0 * rl
        # Experimental: valid for Allard atmospheres only.
        # Question (see the page): Allard atmospheres only.
        lalfail = alsurfp(teffl, gl, lprt, False)
        x = hcomp[1, m]
        z = hcomp[3, m]
        fxion = [0.0] * 4
        lderiv = False
        ksaha = 1
        latmo = True
        # Use the Allard surface point only if the lookup succeeded; AT and
        # AP otherwise still hold the previous lookup's values.
        if not lalfail:
            t = 10.0**atmprt_mod.at
            p = 10.0**atmprt_mod.ap
            tl = atmprt_mod.at
            pl = atmprt_mod.ap
            (t, p, dl, d, beta, betai, beta14, RMU, amu, emu, eta, qdt, qdp,
             qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
                tl, t, pl, p, 0.0, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                lderiv, latmo, ksaha)
            beta = 1.0 - (c.ca3 * t**4 / p)
            RMU = p * beta / (d * t)
            masschg2_mod.sacc = RMU * (1.5 * math.log(t) - math.log(d))
        else:
            print('MASSLOSS: Allard surface lookup failed; '
                  'accreted-matter entropy not computed')
            # EQSTAT's outputs not yet set in this call, passed in below
            # (any value; Fortran keeps SAVEd values from an earlier call).
            dl = betai = beta14 = amu = emu = eta = 0.0
            qdt = qdp = qcp = dela = qdtt = qdtp = qat = qap = 0.0
            qcpt = qcpp = 0.0

        # Alternate surface pressure and temperature, from Stahler 1988
        # (Fortran lines 209–227).
        if FACC > 0.0:
            dmdt2 = dmdt * c.cmsun / c.csecyr
            p = dmdt2 / c.c4pi * math.sqrt(2.0 * FACC * 10.0**gl / rtot**3)
            t = (dmdt2 * eacc0 * 0.75 / c.c4pi / rtot**2 / c.csig)**0.25
            pl = math.log10(p)
            tl = math.log10(t)
            (t, p, dl, d, beta, betai, beta14, RMU, amu, emu, eta, qdt, qdp,
             qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
                tl, t, pl, p, dl, d, x, z, beta, betai, beta14, fxion, RMU,
                amu, emu, eta, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
                qcpt, qcpp, lderiv, latmo, ksaha)
            beta = 1.0 - (c.ca3 * t**4 / p)
            RMU = p * beta / (d * t)
            # SACC2 is not used and SCEN is zeroed: unfinished code (see the
            # page).
            sacc2 = RMU * (1.5 * math.log(t) - math.log(d))  # noqa: F841  (unused in the Fortran too)
            masschg2_mod.scen = 0.0
        else:
            masschg2_mod.scen = 0.0

    # Mass loss or accretion (MDOT, with its 23 arguments) (Fortran line
    # 240).
    if ldomdot:
        (delts, hstot, lnew, RMU, smass, dmdt, senvold) = mdot(
            delts, hcomp, hd, hjm, hp, hr, hs, hs1, hs2, hstot, ht, jenv,
            lnew, m, omega, RMU, rtot, smass, dmdt, EACC, ETHAV, sumdm,
            senvold)
    return delts, hstot, lnew, smass, senvold, lnewfit  # label 9999
