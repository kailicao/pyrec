"""Python version of massloss.f (YREC 09c3b9b): the driver for mass loss or
accretion in a step."""

import math

from pyrec import atmprt_mod
from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import masschg2_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.alsurfp import alsurfp
from pyrec.eqstat import eqstat
from pyrec.mdot import mdot

# MHP 5/02: efficiency factor for the thermal energy content of accreted
# matter (DATA FACC/1.0D0/).
FACC = 1.0

# SAVEd locals passed to MDOT that are set only for accretion: for a loss
# MDOT gets the values of an earlier call (undefined on the first; 0.0 here).
# Question (see the page): arguments with no effect.
EACC = 0.0
ETHAV = 0.0
RMU = 0.0


def massloss(bl, dage, delts, hcomp, hd, hjm, hp, hr, hs, hs1, hs2, hstot,
             ht, jenv, lnew, m, omega, smass, teffl, senvold, lnewfit):
    """Set the mass-change rate (constant, or a Reimers formula), compute
    the convective envelope's mass and, for accretion, the energetics and
    entropies of the envelope and of accreted matter, then call MDOT to
    change the model's mass. dmdt > 0 is accretion, < 0 a loss.
    Experimental, and valid for Allard atmospheres only (llp 06/15/2009).

    Translates: massloss.f#MASSLOSS (09c3b9b)

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
    # Initialize mass loss at the default rate.
    dmdt = nml.dmdt0
    if nml.lmdot:
        ldomdot = True
    else:
        ldomdot = False
        lnewfit = False
        return delts, hstot, lnew, smass, senvold, lnewfit  # label 9999
    # Global quantities in cgs: radius, mass, age.
    rl = 0.5 * (bl + c.clsunl - c.c4pil - c.csigl - 4.0 * teffl)
    rtot = 10.0**rl
    stot = smass * c.cmsun
    agetot = dage * 1.0e9 * c.csecyr  # unused
    # Reimers formula, MDOT = K*L/g/R (overwrites the constant rate).
    if ldomdot and nml.lreimer:
        gsurf = 10.0**c.cgl * stot / rtot**2
        # Question (see the page): cgs, but MDOT reads Msun/yr.
        dmdt = nml.creim * 10.0**(bl + c.clsunl) / gsurf / rtot
    # Convective overturn time (computed earlier, in MIXCZ's call).
    if jenv < m:
        # Question (see the page): a debug line.
        print(ovrtrn_mod.taucz / c.csecyr, rtot / c.crsun)
        deuter_mod.jcz = jenv
    else:
        ovrtrn_mod.taucz = 0.0
        deuter_mod.jcz = m
    # MHP 8/10: the mass of the convective envelope.
    sumdm = 0.0
    for j in range(jenv, m + 1):
        sumdm = sumdm + hs2[j]
    if dmdt > 0.0:
        lnewfit = True
        t = 10.0**ht[jenv]
        p = 10.0**hp[jenv]
        d = 10.0**hd[jenv]
        beta = 1.0 - (c.ca3 * t**4 / p)
        RMU = p * beta / (d * t)
        # Kinetic energy per gram of accreted matter (virial theorem).
        EACC = FACC * stot * (10.0**c.cgl) / rtot
        eacc0 = EACC
        # Mass-weighted means over the convective envelope.
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
        lprt = False
        gl = c.cgl + hstot - 2.0 * rl
        # Experimental: valid for Allard atmospheres only.
        # Question (see the page): lalfail is not checked.
        lalfail = alsurfp(teffl, gl, lprt, False)
        t = 10.0**atmprt_mod.at
        p = 10.0**atmprt_mod.ap
        tl = atmprt_mod.at
        pl = atmprt_mod.ap
        x = hcomp[1, m]
        z = hcomp[3, m]
        fxion = [0.0] * 4
        lderiv = False
        ksaha = 1
        latmo = True
        (t, p, dl, d, beta, betai, beta14, RMU, amu, emu, eta, qdt, qdp,
         qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
            tl, t, pl, p, 0.0, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            lderiv, latmo, ksaha)
        beta = 1.0 - (c.ca3 * t**4 / p)
        RMU = p * beta / (d * t)
        masschg2_mod.sacc = RMU * (1.5 * math.log(t) - math.log(d))
        # Alternate surface pressure and temperature, from Stahler 1988.
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
            # Question (see the page): SACC2 is not used, SCEN zeroed.
            sacc2 = RMU * (1.5 * math.log(t) - math.log(d))
            masschg2_mod.scen = 0.0
        else:
            masschg2_mod.scen = 0.0
    # Mass loss or accretion (the call matches MDOT's 23 arguments since
    # 09c3b9b).
    if ldomdot:
        (delts, hstot, lnew, RMU, smass, dmdt, senvold) = mdot(
            delts, hcomp, hd, hjm, hp, hr, hs, hs1, hs2, hstot, ht, jenv,
            lnew, m, omega, RMU, rtot, smass, dmdt, EACC, ETHAV, sumdm,
            senvold)
    return delts, hstot, lnew, smass, senvold, lnewfit  # label 9999
