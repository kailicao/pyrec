"""Python version of mix.f (YREC 09c3b9b): the composition over one time
step (burning, mixing, semiconvection, settling, deuterium)."""

import math

import numpy as np

from pyrec import burn_mod
from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.convec import convec
from pyrec.dburn import dburn
from pyrec.deutrate import deutrate
from pyrec.eqburn import eqburn
from pyrec.grsett import grsett
from pyrec.kemcom import kemcom
from pyrec.microdiff import microdiff
from pyrec.rates import rates
from pyrec.sconvec import sconvec

JSON = dimensions_mod.json

# The Fortran's local arrays (SAVE): the burning rates, d ln P/dr, the zone
# lists and the convection flags used for mixing (LCC).
hr1 = np.zeros(JSON + 1)
hr2 = np.zeros(JSON + 1)
hr3 = np.zeros(JSON + 1)
hr4 = np.zeros(JSON + 1)
hr5 = np.zeros(JSON + 1)
hr6 = np.zeros(JSON + 1)
hr7 = np.zeros(JSON + 1)
hr8 = np.zeros(JSON + 1)
hr9 = np.zeros(JSON + 1)
hr10 = np.zeros(JSON + 1)
hr11 = np.zeros(JSON + 1)
hr12 = np.zeros(JSON + 1)
hr13 = np.zeros(JSON + 1)
hf1 = np.zeros(JSON + 1)
hf2 = np.zeros(JSON + 1)
hqpr = np.zeros(JSON + 1)
xsum = np.zeros(16)
mrzone = np.zeros((14, 3), dtype=int)
lcc = np.zeros(JSON + 1, dtype=bool)


def mix(delts, hcomp, hd, hl, hp, hr, hs, hs1, hs2, hstot,
        ht, itlvl, lc, m, ddage, jcore, jenv, mxzone, mxzon0, teffl, cmixl):
    """Update the composition hcomp over the step delts (s): find the
    convection zones, compute the burning rates, burn implicitly (KEMCOM)
    or, at itlvl = 1, explicitly (EQBURN), treat semiconvection, mix the
    convection zones, settle (MICRODIFF or GRSETT), clean up the abundances
    and burn deuterium.

    Translates: mix.f#MIX (09c3b9b)

    Arguments: the model (hd, hl, hp, hr, hs = log m, hs1 = m, hs2 = shell
    masses, log total mass hstot, ht, convection flags lc, m shells), the
    composition hcomp[k, j], the iteration level itlvl, teffl and the
    mixing length cmixl. Arrays have index 0 unused. Sets in place: hcomp,
    the zone lists mxzone and mxzon0, burn_mod.hcompm (the rates) and
    deuter_mod.drate; may set nmlpar_mod.ldify to False.

    Returns (ddage, jcore, jenv): the step in years and the edges of the
    central and surface convection zones.
    """
    c = constants_mod
    nml = nmlpar_mod
    drate = deuter_mod.drate
    # NSPEC: the number of species tracked.
    if nml.lexcom:
        nspec = 15
    else:
        nspec = 11
    # DDAGE: the step in years.
    ddage = delts / c.csecyr
    # Deep mixing: DPENV mixes the star to the mass fraction DPENV.
    # Question (see the page): from the centre, though the comment says from
    # the surface.
    if nml.dpenv < 1.0 and itlvl > 1:
        hsmax = nml.dpenv * math.exp(c.cln * hstot)
        for i in range(1, m + 1):
            if hs1[i] > hsmax:
                for j in range(i, m + 1):
                    lcc[j] = lc[j]
                break  # GO TO 5
            else:
                lcc[i] = True
        # A variant using a fraction of the peak luminosity is commented out in
        # the Fortran (MHP 1/95).
    else:
        for j in range(1, m + 1):
            lcc[j] = lc[j]
    # Convection zones with and without overshoot.
    jcore, jenv, nrzone, nzone, nzone0 = convec(
        hcomp, hd, hp, hr, hs, ht, lcc, m, mrzone, mxzone, mxzon0,
        jcore, jenv, 0, 0, 0)

    # Burning rates (HR1 - HR13, HF1, HF2), down to the cutoff TCUT(1).
    if nml.lmdot and nml.dmdt0 > 0.0:
        xh2test = max(hcomp[12, m], nml.compacc[12])
    else:
        xh2test = hcomp[12, m]
    i_cut = m + 1  # label 20: the first shell below the cutoff
    for i in range(1, m + 1):  # DO 10
        if ht[i] <= nml.tcut[1]:
            i_cut = i
            break
        dl = hd[i]
        tl = ht[i]
        x = hcomp[1, i]
        y = hcomp[2, i]
        xhe3 = hcomp[4, i]
        xc12 = hcomp[5, i]
        xc13 = hcomp[6, i]
        xn14 = hcomp[7, i]
        xo16 = hcomp[9, i]
        xo18 = hcomp[11, i]
        rates(dl, tl, x, y, xhe3, xc12, xc13, xn14, xo16, xo18, i,
              hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11,
              hr12, hr13, hf1, hf2)
        # The deuterium burning rate (MHP 5/02).
        if xh2test > 1.0e-11:
            deutrate(dl, tl, x, i, itlvl)
        else:
            drate[i] = 0.0
    for ii in range(i_cut, m + 1):  # DO 21: zero rates in cooler shells
        for arr in (hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10,
                    hr11, hr12, hr13, hf1, hf2):
            arr[ii] = 0.0
        drate[ii] = 0.0

    # Implicit burning: shell by shell in radiative zones (until the
    # cutoff), zone by zone in convection zones (which KEMCOM homogenizes).
    done = False
    for k in range(1, nrzone + 1):  # DO 40
        for j in range(mrzone[k, 1], mrzone[k, 2] + 1):
            if ht[j] <= nml.tcut[1]:
                done = True  # GOTO 45
                break
            ibegin = j
            iend = j
            kemcom(ht, ibegin, iend, hr1, hr2, hr3, hr4, hr5, hr6,
                   hr7, hr8, hr10, hr11, hr12, hf1, hs2, hcomp, ddage)
        if done:
            break
    for k in range(1, nzone + 1):  # DO 50
        ibegin = mxzone[k, 1]
        iend = mxzone[k, 2]
        kemcom(ht, ibegin, iend, hr1, hr2, hr3, hr4, hr5, hr6,
               hr7, hr8, hr10, hr11, hr12, hf1, hs2, hcomp, ddage)
    hcompm = burn_mod.hcompm
    for i in range(1, m + 1):
        for k, arr in enumerate((hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8,
                                 hr9, hr10, hr11, hr12, hr13, hf1, hf2),
                                start=1):
            hcompm[k, i] = arr[i]

    # ITLVL = 1: explicit H and He burning with equilibrium He3 and CN,
    # overwriting X (H burning) and Z, C12, O16 (He burning).
    # Question (see the page): DT in Gyr.
    if itlvl == 1:
        dt = ddage * 1.0e-9
        done = False
        for k in range(1, nrzone + 1):
            for j in range(mrzone[k, 1], mrzone[k, 2] + 1):
                if ht[j] <= nml.tcut[1]:
                    done = True  # GOTO 60
                    break
                ibegin = j
                iend = j
                (dcdt, dodt, dxdt, dydt, xc12, xo16, x, z) = eqburn(
                    hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr10, hr12, hs2, ht,
                    ibegin, iend, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
                if dxdt != 0.0:
                    hcomp[1, j] = x + dxdt * dt
                if dydt != 0.0:
                    hcomp[3, j] = z - dydt * dt
                    hcomp[5, j] = xc12 + dcdt * dt
                    hcomp[9, j] = xo16 + dodt * dt
            if done:
                break
        # label 60. Convection zones.
        for k in range(1, nzone + 1):
            ibegin = mxzone[k, 1]
            iend = mxzone[k, 2]
            (dcdt, dodt, dxdt, dydt, xc12, xo16, x, z) = eqburn(
                hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr10, hr12, hs2, ht,
                ibegin, iend, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            if dxdt != 0.0:
                x = x + dxdt * dt
                for i in range(ibegin, iend + 1):
                    hcomp[1, i] = x
            if dydt != 0.0:
                z = z - dydt * dt
                xc12 = xc12 + dcdt * dt
                xo16 = xo16 + dodt * dt
                for i in range(ibegin, iend + 1):
                    hcomp[3, i] = z
                    hcomp[5, i] = xc12
                    hcomp[9, i] = xo16

    # Semiconvection, if applicable.
    if nml.lsemic:
        if itlvl > 1:
            nzone = sconvec(delts, hcomp, hd, hl, hp, hr, hs, ht, m, mxzone,
                            nzone, teffl, cmixl)
    if nml.lsemic or itlvl == 1:
        # Mix each convection zone (shells I1..I2) to its mass-weighted mean.
        for j in range(1, nzone + 1):  # DO 100
            i1 = mxzone[j, 1]
            i2 = mxzone[j, 2]
            if i1 != 1 and i1 >= i2:
                continue
            wtsum = 0.0
            for i in range(1, nspec + 1):
                xsum[i] = 0.0
            for jj in range(i1, i2 + 1):
                wtsum = wtsum + hs2[jj]
                for i in range(1, nspec + 1):
                    xsum[i] = xsum[i] + hcomp[i, jj] * hs2[jj]
            for i in range(1, nspec + 1):
                xsum[i] = xsum[i] / wtsum
            for jj in range(i1, i2 + 1):
                for i in range(1, nspec + 1):
                    hcomp[i, jj] = xsum[i]
    # The mixed regions.
    if nzone >= 1:
        zones = "".join(f"({mxzone[i, 1]:5d},{mxzone[i, 2]:5d}) "
                        for i in range(1, nzone + 1))
        iounits_mod.write(iounits_mod.ishort,
                          " ZONES MIXED IN ORDER--" + zones)  # format 110

    # Microscopic diffusion (Bahcall and Loeb) (MHP 12/89).
    if nml.ldify:
        _settle(delts, hcomp, hd, hp, hr, hs, hs1, hstot, ht, lc, m, jenv)

    # Clean up small anomalies (negative abundances ...).
    for i in range(1, m + 1):  # DO 180
        for j in range(1, nspec + 1):
            hcomp[j, i] = max(hcomp[j, i], 0.0)
            hcomp[j, i] = min(hcomp[j, i], 1.0)
        hcomp[3, i] = min(hcomp[3, i], 1.0 - hcomp[1, i] - hcomp[4, i])
        hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]
    # Reset JENV, JCORE for deep mixing, with the true flags (MHP 1/95).
    if nml.dpenv < 1.0 and itlvl > 1:
        jcore, jenv, nrzone, nzone, nzone0 = convec(
            hcomp, hd, hp, hr, hs, ht, lc, m, mrzone, mxzone, mxzon0,
            jcore, jenv, nrzone, nzone, nzone0)
    # Deuterium burning (MHP 5/02).
    if nml.lexcom:
        dt = ddage * 1.0e-9
        done = False
        for k in range(1, nrzone + 1):
            for j in range(mrzone[k, 1], mrzone[k, 2] + 1):
                if ht[j] <= nml.tcut[1]:
                    done = True  # GOTO 190
                    break
                ibegin = j
                iend = j
                dburn(ibegin, iend, m, hs2, hcomp, dt)
            if done:
                break
        for k in range(1, nzone + 1):
            ibegin = mxzone[k, 1]
            iend = mxzone[k, 2]
            dburn(ibegin, iend, m, hs2, hcomp, dt)
    return ddage, jcore, jenv


def _settle(delts, hcomp, hd, hp, hr, hs, hs1, hstot, ht, lc, m, jenv):
    # The IF(LDIFY) block, labels 130-170: gravitational settling in
    # NSTEP substeps limited by the helium settling time at the outer
    # boundary. Reads the module-level LCC.
    c = constants_mod
    nml = nmlpar_mod
    # Question (see the page): settling stops for good below XMIN.
    if hcomp[1, 1] < nml.xmin:
        nml.ldify = False
        return  # GOTO 170
    stot = math.exp(c.cln * hstot)
    for i in range(1, m + 1):  # DO 130: d ln P/dr = -G m rho / (r^2 P)
        hqpr[i] = -math.exp(c.cln * (hd[i] + c.cgl + hs[i] - 2.0 * hr[i]
                                     - hp[i]))
    # The outer boundary: the surface convection zone, or the first point
    # (inward from JENV) where Y exceeds YMIN.
    if not lc[m]:
        # A warning only (DBG 2/92).
        iounits_mod.write(iounits_mod.ishort,
                          " NO SURFACE CZ - DIFFUSION NOT MEANINGFUL")
        return
    for j in range(jenv, 0, -1):  # DO 140
        if hcomp[2, j] > nml.ymin:
            break  # GOTO 150
    else:
        return  # Y < YMIN for the whole star
    # FM: the mass fraction above the outer point; TSCALE: the helium
    # settling time there (Michaud et al. 1984).
    fm = (stot - hs1[j]) / stot
    tscale = 4.348e21 * c.csecyr * fm / math.exp(c.cln * 1.5 * ht[j])
    dtmax = nml.dt_gs * tscale
    # Question (see the page): MOD's arguments swapped; NSTEP can be 0.
    nstep = int(delts / dtmax)
    if math.fmod(dtmax, delts) != 0.0:
        nstep = nstep + 1
    dt = delts / float(nstep)
    # Question (see the page): LCC, not the mixed-zone flags.
    for ns in range(1, nstep + 1):  # DO 160
        if nml.lnewdif:
            dt, stot = microdiff(
                dt, hcomp, hqpr, hr, hd, hs1, ht, lcc, m, stot)
        else:
            dt, stot = grsett(
                dt, hcomp, hqpr, hr, hd, hs1, ht, lcc, m, stot)
