"""Python version of mdot.f (YREC 7be4133): mass loss or accretion over a
step."""

import math
import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import disk_mod
from pyrec import masschg2_mod
from pyrec import nmlpar_mod

# DATA DMSUN,DMWSUN,DMWMAX/3*0.0D0/: the scaled-solar-wind parameters, set
# nowhere else (formerly COMMON/MASSCHG3/), so the branch below is inert.
DMSUN = 0.0
DMWSUN = 0.0
DMWMAX = 0.0


def mdot(delts: float, hcomp: np.ndarray, hd: np.ndarray, hjm: np.ndarray, hp: np.ndarray, hr: np.ndarray,
         hs: np.ndarray, hs1: np.ndarray, hs2: np.ndarray, hstot: float, ht: np.ndarray, jenv: int, lnew: bool, m: int,
         omega: np.ndarray, rmu: float, rtot: float, smass: float, dmdt: float, eacc: float, ethav: float, sumdm: float,
         senvold: float) -> tuple[float, float, bool, float, float, float, float]:
    """Apply a mass change (loss, or accretion if dmdt > 0) to the model
    over a step: limit the step, rescale the mass coordinates, total mass
    and shell masses, remove angular momentum with lost mass, and mix
    accreted matter into the surface convection zone.

    Translates: mdot.f#MDOT (7be4133)

    Arguments: the step delts (s); the model hcomp[species, shell], hd, hjm,
    hp, hr, hs, hs1, hs2, ht, omega (arrays, index 0 unused, changed in
    place); log total mass hstot; jenv, the base of the surface convection
    zone; m, the number of shells; rmu; rtot, the radius; smass (Msun);
    dmdt (Msun/yr); eacc, ethav (no effect); sumdm, the zone's mass.

    Returns (delts, hstot, lnew, rmu, smass, dmdt, senvold): the possibly
    shortened step, the new log total mass, lnew = True (recompute the
    surface boundary condition), rmu (recomputed for accretion), the new
    mass, dmdt (changed only by the inert solar-wind branch) and senvold,
    the old hs[m] - hstot.

    """

    c = constants_mod
    senvold = hs[m] - hstot  # Fortran line 76.

    # Scaled solar wind mass loss (inert: DMWSUN = 0) (MHP 8/10) (Fortran
    # lines 99–104).
    if nmlpar_mod.lsolwind and nmlpar_mod.lrot and DMWSUN > 0.0:
        f1 = (omega[m] / DMWSUN)**2
        f2 = (DMWMAX / DMWSUN)**2
        dmdt = DMSUN * min(f1, f2)
        print(omega[m], dmdt)

    # Convert from solar masses per year to g/s (Fortran lines 111–116).
    rate = abs(dmdt) * c.cmsun / c.csecyr
    dmczt = sumdm

    # Mass of the surface convection zone below the fitting point (Fortran
    # lines 118–123).
    # Suspected bug (see the page): a hard stop with code 911 only.
    if jenv == m:
        raise SystemExit(911)
    if jenv > 1:
        dmcz = hs1[m] - hs1[jenv]
    else:
        dmcz = hs1[m]

    # Maximum step: not too much mass from the zone (FCZDMDT) or from the
    # star (FTOTDMDT) (Fortran lines 128–145).
    # A zero rate sets no limit.
    if rate > 0.0:
        dttot = nmlpar_mod.ftotdmdt * smass * c.cmsun / rate
        dtcz = nmlpar_mod.fczdmdt * dmcz / rate
        dtmax = min(dttot, dtcz)
    else:
        dtmax = delts
    # At most half the mass beyond the fitting point.
    if rate > 0.0:
        dtenv = 0.5 * (10.0**hstot - 10.0**hs[m]) / rate
        dtmax = min(dtmax, dtenv)
    if dtmax < delts:
        print(f" TIMESTEP REDUCED FOR MASS LOSS - OLD,NEW"
              f"{delts:12.3E}{dtmax:12.3E}")  # format 10
        delts = dtmax

    # Turn accretion off once the disk is exhausted (MHP 8/10) (Fortran
    # lines 147–157).
    # Suspected bug (see the page): a negative step if sage > tdisk.
    if nmlpar_mod.ldisk and dmdt > 0.0:
        test = disk_mod.sage + 1.0e-9 * delts / c.csecyr
        if test > nmlpar_mod.tdisk:
            delts = (nmlpar_mod.tdisk - disk_mod.sage) * 1.0e9 * c.csecyr
            ldone = True
        else:
            ldone = False
    else:
        ldone = False

    # Mass change in g (Fortran line 159).
    delm = dmdt * c.cmsun * delts / c.csecyr

    # Accretion: mean thermal energy and structure changes of the zone
    # (Fortran lines 162–202).
    if delm > 0.0:
        # ethnew and fac are not used (see the page).
        ethnew = (eacc * delm + ethav * sumdm) / (delm + sumdm)
        fac = math.log10(ethnew / ethav)  # noqa: F841  (unused in the Fortran too)
        t = 10.0**ht[jenv]
        p = 10.0**hp[jenv]
        d = 10.0**hd[jenv]
        beta = 1.0 - (c.ca3 * t**4 / p)
        rmu = p * beta / (d * t)
        # SCEN = 0, so these changes vanish (unfinished code, see the page).
        # Question (see the page): SCEN alone, or SACC - SCEN?
        dlogen = (delm / dmcz) * masschg2_mod.scen / rmu
        dlnm = 0.0
        dlogr = (c.cc23 * dlogen - c.cc13 * dlnm) / c.cln
        masschg2_mod.dlogt = dlnm / c.cln - dlogr
        masschg2_mod.dlogp = 2.0 * dlnm / c.cln - 4.0 * dlogr
        if jenv == 1:
            for j in range(jenv, m + 1):
                hr[j] = hr[j] + dlogr
                hp[j] = hp[j] + masschg2_mod.dlogp
                ht[j] = ht[j] + masschg2_mod.dlogt
        else:
            hp[jenv] = hp[jenv] + masschg2_mod.dlogp
            ht[jenv] = ht[jenv] + masschg2_mod.dlogt
            rb = 10.0**hr[jenv]
            facr = 10.0**dlogr
            for j in range(jenv + 1, m + 1):
                rinit = 10.0**hr[j]
                rnew = rb + facr * (rinit - rb)
                hr[j] = math.log10(rnew)

    # Rotation and mass loss: the lost matter carries its local specific
    # angular momentum (a thermal wind) (Fortran lines 206–236).
    if nmlpar_mod.lrot and delm < 0.0:
        surfi = 2.0 / 3.0 * rtot**2
        walpcz = nmlpar_mod.walpcz
        if walpcz >= 0.0:
            # Solid-body zone.
            delj = omega[m] * surfi * delm
        elif walpcz <= -2.0:
            # Constant J/M.
            delj = hjm[m] * delm
        else:
            # Extrapolate with the law WCZ imposes, Omega = C*r**WALPCZ.
            wsurf = omega[m] * (rtot / 10.0**hr[m])**walpcz
            delj = wsurf * surfi * delm
        # Total starting angular momentum of the zone.
        cjtot = hjm[jenv] * hs2[jenv]
        for j in range(jenv + 1, m + 1):
            cjtot = cjtot + hjm[j] * hs2[j]
        cjnew = cjtot + delj
        print(f"{delm:13.4E}{delj:13.4E}{omega[m]:13.4E}{cjnew:13.4E}"
              f"{cjtot:13.4E}{dmczt:13.4E}")  # format 22
        # J/M is scaled below, once the new shell masses are known.

    # Remove mass in the convection zone: compress toward its base
    # (hs1: shell centres in g; hs: log10 of hs1) (Fortran lines 240–258).
    # A fully convective star (dmcz = hs1[m]) is compressed about the center.
    if delm < 0.0:
        fx = (dmcz + delm) / dmcz
        if jenv == 1:
            for i in range(1, m + 1):
                hs1[i] = fx * hs1[i]
                hs[i] = math.log10(hs1[i])
        else:
            for i in range(jenv + 1, m + 1):
                hs1[i] = hs1[jenv] + fx * (hs1[i] - hs1[jenv])
                hs[i] = math.log10(hs1[i])

    # New total mass in solar units (smass) and log g (hstot, stotal)
    # (Fortran lines 262–270).
    smass = smass + delm / c.cmsun
    masschg2_mod.smass0 = smass
    dmsun0 = delm / c.cmsun
    print(f"MASS LOSS APPLIED - NEW M,DEL M{smass:19.10E}{dmsun0:19.10E}")
    hstot1 = 10.0**hstot
    hstot2 = hstot1 + delm
    hstot = math.log10(hstot2)
    comp_mod.stotal = hstot

    # Shell masses hs2 (g) (Fortran lines 276–280).
    # With jenv = 1 the innermost shell runs from the center to the
    # midpoint of hs1[1] and hs1[2] (the hs1[0] = -hs1[1] convention).
    if jenv == 1:
        hs2[1] = 0.5 * (hs1[1] + hs1[2])
    for i in range(max(jenv, 2), m):
        hs2[i] = 0.5 * (hs1[i + 1] - hs1[i - 1])
    hs2[m] = hstot2 - 0.5 * (hs1[m] + hs1[m - 1])

    # Scale J/M on the new shell masses, so the zone holds exactly cjnew
    # (Fortran lines 282–291).
    if nmlpar_mod.lrot and delm < 0.0:
        cjsum = 0.0
        for j in range(jenv, m + 1):
            cjsum = cjsum + hjm[j] * hs2[j]
        ratio = cjnew / cjsum
        for j in range(jenv, m + 1):
            hjm[j] = hjm[j] * ratio

    # Fortran lines 293–295.
    comp_mod.senv = hs[m] - hstot
    lnew = True

    # Accretion: remix the zone with the accreted composition (Fortran
    # lines 300–321).
    if delm > 0.0:
        for i in range(1, 12):
            # The zone's own sum(X dm), over the shell masses before the
            # accreted delm was added to hs2[m].
            sumx = 0.0
            summ = 0.0
            for j in range(jenv, m + 1):
                dmj = hs2[j]
                if j == m:
                    dmj = dmj - delm
                sumx = sumx + hcomp[i, j] * dmj
                summ = summ + dmj
            abund = (sumx + nmlpar_mod.compacc[i] * delm) / (summ + delm)
            for j in range(jenv, m + 1):
                hcomp[i, j] = abund
        deuter_mod.fmassacc = delm

    # Fortran line 323.
    if ldone:
        nmlpar_mod.lmdot = False
    return delts, hstot, lnew, rmu, smass, dmdt, senvold
