"""Python version of getw.f (YREC development, 09c3b9b): the evolution of
the internal rotation over one model step."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import disk_mod
from pyrec import liov_mod
from pyrec import mdphy_mod
from pyrec import newrat_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldphy_mod
from pyrec import oldrat_mod
from pyrec import oldrot2_mod
from pyrec import oldrot_mod
from pyrec import ovrtrn_mod
from pyrec import quadru_mod
from pyrec import temp_mod
from pyrec.bursmix import bursmix
from pyrec.convec import convec
from pyrec.getrot import getrot
from pyrec.gettau import gettau
from pyrec.liburn import liburn
from pyrec.liburn2 import liburn2
from pyrec.lirate88 import lirate88
from pyrec.midmod import midmod
from pyrec.mixcz import mixcz
from pyrec.mwind import mwind
from pyrec.ovrot import ovrot
from pyrec.physic import physic
from pyrec.seculr import seculr
from pyrec.wczimp import wczimp

JSON = dimensions_mod.json

# GETW has SAVE: its local arrays are module-level. HCOMPS (the composition
# at the start of the step, filled by MIDMOD and read by SECULR) is owned
# by GETW.
hcomps = np.zeros((16, JSON + 1))
lcz = np.zeros(JSON + 1, dtype=bool)
hjmsav = np.zeros(JSON + 1)
eta2m = np.zeros(JSON + 1)
hdm = np.zeros(JSON + 1)
hgm = np.zeros(JSON + 1)
him = np.zeros(JSON + 1)
hlm = np.zeros(JSON + 1)
hpm = np.zeros(JSON + 1)
hrm = np.zeros(JSON + 1)
htm = np.zeros(JSON + 1)
lcm = np.zeros(JSON + 1, dtype=bool)
omegam = np.zeros(JSON + 1)
qiwm = np.zeros(JSON + 1)
r0m = np.zeros(JSON + 1)
lczm = np.zeros(JSON + 1, dtype=bool)
cod2 = np.zeros(JSON + 1)
hv = np.zeros(JSON + 1)
mxzone = np.zeros((13, 3), dtype=int)
mxzon0 = np.zeros((13, 3), dtype=int)
mrzone = np.zeros((14, 3), dtype=int)
nrzone = 0
nzone = 0


def _print_quad(bl, hl, hr, teffl, hstot, m, quad):
    """Label 9999: with LPRT0, print the surface quadrupole term.

    Returns BL, which is overwritten (Question 2 on the page).
    """
    c = constants_mod
    bl = math.log10(hl[m])
    rl = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
    fx = math.exp(c.cln * 3.0 * (hr[m] - rl))
    qsurf = fx * quad[m]
    psurf = math.exp(c.cln * (c.cgl + hstot - rl))
    print(f" QUAD {qsurf:12.3E} PHIS {psurf:12.3E} 3/2 QUAD/G "
          f"{-1.5 * qsurf / psurf:12.3E}")                 # format 9911
    return bl


def getw(bl, delts, dwmax, fp, ft, hcomp, hd, hjm, hl, hp, hr, hs,
         hs1, hs2, hstot, ht, lc, ljdot, m, smass,
         teffl, eta2, hg, hi, omega, qiw, r0, jenv0, cmixl, lprt0):
    """Evolve the angular momentum distribution over one model step.

    Enforces local conservation of angular momentum in radiative regions
    and the rotation law in convective regions; with instabilities
    (LINSTB), solves for the transport of angular momentum and the
    associated mixing in substeps, with the wind loss (LJDOT0) and, with
    LEXCOM, lithium and beryllium burning (Pinsonneault et al. 1989).

    Translates: getw.f#GETW (09c3b9b)

    Parameters
    ----------
    bl : float
        log10 L/Lsun; overwritten when ``lprt0`` is set.
    delts : float
        The model time step (s).
    dwmax : float
        Largest fractional change of Omega in the last step; 0 is replaced
        by DTWIND.
    fp, ft, hcomp, hd, hl, hp, hr, hs, hs1, hs2, ht, lc : arrays
        The model at the end of the step.
    hjm : array
        Specific angular momentum per shell; updated.
    hstot, smass, teffl : float
    ljdot : bool
        Wind on; set to LJDOT0 in the path without instabilities.
    m : int
    eta2, hg, hi, omega, qiw, r0 : arrays
        Rotation quantities; updated.
    jenv0 : int
        Base of the convective envelope at the last burning step.
    cmixl : float
        Mixing-length parameter, passed to PHYSIC and GETTAU.
    lprt0 : bool
        Print the surface quadrupole term.

    Returns
    -------
    bl, dwmax, ljdot, jenv0
    """
    global nrzone, nzone
    c = constants_mod
    nml = nmlpar_mod
    o2 = oldrot2_mod
    # Question (see the page): LBURS is never changed, so BURSMIX and
    # LIBURN2 are never used.
    lburs = False
    # The number of species tracked.
    nspec = 11
    if nml.lexcom:
        nspec = 15
    # Shells convective for angular momentum (convection plus overshoot).
    nrzone, nzone = ovrot(hcomp, hd, hp, hr, hs, ht, lc, m, lcz, mrzone,
                          mxzone, nrzone, nzone)
    # Disk locking.
    lock = False
    if nml.ldisk and disk_mod.sage <= nml.tdisk:
        lock = True

    if not nml.linstb:
        # Without instabilities.
        omegas = omega[m]
        if nml.lsolid or nml.impjmod == 1:
            lcz[1:m + 1] = True
        getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
        if not lock and nml.ljdot0 and lc[m]:
            hicz = 0.0
            if nml.lsolid or nml.impjmod == 1:
                imax = 1
            else:
                imax = mxzone[nzone, 1]
            for i in range(imax, m + 1):
                hicz = hicz + hi[i]
            hstop = math.exp(c.cln * hstot)
            if imax > 1:
                hsbot = 0.5 * (hs1[imax] + hs1[imax - 1])
            else:
                hsbot = 0.0
            ljdot = nml.ljdot0
            # Timestep-average loss rate.
            ovrtrn_mod.fracstep = 0.5
            mwind(bl, delts, hsbot, hstop, imax, m, ljdot, omegas,
                  smass, teffl, hicz, hjm)
            iend = m
            wczimp(hd, hjm, hr, hs, hs2, imax, iend, eta2, hi, omega,
                   qiw, r0, m)
        if lprt0:
            bl = _print_quad(bl, hl, hr, teffl, hstot, m, quadru_mod.quad)
        return bl, dwmax, ljdot, jenv0

    getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
    lskip = not nml.linstb and not ljdot
    if lskip or delts <= 0.0:
        if lprt0:
            bl = _print_quad(bl, hl, hr, teffl, hstot, m, quadru_mod.quad)
        return bl, dwmax, ljdot, jenv0

    # The number of diffusion substeps.
    if dwmax == 0.0:
        dwmax = nml.dtwind
    nstep = int(dwmax / nml.dtdif)
    if math.fmod(dwmax, nml.dtdif) != 0.0:
        nstep = nstep + 1
    nstep2 = int(nml.dtwind / nml.dtdif)
    if math.fmod(nml.dtwind, nml.dtdif) != 0.0:
        nstep2 = nstep2 + 1
    nstep = min(nstep, nstep2)
    dt = delts / float(nstep)

    # Physics and the overturn time, once per step.
    physic(fp, ft, hcomp, hd, hg, hl, hp, hr, hs, ht, lc, m, teffl, cmixl)
    # Question (see the page): ENVR is a scalar local standing for RBCZ.
    envr = 0.0
    envr = gettau(hcomp, hr, hp, hd, hs1, ht, fp, ft, teffl,
                  hstot, bl, m, lc, envr, cmixl)
    sumdt = 0.0
    lfirstit = True

    while True:                                          # label 30
        iredo = 0
        while True:                                      # label 40
            test = delts - delts * 1.0e-6
            if sumdt + dt >= test:
                dt = delts - sumdt
                sumdt = delts
            else:
                sumdt = sumdt + dt
            # Time change of theta, for the d theta/dt term of the circulation (MHP 06/02).
            if lfirstit:
                o2.thp[1] = o2.tho[1]
                o2.wmst[1] = oldrot_mod.wold[1]
                for i in range(2, m + 1):
                    o2.qwrmst[i] = o2.qwrst[i]
                    o2.thp[i] = o2.tho[i]
                    o2.wmst[i] = oldrot_mod.wold[i]
            else:
                o2.wmst[1] = omega[1]
                for i in range(2, m + 1):
                    wmid = 0.5 * (omega[i] + omega[i - 1])
                    dr = 10.0 ** hrm[i] - 10.0 ** hrm[i - 1]
                    qw = (omega[i] - omega[i - 1]) / dr
                    o2.thp[i] = o2.thm[i] * wmid * qw
                    o2.qwrmst[i] = qw
                    o2.wmst[i] = omega[i]
            fx = sumdt / delts
            ovrtrn_mod.fracstep = fx
            # The structure in the middle of the substep.
            (fx, hicz, hsbot, hstop, imin, imax, lallcz, lczsur,
             nrzone, nzone) = midmod(
                delts, dt, fx, hcomp, hd, hg, hjm, hl, hp, hr, hs, hs1,
                hs2, hstot, ht, lc, lfirstit, m, eta2m, hdm, hgm, 0.0, him,
                hlm, hpm, hrm, 0.0, 0.0, htm, 0, 0, False, lcm,
                lczm, False, omegam, r0m, qiwm, mrzone, mxzone, nrzone,
                nzone, hcomps)
            # The wind, here only for a fully convective star (or without
            # instabilities); otherwise SECULR calls MWIND.
            if not lock and ljdot and lczsur:
                omegas = omegam[m]
                if lallcz or not nml.linstb:
                    mwind(bl, dt, hsbot, hstop, imax, m, ljdot, omegas,
                          smass, teffl, hicz, hjm)
                    iend = m
                    wczimp(hdm, hjm, hrm, hs, hs2, imax, iend, eta2m, him,
                           omegam, qiwm, r0m, m)
            # Instabilities in radiative regions.
            if nml.linstb and not lallcz:
                limix = False       # set but no longer used (G Somers 6/14)
                mixcz(hcomp, hs2, lczm, m)
                for i in range(1, m + 1):                # label 60
                    hjmsav[i] = hjm[i]
                dt, lredo, iredo, hicz, omegas, lok = seculr(
                    dt, hdm, hgm, him, hlm, hpm, hrm, hs, hs1, hs2,
                    htm, hjmsav, lczm, m, omegam,
                    hjm, eta2m, qiwm, r0m, hcomp, imin, imax, bl,
                    hstot, smass, teffl, False, iredo,
                    hicz, hsbot, hstop, omegas, lczsur,
                    cod2, hv, False, hcomps)
                if lredo:
                    # Redo the substep (DT already halved): restore j, the
                    # composition and mu.
                    fx = 2.0 * dt / delts
                    sumdt = sumdt - 2.0 * dt
                    for i in range(1, m + 1):            # label 80
                        hjm[i] = hjmsav[i]
                        mdphy_mod.amum[i] = (mdphy_mod.amum[i]
                                             - fx * (temp_mod.hamu[i]
                                                     - oldphy_mod.amuo[i]))
                        for j in range(1, nspec + 1):    # label 70
                            hcomp[j, i] = oldmod_mod.hcompp[j, i]
                    continue                             # GOTO 40
            break

        # The composition at the start of the next substep.
        if not lburs:
            for i in range(1, m + 1):
                for j in range(1, 12):
                    oldmod_mod.hcompp[j, i] = hcomp[j, i]
        lfirstit = False

        # Light elements.
        if nml.lexcom:
            jcore1, jenv1, nrzone, nzone, nzone0 = convec(
                hcomp, hdm, hpm, hrm, hs, htm, lcm, m, mrzone, mxzone,
                mxzon0, 0, 0, nrzone, nzone, 0)
            # Overshoot depth alpha_e H_P, H_P = P r**2 / (rho G m).
            if nml.lovste and lcm[m]:
                liov_mod.pscal1 = nml.alphae * math.exp(
                    c.clndp * (hpm[jenv1] + 2.0 * hrm[jenv1] - hdm[jenv1]
                               - c.cgl - hs[jenv1]))
            else:
                liov_mod.pscal1 = 0.0
            lirate88(hcomp, hdm, htm, m, 2)
            if not nml.linstb or not lburs or lallcz:
                liburn(dt, hcomp, hrm, hs1, hs2, htm, jenv1, jenv0, m)
                jenv0 = jenv1
                liov_mod.pscal0 = liov_mod.pscal1
                for i in range(1, m + 1):                # label 155
                    oldrat_mod.rli60[i] = newrat_mod.rli6[i]
                    oldrat_mod.rli70[i] = newrat_mod.rli7[i]
                    oldrat_mod.rbe90[i] = newrat_mod.rbe9[i]
            else:
                liburn2(dt, hcomp, hrm, hs1, hs2, htm, jenv1, jenv0, m)
            limix = False
            mixcz(hcomp, hs2, lczm, m)
            # Zero out low abundances.
            for i in range(1, m + 1):                    # labels 100, 90
                for j in range(12, 16):
                    if hcomp[j, i] < 1.0e-24:
                        hcomp[j, i] = 0.0
            if not nml.linstb or lallcz:
                oldmod_mod.hcompp[12:16, 1:m + 1] = hcomp[12:16, 1:m + 1]
            elif not lburs:
                oldmod_mod.hcompp[12:16, 1:m + 1] = hcomp[12:16, 1:m + 1]

        # Bulirsch-Stoer mixing plus burning (never reached: LBURS false).
        if nml.linstb and not lallcz:
            if lburs:
                jenv0 = bursmix(cod2, dt, hcomp, hdm, hlm, hpm, hrm, hs, hs1,
                                hs2, hstot, htm, hv, imax, imin, jenv0, jenv1,
                                lok, lczm, m, mrzone, mxzone, nrzone, nzone)

        if not sumdt < delts:
            break                                        # else GOTO 30

    for i in range(1, m + 1):                            # label 110
        omega[i] = omegam[i]
    getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
    if lprt0:                                            # label 9999
        bl = _print_quad(bl, hl, hr, teffl, hstot, m, quadru_mod.quad)
    return bl, dwmax, ljdot, jenv0
