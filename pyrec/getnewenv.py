"""Python version of getnewenv.f (YREC development, 7be4133): move the
outer fitting point outward by turning part of the envelope into new
interior points."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import dimensions_mod
from pyrec import envstruct_mod
from pyrec import nmlpar_mod
from pyrec.envint import envint
from pyrec.getrot import getrot

JSON = dimensions_mod.json

# Locals kept between calls (SAVE). OMEGA is local, not an argument
# (Fortran lines 74–78).
omega = np.zeros(JSON + 1)
dum1 = np.zeros(5)
dum2 = np.zeros(4)
dum3 = np.zeros(4)
dum4 = np.zeros(4)


def getnewenv(senvold: float, hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray,
              hs: np.ndarray, hs1: np.ndarray, hs2: np.ndarray, hstot: float, ht: np.ndarray, lc: np.ndarray,
              eta2: np.ndarray, hi: np.ndarray, hjm: np.ndarray, qiw: np.ndarray, r0: np.ndarray, hkerot: np.ndarray,
              bl: float, sjtot: float, skerot: float, teffl: float, m: int, lnew: bool, cmixl: float) -> tuple[float,
              float, float, int, bool]:
    """Integrate the envelope from the surface, match it to the last
    interior point, add the envelope points below the desired fitting
    point (log m_fit/M = senvold) as new shells, and reset rotation for
    them (method on the page). Called by the main program when LNEWFIT
    is set.

    Translates: getnewenv.f#GETNEWENV (7be4133)

    Arguments: senvold (the desired log(m_fit/M), capped here at -1e-12);
    the model arrays hcomp, hd, hl, hp, hr, hs, hs1, hs2, ht, lc and the
    rotation arrays eta2, hi, hjm, qiw, r0, hkerot, all changed in place;
    hstot, bl, teffl; m (points in the model); cmixl (the mixing length,
    passed to envint); sjtot, skerot and lnew are only set.

    Returns (senvold, sjtot, skerot, m, lnew).

    """

    c = constants_mod
    nml = nmlpar_mod
    cmp = comp_mod
    env = envstruct_mod
    dh = debhu_mod

    # Fortran lines 82–86.
    if nml.lexcom:
        jend = 15
    else:
        jend = 11

    # The envelope mass is at least 1e-12 (Fortran lines 87–91).
    senvold = min(senvold, -1.0e-12)
    mm = m
    senv1 = cmp.senv

    # The envelope integration's step sizes, set to HPTTOL(8) for now
    # (Fortran lines 92–103).
    emax0 = nml.envmax
    emin0 = nml.envmin
    ebeg0 = nml.envbeg
    nml.envmax = nml.hpttol[8]
    nml.envmin = nml.hpttol[8]
    nml.envbeg = nml.hpttol[8]
    lsbc0 = False
    lprt = True
    katm = 0
    kenv = 0
    ksaha = 0

    # Surface luminosity, radius, gravity and composition, no rotational
    # distortion (Fortran lines 104–132).
    b = math.exp(c.cln * bl)
    rl = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
    gl = c.cgl + cmp.stotal - rl - rl
    x = hcomp[1, m]
    z = hcomp[3, m]
    # No rotational distortion.
    fpl = 1.0
    ftl = 1.0
    # Do not overshoot the pressure of the outer model point.
    plim = hp[m]
    lpulpt = False
    ixx = 0
    if nml.ldh:
        # The Debye-Hueckel composition, from the surface shell.
        dh.xxdh = hcomp[1, m]
        dh.yydh = hcomp[2, m] + hcomp[4, m]
        dh.zzdh = hcomp[3, m]
        dh.zdh[1] = hcomp[5, m] + hcomp[6, m]
        dh.zdh[2] = hcomp[7, m] + hcomp[8, m]
        dh.zdh[3] = hcomp[9, m] + hcomp[10, m] + hcomp[11, m]

    # Integrate down to the current fitting point; with spots, at the
    # ambient temperature (Fortran lines 133–146).
    idum = 0
    if lc[m] and nml.spotf != 0.0 and nml.spotx != 1.0:
        ateffl = teffl - 0.25 * math.log10(
            nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
    else:
        ateffl = teffl
    rl, idum, katm, kenv, ksaha = envint(
        b, fpl, ftl, gl, hstot, ixx, lprt, lsbc0, plim, rl, ateffl, x, z,
        dum1, idum, katm, kenv, ksaha, dum2, dum3, dum4, lpulpt, cmixl)

    # Restore the step sizes (Fortran lines 147–153).
    nml.envmax = emax0
    nml.envmin = emin0
    nml.envbeg = ebeg0
    cmp.senv = senvold
    if m + env.numenv >= JSON:
        raise SystemExit("STOP 9999")

    # Match the envelope to the last interior point, shift the points
    # down by one, and drop those closer than 1e-10 in log m (Fortran
    # lines 154–182).
    dp = hp[m] - env.envp[1]
    dd = hd[m] - env.envd[1]
    dt = ht[m] - env.envt[1]
    dr = hr[m] - env.envr[1]
    test = hs[m]
    numnew = 0
    for j in range(1, env.numenv):
        test2 = cmp.stotal + env.envs[j + 1]
        if test2 - test > 1.0e-10:
            env.envd[j] = env.envd[j + 1] + dd
            env.envp[j] = env.envp[j + 1] + dp
            env.envr[j] = env.envr[j + 1] + dr
            env.envs[j] = env.envs[j + 1]
            env.envt[j] = env.envt[j + 1] + dt
            env.envx[j] = env.envx[j + 1]
            env.envz[j] = env.envz[j + 1]
            env.lcenv[j] = env.lcenv[j + 1]  # MHP 10/30/26
            test = test2
            numnew = numnew + 1
    env.numenv = numnew

    # The new points (Fortran lines 183–258).
    senv = cmp.senv
    stotal = cmp.stotal
    for j in range(m + 1, m + env.numenv + 1):
        jm = j - m
        # Luminosity held constant.
        hl[j] = hl[m]
        if env.envs[jm] < senv:
            hd[j] = env.envd[jm]
            hp[j] = env.envp[jm]
            hr[j] = env.envr[jm]
            hs[j] = env.envs[jm] + stotal
            ht[j] = env.envt[jm]
            hcomp[1, j] = env.envx[jm]
            hcomp[3, j] = env.envz[jm]
            for k in range(4, jend + 1):
                hcomp[k, j] = hcomp[k, m]
            hcomp[2, j] = 1.0 - hcomp[1, j] - hcomp[3, j] - hcomp[4, j]
            lc[j] = env.lcenv[jm]
        else:
            # Above the new fitting point: interpolate linearly in log m
            # to it, set the number of points and leave the loop.
            if jm == 1:
                # Between the last interior and the first envelope point.
                x0 = hs[m]
                x1 = stotal + senv
                x2 = env.envs[jm] + stotal
                if x2 - x0 < 1.0e-14:
                    raise SystemExit("STOP 9998")
                fx = (x1 - x0) / (x2 - x0)
                hd[j] = hd[m] + fx * (env.envd[jm] - hd[m])
                hp[j] = hp[m] + fx * (env.envp[jm] - hp[m])
                hr[j] = hr[m] + fx * (env.envr[jm] - hr[m])
                hs[j] = x1
                ht[j] = ht[m] + fx * (env.envt[jm] - ht[m])
                # Suspected bug (see the page): the sign, a + f(a - b).
                hcomp[1, j] = hcomp[1, m] + fx * (hcomp[1, m] - env.envx[jm])
                hcomp[3, j] = hcomp[3, m] + fx * (hcomp[3, m] - env.envz[jm])
                for k in range(4, jend + 1):
                    hcomp[k, j] = hcomp[k, m]
                hcomp[2, j] = 1.0 - hcomp[1, j] - hcomp[3, j] - hcomp[4, j]
                lc[j] = bool(env.lcenv[jm] or lc[m])
            else:
                # Between the last two envelope points.
                x0 = env.envs[jm - 1] + stotal
                x1 = stotal + senv
                x2 = env.envs[jm] + stotal
                if x2 - x0 < 1.0e-14:
                    raise SystemExit("STOP 9998")
                fx = (x1 - x0) / (x2 - x0)
                hd[j] = env.envd[jm - 1] + fx * (env.envd[jm] - env.envd[jm - 1])
                hp[j] = env.envp[jm - 1] + fx * (env.envp[jm] - env.envp[jm - 1])
                hr[j] = env.envr[jm - 1] + fx * (env.envr[jm] - env.envr[jm - 1])
                hs[j] = x1
                ht[j] = env.envt[jm - 1] + fx * (env.envt[jm] - env.envt[jm - 1])
                hcomp[1, j] = (env.envx[jm - 1]
                               + fx * (env.envx[jm] - env.envx[jm - 1]))
                hcomp[3, j] = (env.envz[jm - 1]
                               + fx * (env.envz[jm] - env.envz[jm - 1]))
                for k in range(4, jend + 1):
                    hcomp[k, j] = hcomp[k, m]
                hcomp[2, j] = 1.0 - hcomp[1, j] - hcomp[3, j] - hcomp[4, j]
                lc[j] = bool(env.lcenv[jm] or env.lcenv[jm - 1])
            m = j
            break                                  # go to 587
    else:
        # Every envelope point is below the fitting point: the boundary
        # is at the photosphere.
        m = m + env.numenv

    # label 587: unlogged masses (hs1) and shell masses (hs2) (Fortran
    # lines 259–268).
    for j in range(mm, m + 1):
        hs1[j] = 10.0**hs[j]
    for j in range(mm, m):
        hs2[j] = 0.5 * (hs1[j + 1] - hs1[j - 1])
    hsbase = 0.5 * (hs1[m] + hs1[m - 1])
    hs2[m] = 10.0**hstot - hsbase

    # Rotation of the new points (Fortran lines 269–314).
    if nml.lrot:
        # First guess: omega = C * r**walpcz, the law WCZ imposes.
        if nml.walpcz <= -2.0:
            # Constant j.
            for j in range(mm + 1, m + 1):
                hjm[j] = hjm[mm]
                hi[j] = c.cc23 * 10.0**(2.0 * hr[j])
                omega[j] = hjm[j] / hi[j]
        elif nml.walpcz >= 0.0:
            # Solid body.
            # Suspected bug (see the page): omega is local; omega[mm] is
            # stale.
            for j in range(mm + 1, m + 1):
                omega[j] = omega[mm]
                hi[j] = c.cc23 * 10.0**(2.0 * hr[j])
                hjm[j] = omega[mm] * hi[j]
        else:
            # General case (MHP 10/29/26). Suspected bug (see the page):
        # omega[mm] here too.
            w0 = omega[mm] / 10.0**(hr[mm] * nml.walpcz)
            for j in range(mm + 1, m + 1):
                omega[j] = w0 * 10.0**(hr[j] * nml.walpcz)
                hi[j] = c.cc23 * 10.0**(2.0 * hr[j])
                hjm[j] = omega[j] * hi[j]
        getrot(hd, hjm, hr, hs, hs2, lc, m, eta2, hi, omega, qiw, r0)
        # Angular momentum and rotational kinetic energy.
        sumj = 0.0
        sumke = 0.0
        for i in range(1, m + 1):                  # loop 550
            hj = hjm[i] * hs2[i]
            hkerot[i] = 0.5 * omega[i] * hj
            sumj = sumj + hj
            sumke = sumke + hkerot[i]
        sjtot = sumj
        skerot = sumke

    # Question (see the page): written to standard output (Fortran lines
    # 315–324).
    print(" NEW INTERIOR POINTS FROM CHANGE IN ENVELOPE MASS")
    print(" J,LOG RHO, LOG L, LOG P, LOG R, LOG M, LOG T, CONV T/F")
    for j in range(mm, m + 1):
        print(f"{j:5d}{hd[j]:16.8E}{hl[j]:16.8E}{hp[j]:16.8E}{hr[j]:16.8E}"
              f"{hs[j] - stotal:16.8E}{ht[j]:16.8E} {'T' if lc[j] else 'F'}")
    lnew = True
    print("     ***** NEW ENVELOPE MASS CALCULATED *****")
    print(f"        OLD SENV {senv1:22.13E}  NEW SENV{cmp.senv:22.13E}")

    return senvold, sjtot, skerot, m, lnew
