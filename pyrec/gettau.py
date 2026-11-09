"""Python version of gettau.f (YREC development, 7be4133): the convective
overturn timescale near the base of the surface convection zone."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import envstruct_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec import scrtch_mod
from pyrec.envint import envint
from pyrec.tauint import tauint
from pyrec.tauintnew import tauintnew

JSON = dimensions_mod.json

# GETTAU has SAVE: its work arrays are module-level.
hcompf = np.zeros((16, JSON + 1))
hrf = np.zeros(JSON + 1)
hpf = np.zeros(JSON + 1)
hdf = np.zeros(JSON + 1)
hsf = np.zeros(JSON + 1)
htf = np.zeros(JSON + 1)
hgf = np.zeros(JSON + 1)
hvf = np.zeros(JSON + 1)
delf1 = np.zeros(JSON + 1)
delf2 = np.zeros(JSON + 1)
lcf = np.zeros(JSON + 1, dtype=bool)
dum1 = np.zeros(5)
dum2 = np.zeros(4)
dum3 = np.zeros(4)
dum4 = np.zeros(4)
lallcz = False


def gettau(hcomp: np.ndarray, hr: np.ndarray, hp: np.ndarray, hd: np.ndarray, hs1: np.ndarray, ht: np.ndarray,
           fp: np.ndarray, ft: np.ndarray, teffl: float, hstot: float, bl: float, m: int, lc: np.ndarray, rbcz: float,
           cmixl: float) -> float:
    """Convective overturn time one pressure scale height above the base of
    the surface convection zone ("G Somers, 3/17").

    Uses the interior model, the envelope integration, or the two stitched
    together. Sets ``ovrtrn_mod.taucz`` and ``ovrtrn_mod.pphot``.

    Translates: gettau.f#GETTAU (7be4133)

    Parameters
    ----------
    hcomp, hr, hp, hd, hs1, ht, fp, ft, lc : arrays
        The model (composition, log r, log P, log rho, mass, log T, the
        form factors, convective flags).
    teffl, hstot, bl : float
        log Teff, log total mass, log L.
    m : int
        Number of model points.
    rbcz : float
        Radius of the convection zone's base (solar radii); set by
        TAUINTNEW, left unchanged by the old method.
    cmixl : float
        Mixing-length parameter, passed on to ENVINT.

    Returns
    -------
    rbcz

    """

    global lallcz
    cln = constants_mod.cln
    nml = nmlpar_mod

    # Fortran lines 79–80.
    ovrtrn_mod.taucz = 0.0
    ovrtrn_mod.pphot = 0.0

    # Base of the surface convection zone, scanning inward (Fortran lines
    # 88–93).
    for i in range(m - 1, 0, -1):                       # label 71
        if not lc[i]:
            break                                       # GOTO 81
    else:
        lallcz = True        # unused, as in the Fortran
        i = 0
    imax = i + 1                                        # label 81

    # Keep the envelope's structure if one scale height above the base may
    # lie in the envelope (Fortran lines 95–103).
    # Question (see the page): "at least 3 PSCAS".
    chkprs = hp[imax] - hp[m]
    nml.lcalcenv = chkprs < 1.0 or not nml.lnewtcz

    # Integrate the envelope (Fortran lines 107–132).
    fpl = fp[m]
    ftl = ft[m]
    idum = 0
    b = math.exp(cln * bl)
    katm = 0
    kenv = 0
    ksaha = 0
    ixx = 0
    lprt = True
    lsbc0 = False
    lpulpt = False
    x = hcomp[1, m]
    z = hcomp[3, m]
    rl = 0.5 * (bl + constants_mod.clsunl - 4.0 * teffl - constants_mod.c4pil
                - constants_mod.csigl)
    gl = constants_mod.cgl + hstot - rl - rl
    plim = hp[m]
    spotf = nmlpar_mod.spotf
    spotx = nmlpar_mod.spotx
    if lc[m] and spotf != 0.0 and spotx != 1.0:
        ateffl = teffl - 0.25 * math.log10(spotf * spotx ** 4.0 + 1.0 - spotf)
    else:
        ateffl = teffl
    rl, idum, katm, kenv, ksaha = envint(
        b, fpl, ftl, gl, hstot, ixx, lprt, lsbc0, plim, rl, ateffl, x, z,
        dum1, idum, katm, kenv, ksaha, dum2, dum3, dum4, lpulpt, cmixl)

    # The old method's TAUCZ, set in ENVINT, ends here (Fortran line 140).
    if ovrtrn_mod.taucz == 0.0:
        # The interior, with g = G m / r**2, velocities and gradients from
        # the last physics pass (Fortran lines 143–159).
        # Question (see the page): stale interior velocities.
        mm = m
        for i in range(1, mm + 1):
            hcompf[1:16, i] = hcomp[1:16, i]
            hrf[i] = hr[i]
            hpf[i] = hp[i]
            hdf[i] = hd[i]
            hsf[i] = hs1[i]
            hgf[i] = hs1[i] * math.exp(cln * (constants_mod.cgl - 2.0 * hr[i]))
            htf[i] = ht[i]
            hvf[i] = scrtch_mod.svel[i]
            delf1[i] = scrtch_mod.sdel[1, i]
            delf2[i] = scrtch_mod.sdel[3, i]
            lcf[i] = lc[i]

        # The envelope, points 2 to NUMENV (Fortran lines 161–182).
        if nml.lcalcenv and nml.lnewtcz:
            es = envstruct_mod
            mm = m + es.numenv - 1
            for i in range(m + 1, mm + 1):
                k = i - m + 1
                # HCOMP(J,M+1) is one point past the model; hcompf is read
                # nowhere (see the page).
                hcompf[1:16, i] = hcomp[1:16, m + 1]
                hrf[i] = es.envr[k]
                hpf[i] = es.envp[k]
                hdf[i] = es.envd[k]
                hsf[i] = math.exp(cln * (hstot + es.envs[k]))
                hgf[i] = hsf[i] * math.exp(cln * (constants_mod.cgl - 2.0 * hrf[i]))
                htf[i] = es.envt[k]
                hvf[i] = es.evels[k]
                delf1[i] = es.edels[1, k]
                delf2[i] = es.edels[2, k]
                lcf[i] = es.lcenv[k]

        # The new or the old method (Fortran lines 184–194).
        if nml.lnewtcz:
            rbcz = tauintnew(hsf, lcf, hrf, hpf, hdf, hgf, mm, m, hvf,
                             delf1, delf2, rbcz)
        else:
            tauint(hsf, lcf, hrf, hpf, hdf, hgf, mm, hvf, delf1, delf2)

    # Fortran lines 196–197 (label 100).
    nml.lcalcenv = True
    return rbcz
