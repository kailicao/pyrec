"""Python version of calcad.f (YREC 7be4133): acoustic depths of the star
and of its convection zone."""

import math

import numpy as np

from pyrec import acdpth_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import envstruct_mod
from pyrec import eos06misc_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.boole import boole
from pyrec.envint import envint
from pyrec.eqstat2 import eqstat2
from pyrec.esac06 import esac06
from pyrec.splinj import splinj
from pyrec.splint import splint

JSON = dimensions_mod.json

# CALCAD has SAVE. Kept between calls because later calls can read them:
# ICZ (when no pair of points brackets the base of the convection zone, or
# when JENV >= M) and the convection-zone profile written with LCLCD.
icz = 0
starrz = np.zeros(JSON + 1)
starcz = np.zeros(JSON + 1)
# The step sizes saved on entry and never restored. Suspected bug (see the
# page).
abeg0 = amin0 = amax0 = ebeg0 = emin0 = emax0 = 0.0


def calcad(hr: np.ndarray, rczl: float, m: int, hd: np.ndarray, hp: np.ndarray, ht: np.ndarray, bl: float,
           fp: np.ndarray, ft: np.ndarray, hstot: float, teffl: float, hcomp: np.ndarray, dage: float, jenv: int,
           cmixl: float) -> None:
    """Join the interior model and a fresh envelope integration into one
    profile, compute 1/c_s from the equation of state, and integrate it to
    get the acoustic depth of the star (TNORM) and of its convection zone
    (TCZ), and their ratio TAUCZN.

    Translates: calcad.f#CALCAD (7be4133)

    Arguments: the interior model hr, hd, hp, ht (log r, log rho, log P,
    log T; arrays, index 0 unused); rczl, log r_CZ / R_sun (RCZL(1) in the
    Fortran, a scalar here); m, the last interior point; bl, log L; fp, ft,
    the rotational distortion; hstot; teffl, log Teff; hcomp[species,
    shell]; dage, the age; jenv, the base of the convection zone; cmixl,
    the mixing length (passed on to ENVINT).

    Returns nothing: the results go to acdpth_mod (tnorm, tcz, tauczn,
    tatmos, deladj).

    """

    global icz, abeg0, amin0, amax0, ebeg0, emin0, emax0
    c = constants_mod
    nml = nmlpar_mod
    acd = acdpth_mod
    env = envstruct_mod
    eos = eos06misc_mod.eos

    starr = np.zeros(JSON + 1)
    start = np.zeros(JSON + 1)
    stard = np.zeros(JSON + 1)
    starx = np.zeros(JSON + 1)
    starp = np.zeros(JSON + 1)
    starz = np.zeros(JSON + 1)
    starc = np.zeros(JSON + 1)
    gm1 = np.zeros(JSON + 1)
    hood1 = np.zeros(8)
    hood2 = np.zeros(8)
    d = np.zeros(8)
    dum1 = np.zeros(5)
    dum2 = np.zeros(4)
    dum3 = np.zeros(4)
    dum4 = np.zeros(4)
    fxionj = np.zeros(4)

    # Initialize the results (Fortran lines 180–186).
    acd.tauczn = 0.0
    acd.tcz = 0.0
    acd.tnorm = 0.0
    acd.whe = 0.0
    acd.tauhe = 0.0
    acd.tatmos = 0.0

    # The envelope, integrated with fixed step sizes (Fortran lines
    # 196–207). Suspected bug (see the page): the step sizes are never
    # restored.
    abeg0, amin0, amax0 = nml.atmbeg, nml.atmmin, nml.atmmax
    ebeg0, emin0, emax0 = nml.envbeg, nml.envmin, nml.envmax
    nml.atmbeg = nml.atmmin = nml.atmmax = nml.atmstp
    nml.envbeg = nml.envmin = nml.envmax = nml.envstp

    # The arguments of ENVINT (Fortran lines 208–223).
    idumj = 0
    bj = math.exp(c.cln * bl)
    fplj = fp[m]
    ftlj = ft[m]
    katmj = 0
    kenvj = 0
    ksahaj = 0
    lpulptj = True
    ixxj = 0
    lprtj = True
    lsbc0j = False
    xj = hcomp[1, m]
    zj = hcomp[3, m]
    rlj = 0.5 * (bl + c.clsunl - 4.0 * teffl - c.c4pil - c.csigl)
    glj = c.cgl + comp_mod.stotal - rlj - rlj
    plimj = hp[m]

    # For a spotted star, the ambient temperature (Fortran lines 228–232).
    # Question (see the page): spots only when JENV = M.
    if jenv == m and nml.spotf != 0.0 and nml.spotx != 1.0:
        ateffl = teffl - 0.25 * math.log10(
            nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
    else:
        ateffl = teffl

    # The envelope integration, which fills envstruct_mod (Fortran lines
    # 233–235).
    rlj, idumj, katmj, kenvj, ksahaj = envint(
        bj, fplj, ftlj, glj, hstot, ixxj, lprtj, lsbc0j, plimj, rlj, ateffl,
        xj, zj, dum1, idumj, katmj, kenvj, ksahaj, dum2, dum3, dum4,
        lpulptj, cmixl)
    numenv = env.numenv

    # One profile: interior points 1..M, then envelope points 2..NUMENV,
    # in linear units, T in 10^6 K (labels 34, 35; Fortran lines 241–256).
    for i in range(1, m + 1):
        starr[i] = math.exp(c.cln * hr[i])
        start[i] = math.exp(c.cln * ht[i]) / 1.0e6
        stard[i] = math.exp(c.cln * hd[i])
        starp[i] = math.exp(c.cln * hp[i])
        starx[i] = hcomp[1, i]
        # Suspected bug (see the page): the envelope's Z at every interior
        # point.
        starz[i] = hcomp[3, m]
    for i in range(1, numenv):
        starr[m + i] = math.exp(c.cln * env.envr[i + 1])
        start[m + i] = math.exp(c.cln * env.envt[i + 1]) / 1.0e6
        stard[m + i] = math.exp(c.cln * env.envd[i + 1])
        starp[m + i] = math.exp(c.cln * env.envp[i + 1])
        starx[m + i] = hcomp[1, m]
        starz[m + i] = hcomp[3, m]

    # Inverse sound speed and the adiabatic gradient, from the OPAL 2006
    # EOS (label 36; Fortran lines 261–275).
    if nml.lopale06:
        iorder = 9
        irad = 1
        for i in range(1, m + numenv):
            # The Fortran's alternate return (*999) lands on the next
            # statement, so a failed call goes on with the previous EOS.
            # Suspected bug (see the page): EOS errors ignored.
            esac06(starx[i], start[i], stard[i], iorder, irad)
            starc[i] = 1.0 / math.sqrt(eos[8] * starp[i] / stard[i])
            acd.deladj[i] = 1.0 / eos[9]
            gm1[i] = eos[8]

    # Otherwise from EQSTAT2, with Gamma_1 from Cox and Giuli, chapter 9
    # (label 41; Fortran lines 277–304).
    if not nml.lopale06:
        for i in range(1, m + numenv):
            stx = starx[i]
            stt = start[i] * 1.0e6
            std = stard[i]
            stdl = math.log10(stard[i])
            stz = starz[i]
            sttl = math.log10(stt)
            stp = starp[i]
            stpl = math.log10(stp)
            lcalcj = True  # noqa: F841  (unused in the Fortran too)
            lderij = False
            latmoj = True
            # Suspected bug (see the page): FXIONJ is a scalar in CALCAD but
            # FXION(3) in EQSTAT2, which writes past it; here an array, so
            # that EQSTAT2 can fill it.
            (stt, stp, stdl, std, betaj, betaij, bet14j, rmuj, amuj, emuj,
             etaj, qdtj, qdpj, qcpj, delaj, qdttj, qdtpj, qatj, qapj, qcptj,
             qcppj, ksahaj) = eqstat2(
                sttl, stt, stpl, stp, stdl, std, stx, stz, 0.0, 0.0, 0.0,
                fxionj, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, lderij, latmoj, ksahaj)
            chrh = 1.0 / qdpj
            cht = -chrh * qdtj
            cvj = qcpj - math.exp(c.cln * (stpl - stdl - sttl)) * cht**2 / chrh
            gm1[i] = chrh * qcpj / cvj
            starc[i] = 1.0 / math.sqrt(gm1[i] * starp[i] / stard[i])
            acd.deladj[i] = delaj

    # t_norm: Boole's rule on 1 + 4p points (Fortran lines 310–322).
    # Question (see the page): K leaves out the outermost point.
    k = numenv - 1 + m - 1
    u = (numenv - 1 + m - 1) % 4
    if u != 0:
        nn = k + 4 - u + 1
    else:
        nn = k + 1
    # Suspected bug (see the page): a short profile gives TNORM = 0.
    if k <= 3:
        k = 4
        nn = 1 + 1
    acd.tnorm = boole(starr, starc, k, nn, acd.tnorm)

    # t_CZ: from r_CZ to the surface, if there is a convection zone
    # (Fortran lines 328–385).
    if jenv < m:
        rcz = math.exp(c.cln * (rczl + c.crsunl))
        ans = 0.0
        klo = khi = 0
        flg = 0

        # Find the point below r_CZ and put r_CZ there, with 1/c_s from a
        # spline through 7 points (loop 22; Fortran lines 333–369).
        for i in range(1, k + 1):
            if starr[i] <= rcz and starr[i + 1] >= rcz:
                icz = i
            # Suspected bug (see the page): ICZ may be stale, from an
            # earlier call, and then can match I before the bracket is
            # found.
            if i == icz and flg == 0 and icz > 3:
                flg = 1
                for p in range(1, 8):
                    hood1[p] = starr[icz - 4 + p]
                    hood2[p] = starc[icz - 4 + p]
                splinj(hood1, hood2, d, 7)
                ans, klo, khi = splint(hood1, hood2, 7, d, rcz, ans, klo, khi)
                starr[icz] = rcz
                starc[icz] = ans
            if i == icz and flg == 0 and icz <= 3:
                flg = 1
                if icz == 1:
                    rcz = starr[icz]
                    ans = starc[icz]
                else:
                    for p in range(1, 8):
                        hood1[p] = starr[icz - 2 + p]
                        hood2[p] = starc[icz - 2 + p]
                    splinj(hood1, hood2, d, 7)
                    ans, klo, khi = splint(hood1, hood2, 7, d, rcz, ans,
                                           klo, khi)
                    starr[icz] = rcz
                    starc[icz] = ans

        # The profile from r_CZ outward (loop 27; Fortran lines 370–373).
        for i in range(1, k - icz + 2):  # label 27
            starrz[i] = starr[icz - 1 + i]
            starcz[i] = starc[icz - 1 + i]

        # Boole's rule from r_CZ to the surface (Fortran lines 375–382).
        u = (k - icz + 1) % 4
        v = k - icz + 1
        if u != 0:
            nn = k - icz + 1 + 4 - u + 1
        else:
            nn = k - icz + 1 + 1
        acd.tcz = boole(starrz, starcz, v, nn, acd.tcz)
    else:
        acd.tcz = 0.0

    # The normalized acoustic depth (Fortran lines 440–441). Suspected bug
    # (see the page): TNORM = 0 for a short profile; Python raises
    # ZeroDivisionError where the Fortran gives Inf or NaN.
    acd.tauczn = 0.0
    acd.tauczn = acd.tcz / acd.tnorm

    # The line to ISHORT (Fortran lines 444–446). Suspected bug (see the
    # page): format 67. Question (see the page): TATMOS is always 0.
    iounits_mod.write(
        iounits_mod.ishort,
        f" Acoustic depth to CZ:{acd.tcz:14.8f}  Acoustic depth to center"
        f"{acd.tnorm:13.7f}  Normalized taucz:{acd.tauczn:11.9f}"
        f"Acoustic depth of atmopshere:{acd.tatmos:16.8f}  ")

    # With LCLCD, the full profile to ICLCD (Fortran lines 451–486).
    if acd.lclcd:
        iounits_mod.write(
            iounits_mod.iclcd,
            f" Number of points in envelope:{numenv:5d}  "
            f"Number of points in interior:{m:5d}  Index near Rcz:{icz:5d}  ")
        for i in range(1, m + numenv):  # label 1505
            if i <= k - icz + 1:
                row = (dage, starr[i], starc[i], starrz[i], starcz[i],
                       acd.deladj[i], gm1[i], starp[i], start[i], stard[i],
                       starx[i])
            else:
                row = (dage, starr[i], starc[i], 0.0, 0.0, acd.deladj[i],
                       gm1[i], starp[i], start[i], stard[i], starx[i])
            iounits_mod.write(iounits_mod.iclcd,
                              " " + "".join(f"{x:16.8E}" for x in row))
