"""Python version of sconvec.f (YREC 7be4133): semiconvection at the
edges of convection zones (Castellani et al. 1971), and the merging of
zones that come to overlap."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import debhu_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.eqstat import eqstat
from pyrec.getopac import getopac
from pyrec.tpgrad import tpgrad

# Local kept between calls (SAVE): the ionization fractions from eqstat.
fxion = np.zeros(4)


def _load_debhu(hcomp: np.ndarray, j: int) -> None:
    """Load shell j's composition for the Debye-Hueckel correction in
    eqstat (the block the Fortran repeats three times)."""

    if nmlpar_mod.ldh:
        debhu_mod.xxdh = hcomp[1, j]
        debhu_mod.yydh = hcomp[2, j] + hcomp[4, j]
        debhu_mod.zzdh = hcomp[3, j]
        debhu_mod.zdh[1] = hcomp[5, j] + hcomp[6, j]
        debhu_mod.zdh[2] = hcomp[7, j] + hcomp[8, j]
        debhu_mod.zdh[3] = hcomp[9, j] + hcomp[10, j] + hcomp[11, j]


def _gradients(tl: float, pl: float, dl: float, rl: float, sl: float, b: float, x: float, z: float, lderiv: bool,
               latmo: bool, ksaha: int, fpl: float, ftl: float, teffl: float,
               cmixl: float) -> tuple[float, float, float, float, float]:
    """EQSTAT, GETOPAC and TPGRAD at one point (repeated three times in
    the Fortran). Returns (dl, amu, emu, dela, delr)."""

    (t, p, dl, d, beta, beta1, beta14, rmu, amu, emu, eta, qdt, qdp,
     qcp, dela, qdtt, qdtp, qat, qap, qcpt, qcpp, ksaha) = eqstat(
        tl, 0.0, pl, 0.0, dl, 0.0, x, z, 0.0, 0.0, 0.0, fxion, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        lderiv, latmo, ksaha)
    o, ol, qod, qot = getopac(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion)
    (dele, delr, qact, qacp, qacr, vel, lconv) = tpgrad(
        tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, 0.0,
        0.0, dela, qdtt, qdtp, qat, qap, 0.0, 0.0, 0.0, qcpt, qcpp,
        0.0, lderiv, False, fpl, ftl, teffl, cmixl)
    return dl, amu, emu, dela, delr


def sconvec(delts: float, hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray,
            hs: np.ndarray, ht: np.ndarray, m: int, mxzone: np.ndarray, nzone: int, teffl: float, cmixl: float) -> int:
    """Test the edges of the first convection zone (L_ONLY_CORE) for
    semiconvection and extend the zone where the region next to it is
    unstable, as far as the propagation speed allows in one step; then
    merge zones that overlap.

    Translates: sconvec.f#SCONVEC (7be4133)

    Arguments: the step delts (s); hcomp[species, shell]; the model's hd
    (log rho; overwritten at the shells tested), hl (L/Lsun), hp, hr, hs,
    ht (logs); m shells; the zones mxzone[k, 1..2] (changed in place) and
    their count nzone; teffl; cmixl, the mixing length (for tpgrad).

    Returns nzone (reduced by mergers).

    """

    c = constants_mod

    # Only the first zone (Fortran lines 89–98).
    # Question (see the page): only the first zone is tested.
    l_only_core = True
    if l_only_core:
        iloop = 1
    else:
        iloop = nzone
    for i in range(1, iloop + 1):                   # loop 210
        # Core, surface or intermediate zone: LFLAGU (test above),
        # LFLAGD (test below) (Fortran lines 101–119).
        if mxzone[i, 1] == 1:
            # Convective core; skip a fully convective star.
            if mxzone[i, 2] == m:
                return nzone
            lflagu = True
            lflagd = False
        elif mxzone[i, 2] == m:
            lflagu = False
            lflagd = True
        else:
            lflagu = True
            lflagd = True

        for k in (1, 2):                            # loop 200
            # Fortran lines 120–134.
            if k == 1 and not lflagd:
                continue
            if k == 2 and not lflagu:
                continue
            # The edge of the zone (JMC) and the first radiative shell
            # outside it (JMR).
            jmc = mxzone[i, k]
            if k == 1:
                jmr = jmc - 1
            if k == 2:
                jmr = jmc + 1
            lderiv = False
            locond = True  # noqa: F841  (unused in the Fortran too)
            latmo = True
            ksaha = 0

            # The structure of JMR with the zone's composition (Fortran
            # lines 135–168).
            b = hl[jmr]
            sl = hs[jmr]
            pl = hp[jmr]
            tl = ht[jmr]
            dl = hd[jmr]
            rl = hr[jmr]
            x = hcomp[1, jmc]
            z = hcomp[3, jmc]
            fpl = 1.0
            ftl = 1.0
            iu = jmc
            _load_debhu(hcomp, jmc)
            dl, amu, emu, dela, delr = _gradients(
                tl, pl, dl, rl, sl, b, x, z, lderiv, latmo, ksaha,
                fpl, ftl, teffl, cmixl)
            # Skip if stable with the zone's composition.
            if delr < dela:
                continue

            # Fortran lines 169–173.
            xintmu = amu + emu
            delr1 = delr
            rmax = 1.0 - dela / delr

            # Repeat with the local composition (Fortran lines 174–194).
            # Question (see the page): hd overwritten.
            x = hcomp[1, jmr]
            z = hcomp[3, jmr]
            iu = jmr
            _load_debhu(hcomp, jmr)
            dl, amu, emu, dela, delr = _gradients(
                tl, pl, dl, rl, sl, b, x, z, lderiv, latmo, ksaha,
                fpl, ftl, teffl, cmixl)
            hd[jmr] = dl

            # Ratio of the gradients with the zone's and the local
            # composition (Fortran lines 195–197).
            fdel = delr1 / delr
            delr0 = delr

            # RMAX: the step times the propagation speed (eq. 5') (Fortran
            # lines 198–202).
            # Question (see the page): the factor 10.
            rmax = (rmax * (hl[jmc] * c.clsun
                            / (10.0 * c.c4pi * math.exp(c.cln * hp[jmc])))
                    * (delts / math.exp(c.cln * (hr[jmc] + hr[jmc]))))

            # The walk's start (Fortran lines 203–217).
            rsum = 0.0
            lmax = False
            if k == 1:
                rval1 = math.exp(c.cln * hr[jmr + 1])
                jbeg = jmr - 1
                jend = 1
                jstep = -1
            elif k == 2:
                rval1 = math.exp(c.cln * hr[jmr - 1])
                jbeg = jmr + 1
                jend = m
                jstep = 1

            # The walk, shell by shell (Fortran lines 218–259).
            stable = False
            for ii in range(jbeg, jend + jstep, jstep):     # loop 32
                # The limit in radius (Fortran lines 219–225).
                if k == 1:
                    rval2 = math.exp(c.cln * hr[ii + 1])
                if k == 2:
                    rval2 = math.exp(c.cln * hr[ii - 1])
                rdiff = rval2 - rval1
                rsum = rsum + (1.0 - xintmu / (amu + emu)) * rdiff
                rval1 = rval2
                if rsum > rmax:
                    break                                   # go to 33

                # The rescaled radiative gradient (Fortran lines 226–254).
                # Suspected bug (see the page): RL is not updated here.
                b = hl[ii]
                sl = hs[ii]
                pl = hp[ii]
                tl = ht[ii]
                dl = hd[ii]
                x = hcomp[1, ii]
                z = hcomp[3, ii]
                iu = ii  # noqa: F841  (unused in the Fortran too)
                _load_debhu(hcomp, ii)
                dl, amu, emu, dela, delr = _gradients(
                    tl, pl, dl, rl, sl, b, x, z, lderiv, latmo, ksaha,
                    fpl, ftl, teffl, cmixl)
                hd[ii] = dl
                if fdel * delr < dela:
                    stable = True
                    break                                   # go to 34
            else:
                # The loop ran out: the edge of the star.
                if k == 1:
                    ii = 0
                if k == 2:
                    ii = m + 1
            if not stable:
                lmax = True                                 # label 33

            # label 34: the new edge (Fortran lines 260–270).
            if k == 1:
                mxzone[i, k] = ii + 1
                jmc2 = ii + 1
            elif k == 2:
                mxzone[i, k] = ii - 1
                jmc2 = ii - 1
            iounits_mod.write(
                iounits_mod.ishort,
                f" CZ OLD EDGE {jmc:3d} EXTENDED TO-{jmc2:3d} LIMIT="
                f"{'T' if lmax else 'F'} RAD.GRADS-IN/OUT"
                f"{delr1:8.4f}{delr0:8.4f}")                # format 601

    # Mergers of zones that overlap after semiconvection (Fortran lines
    # 273–301).
    if nzone <= 1:
        return nzone
    j = 1
    while True:                                             # label 85
        # Question (see the page): overlap only, not adjacency.
        if mxzone[j, 2] > mxzone[j + 1, 1]:
            iounits_mod.write(
                iounits_mod.ishort,
                "  MIXED ZONES MERGED DUE TO SEMICONVECTION\n"
                f"  OLD[{mxzone[j, 1]:3d}-{mxzone[j, 2]:3d}]"
                f"[{mxzone[j + 1, 1]:3d}-{mxzone[j + 1, 2]:3d}]"
                f" NEW[{mxzone[j, 1]:3d}-{mxzone[j + 1, 2]:3d}]")  # 93
            # Suspected bug (see the page): the merged zone keeps zone j+1's
            # top.
            mxzone[j + 1, 1] = mxzone[j, 1]
            for kk in range(j, nzone):
                for j2 in (1, 2):
                    mxzone[kk, j2] = mxzone[kk + 1, j2]
            nzone = nzone - 1
            if j <= nzone - 1:
                continue
            return nzone
        j = j + 1
        if j > nzone - 1:
            return nzone
