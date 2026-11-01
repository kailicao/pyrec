"""Python version of crrect.f (YREC 09c3b9b): the Henyey iterations that
correct the interior model at one level of the time step."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import ctol_mod
from pyrec import debhu_mod
from pyrec import dimensions_mod  # JSON, the array sizes (USEd by the Fortran)
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec.coefft import coefft
from pyrec.fpft import fpft
from pyrec.getrot import getrot
from pyrec.hsolve import hsolve
from pyrec.mix import mix
from pyrec.surfbc import surfbc

# CRRECT has SAVE. Locals kept between calls that matter: JCORE, JENV and
# MXZON0 are set only by MIX (at ITLVL > 2), but JENV is passed to COEFFT at
# every level; at levels 1 and 2 it keeps its value from an earlier call
# (unset at the first call; 0 here).
jcore = 0
jenv = 0
mxzon0 = np.zeros((13, 3), dtype=int)  # MXZON0(12,2), index 0 unused
ddage = 0.0
lflp = False


def _l(flag):
    """Fortran's L format of a logical."""
    return "T" if flag else "F"


def crrect(delts, m, niter, lfini, large, lnew, lreset, lsbc,
           trit, tril, cfenv, ps, ts, rs, ftri, istore, estore, hstot, bl,
           teffl, hd, hha, hhb, hhc, hl, hmax, hp, hpp, hr, hs, hs1, hs2,
           hsrf, ht, htt, hcomp, lc, nhmax, tlumx, latmo, lderiv, lmix,
           locond, qdt, qdp, fp, ft, eta2, omega, r0, itdone, hg, hi, hjm,
           itlvl, lcz, mxzone, qiw, hkerot, hkerot0, cmixl):
    """Henyey iteration for one level of the time step: the surface
    boundary condition (SURFBC), mixing and burning (MIX), the equations
    (COEFFT), their solution (HSOLVE), damped corrections and the
    convergence tests against HTOLER.

    Translates: crrect.f#CRRECT (09c3b9b)

    Arguments: delts, the time step (s); m, the number of shells; niter,
    the maximum number of iterations; lfini (read: converged at the
    previous level); the flags lnew, lreset, lsbc; the envelope triangle
    trit, tril, ps, ts, rs, ftri, istore, estore and the fit cfenv[1..9];
    hstot, bl, teffl; the model hd, hl, hp, hpp, hr, hs, hs1, hs2, ht, htt,
    hcomp[species, shell], lc; the Henyey work arrays hha, hhb, hhc, hsrf,
    hmax, nhmax; tlumx[1..8], the luminosity by source; latmo, lderiv,
    lmix, locond, qdt, qdp (for COEFFT); rotation fp, ft, eta2, omega, r0,
    hg, hi, hjm, lcz, qiw, hkerot, hkerot0; itdone, the iterations done;
    itlvl, the level of iteration; mxzone[12, 2], the convection zones;
    cmixl, the mixing length. Arrays have index 0 unused and are changed
    in place.

    Returns (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
    latmo, lderiv, lmix, locond, qdt, qdp, itdone): converged; corrections
    too large; the flags and triangle data SURFBC updates; the new log L
    and log Teff; the flags and derivatives COEFFT sets; the iteration
    count.
    """
    global jcore, jenv, ddage, lflp
    nml = nmlpar_mod
    io = iounits_mod

    if niter <= 0:
        return (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
                latmo, lderiv, lmix, locond, qdt, qdp, itdone)
    bl = math.log10(hl[m])
    # Zero counters.
    kenv = 0
    katm = 0
    ksaha = 0
    comp_mod.senv = hs[m] - hstot
    # Question (see the page): reads as lnew or (lreset and itlvl == 2).
    if lnew or (lreset and itlvl == 2):
        lsbc = True
    # Find new FP and FT if the model is rotating.
    if nml.lrot and lsbc:
        fpl = fp[m]
        ftl = ft[m]
    else:
        fpl = 1.0
        ftl = 1.0
    # Surface boundary condition: envelope composition, then SURFBC.
    if lsbc:
        x = comp_mod.xnew
        z = comp_mod.znew
        plim = hp[m]
        if nml.ldh:
            debhu_mod.xxdh = hcomp[1, m]
            debhu_mod.yydh = hcomp[2, m] + hcomp[4, m]
            debhu_mod.zzdh = hcomp[3, m]
            debhu_mod.zdh[1] = hcomp[5, m] + hcomp[6, m]
            debhu_mod.zdh[2] = hcomp[7, m] + hcomp[8, m]
            debhu_mod.zdh[3] = hcomp[9, m] + hcomp[10, m] + hcomp[11, m]
        (ftri, istore, lnew, lreset, ksaha, kenv, katm, teffl,
         lflp) = surfbc(trit, tril, cfenv, ps, ts, rs, ftri, istore, estore,
                        lnew, lreset, ksaha, kenv, katm, hstot, bl,
                        teffl, x, z, fpl, ftl, lflp, plim, lc, m, cmixl)
    # Reset the composition to the start of the step, then MIX (7/91).
    if itlvl > 2 and delts > 0.0:
        nspec = 11
        if nml.lexcom:
            nspec = 15
        for i in range(1, m + 1):
            for j in range(1, nspec + 1):
                hcomp[j, i] = oldmod_mod.hcompp[j, i]
        ddage, jcore, jenv = mix(delts, hcomp, hd, hl, hp, hr, hs, hs1, hs2,
                                 hstot, ht, itlvl, lc, m, ddage, jcore, jenv,
                                 mxzone, mxzon0, teffl, cmixl)
    # Question (see the page): the call to OVROT is commented out.
    # Converged at the previous level and the triangle not flipped: done.
    if lsbc and not lflp and lfini:
        return (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
                latmo, lderiv, lmix, locond, qdt, qdp, itdone)
    lfini = False
    for iter_ in range(1, niter + 1):  # loop 100
        # Henyey reduction.
        hmax[1] = 0.0
        hmax[2] = 0.0
        hmax[3] = 0.0
        hmax[4] = 0.0
        # Question (see the page): jenv is stale at levels 1 and 2.
        (latmo, lderiv, lmix, locond, qdt, qdp,
         ksaha) = coefft(delts, m, hd, hha, hhb, hhc, hl, hmax, hp, hpp, hr,
                         hs, hs1, hs2, ht, htt, hcomp, lc, tlumx, latmo,
                         lderiv, lmix, locond, qdt, qdp, ksaha, fp, ft,
                         hkerot, hkerot0, jenv, teffl, cmixl)
        # Renormalize TLUMX (not during the helium flash).
        thl = (tlumx[1] + tlumx[2] + tlumx[3] + tlumx[4] + tlumx[5]
               + tlumx[6] + tlumx[7] + tlumx[8])
        if not nml.lkuthe and thl > 0.0:
            temp = hl[m] / thl
            for j in range(1, 9):
                tlumx[j] = tlumx[j] * temp
        # Residuals of the P and T equations (not the boundary equations).
        # Question (see the page): converged without solving.
        if iter_ > 1:
            if (hmax[1] <= nml.htoler[5, 1] and hmax[2] <= nml.htoler[5, 2]
                    and hmax[3] <= nml.htoler[5, 1]
                    and hmax[4] <= nml.htoler[5, 2]):
                io.write(io.ishort, " R.H.S. BELOW TOLERANCES--P"
                         f"{hmax[1]:9.2E}  T {hmax[2]:9.2E}"
                         f"  R {hmax[3]:9.2E}  L {hmax[4]:9.2E}")  # format 20
                io.write(io.ishort, f"{'':10}MODEL CONVERGED AFTER "
                         f"{itdone + 1:4d}  ITERATIONS")  # format 75
                lfini = True
                return (lfini, large, lnew, lreset, lsbc, ftri, istore, bl,
                        teffl, latmo, lderiv, lmix, locond, qdt, qdp, itdone)
        # Henyey coefficients for the surface, from the envelope fit.
        hsrf[1] = -cfenv[1]
        hsrf[2] = -cfenv[2]
        hsrf[3] = cfenv[1] * hp[m] + cfenv[2] * ht[m] + cfenv[3] - hr[m]
        temp = constants_mod.cln * hl[m]
        hsrf[4] = -temp * cfenv[4]
        hsrf[5] = -temp * cfenv[5]
        hsrf[6] = temp * (cfenv[4] * hp[m] + cfenv[5] * ht[m] + cfenv[6] - bl)
        # Back solution for the corrections.
        hsolve(m, hha, hhb, hsrf)
        # Largest corrections, and where.
        for j in range(1, 5):
            hmax[j] = abs(hhb[j, 1])
            nhmax[j] = 1
        for i in range(2, m + 1):
            test = abs(hhb[1, i])
            if hmax[1] <= test:
                hmax[1] = test
                nhmax[1] = i
            test = abs(hhb[2, i])
            if hmax[2] <= test:
                hmax[2] = test
                nhmax[2] = i
            test = abs(hhb[3, i])
            if hmax[3] <= test:
                hmax[3] = test
                nhmax[3] = i
            # Question (see the page): the mixed absolute/relative L test.
            test = min(abs(hhb[4, i]), abs(hhb[4, i] / hl[i]))
            if hmax[4] <= test:
                hmax[4] = test
                nhmax[4] = i
        hlmax = hmax[4]
        # Converged below HTOLER(.,1); too large above HTOLER(.,2).
        lfini = (hmax[1] < nml.htoler[1, 1] and hmax[2] < nml.htoler[2, 1]
                 and hmax[3] < nml.htoler[3, 1]
                 and hmax[4] < nml.htoler[4, 1])
        large = (hmax[1] > nml.htoler[1, 2] or hmax[2] > nml.htoler[2, 2]
                 or hmax[3] > nml.htoler[3, 2] or hmax[4] > nml.htoler[4, 2])
        for j in range(1, 5):
            ij = nhmax[j]
            hmax[j] = hhb[j, ij]
        if nml.fcorr0 > 0.0:
            ctol_mod.fcorr = min(1.0, ctol_mod.fcorr + nml.fcorri)
        # Helium flash: damping by the size of the L correction.
        fact = ctol_mod.fcorr
        if nml.lkuthe:
            if hlmax <= 5.0e-1:
                fact = 8.0e-1
            if hlmax <= 5.0e-3:
                fact = 1.0
        tlumh = tlumx[1] + tlumx[2] + tlumx[3] + tlumx[4]
        tlumhe = tlumx[5] + tlumx[8]
        if nml.lcorr:
            io.write(io.ishort, f" DEL-L/L  {_l(lfini):>2}{hmax[4]:12.4E}"
                     f"{nml.htoler[4, 1]:12.4E}{'':5}{nhmax[4]:5d}")  # fmt 60
            io.write(io.ishort, " CORR"
                     f"{nhmax[1]:5d}P{hmax[1]:9.2E}{nhmax[2]:5d}T{hmax[2]:9.2E}"
                     f"{nhmax[3]:5d}R{hmax[3]:9.2E}{nhmax[4]:5d}L{hmax[4]:9.2E}"
                     f"  F={fact:5.3f}  E-HY{tlumh:10.3E} HE{tlumhe:10.3E}"
                     f" NU{tlumx[6]:10.3E} G{tlumx[7]:10.3E}")  # format 70
            if lfini:
                io.write(io.ishort, f"{'':10}MODEL CONVERGED AFTER "
                         f"{itdone + 1:4d}  ITERATIONS")  # format 75
        if large:
            io.write(io.iowr, " -----CORRECTIONS EXCEEDED TOLERANCES")
            io.write(io.ishort, " -----CORRECTIONS EXCEEDED TOLERANCES")
            return (lfini, large, lnew, lreset, lsbc, ftri, istore, bl,
                    teffl, latmo, lderiv, lmix, locond, qdt, qdp, itdone)
        # Apply the corrections, damped by fact.
        for i in range(1, m + 1):
            temp = fact * hhb[1, i]
            hp[i] = hp[i] + temp
            hpp[i] = hpp[i] + temp
            temp = fact * hhb[2, i]
            ht[i] = ht[i] + temp
            htt[i] = htt[i] + temp
            hr[i] = hr[i] + fact * hhb[3, i]
            hl[i] = hl[i] + fact * hhb[4, i]
        bl = math.log10(hl[m])
        teffl = cfenv[7] * hp[m] + cfenv[8] * ht[m] + cfenv[9]
        if nml.lrot:
            getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
            fpft(hd, hr, hs, m, omega, eta2, fp, ft, hg, r0)
            for i in range(1, m + 1):
                hj = hjm[i] * hs2[i]
                hkerot[i] = 0.5 * omega[i] * hj
        itdone = itdone + 1
        if lfini:
            break
    return (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
            latmo, lderiv, lmix, locond, qdt, qdp, itdone)
