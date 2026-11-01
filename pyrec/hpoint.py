"""Python version of hpoint.f (YREC 09c3b9b): rezoning the model at the
start of each step."""

import math
import sys

import numpy as np

from pyrec import burn_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod
from pyrec import dwmax_mod
from pyrec import entrop_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldphy_mod
from pyrec import oldrot2_mod
from pyrec import oldrot_mod
from pyrec import scrtch_mod
from pyrec import temp_mod
from pyrec.fpft import fpft
from pyrec.getrot import getrot
from pyrec.momi import momi
from pyrec.osplin import osplin
from pyrec.ovrot import ovrot
from pyrec.physic import physic
from pyrec.splinc import splinc
from pyrec.splintd2 import splintd2
from pyrec.surfopac import surfopac

JSON = dimensions_mod.json

# The Fortran's local arrays. HPOINT has SAVE, so they keep their values
# between calls; they are module-level here for the same reason (the debug
# table can print entries of HKO and HFO left from an earlier call).
omegal = np.zeros(JSON + 1)
hpmax = np.zeros(5)
nflag = np.zeros(102, dtype=int)  # NFLAG(100); two extra so that the flags
                                  # after the gradient tests cannot fall off
lcz = np.zeros(JSON + 1, dtype=bool)
ftx = np.zeros(JSON + 1)
fpx = np.zeros(JSON + 1)
mrzone = np.zeros((14, 3), dtype=int)
mxzone = np.zeros((13, 3), dtype=int)
hfo = np.zeros(JSON + 1)
hko = np.zeros(JSON + 1)
xtab = np.zeros(JSON + 1)
ytab = np.zeros(JSON + 1)
chi = np.zeros(JSON + 1)
y2 = np.zeros(JSON + 1)
ixflag = np.zeros(JSON + 1, dtype=int)
# DATA ID: the species of HCOMPM that are interpolated (index 0 unused).
ID = np.array([0, 1, 2, 4, 5, 6, 7, 9])


def _fmt_row(i, s, p, t, r, l, x, z, o16):
    # 1X,I3,F11.7,F8.4,F8.5,F8.4,1PE9.2,0PF6.3,2F5.3
    return (f" {i:3d}{s:11.7f}{p:8.4f}{t:8.5f}{r:8.4f}{l:9.2E}"
            f"{x:6.3f}{z:5.3f}{o16:5.3f}")


def hpoint(m, hstot, hs, hs1, hs2, ht, hp, hr, lc, hl, hd, hcomp, ps, rs,
           estore, istore, lreset, jxbeg, model, lshell, jcore, jenv, omega,
           eta2, r0, hi, hjm, hkerot, sjtot, skerot, fp, ft, hg, qiw, teffl,
           cmixl):
    """Rezone the model: keep the jumps between neighbouring points in
    mass, L, P, X and Z below the user's limits (HPTTOL), keep
    discontinuities as fixed points, add or drop outer points at the
    envelope temperature limits, interpolate every variable to the new
    points, recompute the shell masses and, for rotating models, the
    rotation and the start-of-step physics.

    Translates: hpoint.f#HPOINT (09c3b9b)

    Arguments: the model (hs = log m, hp, ht, hr, hl, hd, hcomp, lc, and m
    points), hstot (log M), the envelope fit (ps, rs, estore, istore), the
    convection and shell indices (jcore, jenv, jxbeg, lshell), model, the
    rotation arrays (omega, eta2, r0, hi, hjm, hkerot, fp, ft, hg, qiw), the
    totals sjtot and skerot, teffl, and cmixl (the run's mixing length,
    passed on to physic). Arrays are changed in place, hs1 and hs2 set.
    lreset is output only.

    Returns (m, lreset, jcore, jenv, sjtot, skerot).
    """
    y = 0.0
    klo = khi = 0

    # Rezone so that the user's maximum jumps between successive points in
    # P, L, X and Z are not exceeded; flag real discontinuities.
    lreset = False
    # IEND is the number of species followed
    iend = 11
    if nmlpar_mod.lexcom:
        iend = 15

    # Outermost Henyey point colder than the minimum envelope T: drop points
    if ht[m] < nmlpar_mod.tenv0:
        for i in range(m - 1, 0, -1):  # loop 20
            if ht[i] > nmlpar_mod.tenv0:
                iounits_mod.write(iounits_mod.ishort,
                                  f" OUTER POINTS DELETED OLD M ={m:5d}"
                                  f"  NEW M ={i:5d}")
                m = i
                comp_mod.senv = hs[m] - hstot
                lreset = True
                break  # to label 40
        else:
            # The entire model has T < TENV0: stop
            msg = (" " + ">" * 39 + "<" * 40 + "\n ERROR IN HPOINT\n"
                   " MAX. STAR T LESS THAN MINIMUM ENVELOPE T.RUN STOPPED")
            iounits_mod.write(iounits_mod.ishort, msg)
            iounits_mod.write(iounits_mod.iowr, msg)
            sys.exit()
    # Outermost point hotter than the maximum envelope T: add one point
    elif ht[m] > nmlpar_mod.tenv1 and istore != 0:
        m = m + 1
        comp_mod.senv = estore[4]
        hs[m] = hstot + comp_mod.senv
        hp[m] = hp[m - 1] + (estore[1] - ps[istore])
        # Question (see the page): rough guesses for the added point's T, rho.
        ht[m] = ht[m - 1] + (estore[1] - ps[istore]) * 0.250
        hr[m] = hr[m - 1] + (estore[3] - rs[istore])
        hl[m] = hl[m - 1]
        for i in range(1, iend + 1):
            hcomp[i, m] = hcomp[i, m - 1]
        hd[m] = hp[m] - ht[m] - 8.0
        j = m - 1
        # FORMAT 60 holds 5 items: the 10 written make two records
        for k in (j, m):
            iounits_mod.write(iounits_mod.ishort,
                              f" OUTER POINT ADDED{k:5d}{hs[k]:15.10f}  PTR"
                              f"{hp[k]:10.6f}{ht[k]:10.6f}{hr[k]:10.6f}")
        lreset = True

    # Flagged points: the program does not rezone across them
    iflag = 1
    # Edges of the central and surface convection zones
    if jcore > 1:
        nflag[iflag] = jcore
        iflag = iflag + 1
    if jenv < m and jenv > 1:
        nflag[iflag] = jenv
        iflag = iflag + 1
    # Edge of the H shell
    if lshell and jxbeg > 1:
        nflag[iflag] = jxbeg - 1
        iflag = iflag + 1
    too_many = False
    for i in range(2, m + 1):  # loop 120
        if abs(hcomp[1, i] - hcomp[1, i - 1]) > nmlpar_mod.hpttol[3]:
            # X gradient
            nflag[iflag] = i
            iflag = iflag + 1
        elif abs(hcomp[3, i] - hcomp[3, i - 1]) > nmlpar_mod.hpttol[4]:
            # Z gradient
            nflag[iflag] = i
            iflag = iflag + 1
        elif nmlpar_mod.lrot:
            # gradient in log omega
            wltop = math.log10(omega[i])
            wlbot = math.log10(omega[i - 1])
            if abs(wltop - wlbot) > nmlpar_mod.hpttol[12]:
                nflag[iflag] = i
                iflag = iflag + 1
        if iflag >= 100:
            iounits_mod.write(iounits_mod.ishort,
                              " MORE THAN 100 FLAG POINTS-FIRST 100 RETAINED")
            # Question (see the page): this skips the overshoot and
            # fine-zoning flags below.
            too_many = True
            break  # GOTO 130

    if not too_many:
        # PMAX1..PMAX5: maximum steps in log P by region (computed, unused).
        # Question (see the page): PMAX1-PMAX5, JOS, JFINE are not used.
        pmax1 = nmlpar_mod.hpttol[11]
        pmax4 = nmlpar_mod.hpttol[10]
        pmax5 = nmlpar_mod.hpttol[8]
        if not lc[m]:
            jos = m
            jfine = m
        elif jenv == 1:
            jos = 1
            jfine = 1
        else:
            # Base of the overshoot region, if any
            if not nmlpar_mod.lovste:
                i = jenv
                jos = jenv
            else:
                # Question (see the page): ALPHAE used as a depth in log P.
                jos = jenv - 1
                while jos >= 1:  # loop 191
                    if hp[jos] - hp[jenv] > nmlpar_mod.alphae:
                        break  # GOTO 193
                    jos -= 1
                jos = jos + 1  # label 193
                delp = hp[jos] - hp[jenv]
                if delp > 0.0:
                    np_ = int(delp / nmlpar_mod.hpttol[10])
                    if math.fmod(delp, nmlpar_mod.hpttol[10]) != 0.0:
                        np_ = np_ + 1
                    pmax3 = delp / float(np_)
                else:
                    pmax3 = nmlpar_mod.hpttol[10]
                if jos > 1:
                    nflag[iflag] = jos
                    iflag = iflag + 1
            if nmlpar_mod.hpttol[7] == 0.0:
                jfine = jos
            elif jos == 1:
                jfine = 1  # GOTO 198
            else:
                # Base of the finely zoned region
                jfine = jos - 1
                while jfine >= 1:  # loop 195
                    if hp[jfine] - hp[jenv] > nmlpar_mod.hpttol[7]:
                        break  # GOTO 197
                    jfine -= 1
                jfine = jfine + 1  # label 197
                if jfine != jos:
                    delp = hp[jfine] - hp[jos]
                    np_ = int(delp / nmlpar_mod.hpttol[10])
                    if math.fmod(delp, nmlpar_mod.hpttol[10]) != 0.0:
                        np_ = np_ + 1
                    pmax2 = delp / float(np_)
                    if pmax2 == 0.0:
                        pmax2 = nmlpar_mod.hpttol[10]
                    if jfine > 1:
                        nflag[iflag] = jfine
                        iflag = iflag + 1
                # label 198

    nflag[iflag] = m  # label 130
    if iflag != 1:
        # Sort the flag points in ascending order (label 140)
        ldone = False
        while not ldone:
            ldone = True
            for i in range(1, iflag):
                if nflag[i + 1] < nflag[i]:
                    nflag[i], nflag[i + 1] = nflag[i + 1], nflag[i]
                    ldone = False
        # Remove duplicates (label 160). After a removal i still advances,
        # so a third equal flag in a row is kept.
        i = 2
        while True:
            if nflag[i] == nflag[i - 1]:
                if i < iflag:
                    for j in range(i, iflag):
                        nflag[j] = nflag[j + 1]
                iflag = iflag - 1
            i = i + 1
            if not i <= iflag:
                break
    # label 180: FORMAT 185 writes 20 flags per record
    flags = [nflag[j] for j in range(1, iflag + 1)]
    for k in range(0, max(len(flags), 1), 20):
        iounits_mod.write(iounits_mod.ishort, " FLAG-POINTS"
                          + "".join(f"{n:4d}" for n in flags[k:k + 20]))

    # Refloat the points
    if nmlpar_mod.lrot:
        for i in range(1, m + 1):  # loop 190
            if omega[i] > 0.0:
                omegal[i] = math.log10(omega[i])
            else:
                omegal[i] = 0.0
    oldmod_mod.hss[1] = hs[1]
    oldmod_mod.hpo[1] = hp[1]
    oldmod_mod.hlo[1] = hl[1]
    hko[1] = hcomp[1, 1]
    hfo[1] = hcomp[3, 1]
    hlmax = hl[m]
    # Teff saved, and the number of points saved (JVS 04/14, 05/25)
    oldmod_mod.tefflo = teffl
    # Question (see the page): MO is the old number of points.
    oldmod_mod.mo = m
    for i in range(m - 1, 0, -1):
        if hl[i] > hlmax:
            hlmax = hl[i]
    hpmax[1] = nmlpar_mod.hpttol[8]
    hpmax[2] = nmlpar_mod.hpttol[9] * hlmax
    hpmax[3] = nmlpar_mod.hpttol[5]
    hpmax[4] = nmlpar_mod.hpttol[6]
    iiflg = 1
    j = 2
    i0 = 2
    # CHI: the running total of normalized differences in m, L, P
    # (dCHI = dHS/DM + dHL/DL - dHP/DP, the L term only where L increases)
    dm = nmlpar_mod.hpttol[2]
    dl = hpmax[2]
    chi[1] = 1.0
    for j in range(2, m + 1):
        testp = hp[j] - hp[jenv]
        if abs(testp) < nmlpar_mod.hpttol[7]:
            dp = nmlpar_mod.hpttol[10]   # finely zoned region
        elif testp > nmlpar_mod.hpttol[7]:
            dp = nmlpar_mod.hpttol[11]   # below the surface convection zone
        else:
            dp = nmlpar_mod.hpttol[8]    # in the surface convection zone
        if hl[j] > hl[j - 1]:
            dchi = ((hs[j] - hs[j - 1]) / dm + (hl[j] - hl[j - 1]) / dl
                    - (hp[j] - hp[j - 1]) / dp)
        else:
            dchi = (hs[j] - hs[j - 1]) / dm - (hp[j] - hp[j - 1]) / dp
        chi[j] = chi[j - 1] + dchi
    for j in range(1, m + 1):
        xtab[j] = chi[j]
        ytab[j] = hs[j]
    # Spline of log m against CHI
    splinc(xtab, ytab, y2, m)
    chiprev = chi[1]
    jprev = 1
    mnew = 1
    for i in range(1, iflag + 1):
        delchi = chi[nflag[i]] - chiprev
        npt = int(delchi) + 1
        if mnew + npt > JSON:
            # Exceeded the array dimensions: complain and stop
            print(f" DESIRED NUMBER OF POINTS {mnew + npt:6d}EXCEEEDS JSON "
                  f"{JSON:6d}\n RUN STOPPED")
            sys.exit(9999)
        dchi = delchi / float(npt)
        # Assign the new points, evenly in CHI
        for j in range(jprev + 1, jprev + npt + 1):
            x = chiprev + dchi
            y, klo, khi = splintd2(xtab, ytab, m, y2, x, y, klo, khi)
            oldmod_mod.hss[j] = y
            chiprev = x
        mnew = mnew + npt
        jprev = mnew

    # Points from the gradient in X: X at the new points, into HKO
    for j in range(1, m + 1):
        xtab[j] = hs[j]
        ytab[j] = hcomp[1, j]
    splinc(xtab, ytab, y2, m)
    for i in range(2, mnew + 1):
        x = oldmod_mod.hss[i]
        y, klo, khi = splintd2(xtab, ytab, m, y2, x, y, klo, khi)
        hko[i] = y

    # Test on the change in X (only where X increases outward); skipped if
    # HPMAX(3) is zeroed out
    if hpmax[3] > 1.0e-15:
        ixx = 0
        for j in range(mnew, 1, -1):
            if hko[j] - hko[j - 1] > hpmax[3]:
                ixx = ixx + 1
                ixflag[ixx] = j
        if ixx != 0:
            mm = mnew
            # Question (see the page): no check that mm stays within JSON.
            for i in range(1, ixx + 1):
                j = ixflag[i]
                dx = (hko[j] - hko[j - 1]) / hpmax[3]
                nn = int(dx)  # the number of new points needed
                dss = ((oldmod_mod.hss[j] - oldmod_mod.hss[j - 1])
                       / float(nn + 1))
                for k in range(mm + nn, j + nn - 1, -1):
                    oldmod_mod.hss[k] = oldmod_mod.hss[k - nn]
                for k in range(j + nn - 1, j - 2, -1):
                    oldmod_mod.hss[k] = oldmod_mod.hss[k + 1] - dss
                mm = mm + nn
            mnew = mm
    # label 102

    # Points from the gradient in Z: Z at the new points, into HFO
    for j in range(1, m + 1):
        xtab[j] = hs[j]
        ytab[j] = hcomp[3, j]
    splinc(xtab, ytab, y2, m)
    for i in range(2, mnew + 1):
        x = oldmod_mod.hss[i]
        y, klo, khi = splintd2(xtab, ytab, m, y2, x, y, klo, khi)
        hfo[i] = y
    # Test on the change in Z (only where Z decreases outward); skipped if
    # HPMAX(4) is zeroed out
    if hpmax[4] > 1.0e-15:
        ixx = 0
        for j in range(mnew, 1, -1):
            if hfo[j - 1] - hfo[j] > hpmax[4]:
                ixx = ixx + 1
                ixflag[ixx] = j
        if ixx != 0:
            mm = mnew
            for i in range(1, ixx + 1):
                j = ixflag[i]
                dz = (hfo[j - 1] - hfo[j]) / hpmax[4]
                nn = int(dz)  # the number of new points needed
                dss = ((oldmod_mod.hss[j] - oldmod_mod.hss[j - 1])
                       / float(nn + 1))
                for k in range(mm + nn, j + nn - 1, -1):
                    oldmod_mod.hss[k] = oldmod_mod.hss[k - nn]
                for k in range(j + nn - 1, j - 2, -1):
                    oldmod_mod.hss[k] = oldmod_mod.hss[k + 1] - dss
                mm = mm + nn
            mnew = mm
    # label 103

    # Delete new points that are too close together (HDO as scratch)
    j = 1
    oldmod_mod.hdo[j] = oldmod_mod.hss[j]
    for jj in range(2, mnew):  # loop 810
        if oldmod_mod.hss[jj] - oldmod_mod.hdo[j] > nmlpar_mod.hpttol[1]:
            j = j + 1
            oldmod_mod.hdo[j] = oldmod_mod.hss[jj]
    j = j + 1
    oldmod_mod.hdo[j] = oldmod_mod.hss[mnew]
    mnew = j
    for j in range(2, mnew + 1):  # loop 820
        oldmod_mod.hss[j] = oldmod_mod.hdo[j]

    # Outer edge of the convective core and inner edge of the convective
    # envelope on the new points
    if jcore > 1:
        j = 2
        while j <= mnew:  # loop 823
            if oldmod_mod.hss[j] > hs[jcore]:
                break  # GOTO 824
            j += 1
        jcore = j - 1  # label 824
    else:
        jcore = 1
    if jenv < m:
        j = mnew - 1
        while j >= 1:  # loop 825
            if oldmod_mod.hss[j] < hs[jenv]:
                break  # GOTO 826
            j -= 1
        jenv = j + 1  # label 826
    else:
        jenv = mnew
    if jcore > 1:
        for i in range(1, jcore + 1):
            lc[i] = True
        jbeg = jcore + 1
    else:
        jbeg = 1
    if jenv < mnew:
        for i in range(jenv, mnew + 1):
            lc[i] = True
        jend = jenv - 1
    else:
        jend = mnew
    if not (jend < 1 or jbeg > mnew):
        for j in range(jbeg, jend + 1):  # loop 829
            lc[j] = False
    # label 830

    # Oscillatory spline of the old run of each variable (on HS, NTAB = M
    # points) at the new points (HSS, NTOT = MNEW): osplin(xval, yval,
    # xtab, ytab, ntab, ntot). HPO and HTO serve as scratch.
    ntab = m
    ntot = mnew

    def interpolate(arr):
        # copy the old run into HPO, interpolate into HTO, copy back
        for j in range(1, m + 1):
            oldmod_mod.hpo[j] = arr[j]
        osplin(oldmod_mod.hss, oldmod_mod.hto, hs, oldmod_mod.hpo, ntab, ntot)
        for j in range(1, mnew + 1):
            arr[j] = oldmod_mod.hto[j]

    # The entropy terms (loops 904-914) (7/91)
    interpolate(entrop_mod.htt0)
    interpolate(entrop_mod.hpp0)
    interpolate(entrop_mod.hll0)
    interpolate(entrop_mod.hrr0)

    # The composition, and HCOMPP, the composition at the start of the step
    # (needed for composition diffusion in rotating models): loop 850
    for i in range(1, iend + 1):
        interpolate(hcomp[i])
        interpolate(oldmod_mod.hcompp[i])

    # HCOMPM, the changes in composition from nuclear burning (loop 849)
    for i in range(1, 8):
        interpolate(burn_mod.hcompm[ID[i]])

    # The deuterium burning rates, if surface D >= 1e-14 (MHP 05/02)
    if nmlpar_mod.lexcom and hcomp[12, m] >= 1.0e-14:
        interpolate(deuter_mod.drate0)

    # P, R, L, T and rho at the new points
    osplin(oldmod_mod.hss, oldmod_mod.hpo, hs, hp, ntab, ntot)
    osplin(oldmod_mod.hss, oldmod_mod.hro, hs, hr, ntab, ntot)
    osplin(oldmod_mod.hss, oldmod_mod.hlo, hs, hl, ntab, ntot)
    osplin(oldmod_mod.hss, oldmod_mod.hto, hs, ht, ntab, ntot)
    osplin(oldmod_mod.hss, oldmod_mod.hdo, hs, hd, ntab, ntot)

    # Rotating models: omega, j/m, fp, ft, r0 and eta2
    if nmlpar_mod.lrot:
        osplin(oldmod_mod.hss, oldrot_mod.wold, hs, omega, ntab, ntot)
        osplin(oldmod_mod.hss, oldrot_mod.hjx, hs, hjm, ntab, ntot)
        osplin(oldmod_mod.hss, fpx, hs, fp, ntab, ntot)
        osplin(oldmod_mod.hss, ftx, hs, ft, ntab, ntot)
        osplin(oldmod_mod.hss, oldrot_mod.r0x, hs, r0, ntab, ntot)
        osplin(oldmod_mod.hss, oldrot_mod.eta2x, hs, eta2, ntab, ntot)

    # The point distribution, if requested
    if nmlpar_mod.ldebug and nmlpar_mod.npoint < 9999:
        if model % nmlpar_mod.npoint == 0:
            idebug = iounits_mod.idebug
            minmin = min(m, mnew)
            # FORMAT 910 ('1' starts a new page)
            iounits_mod.write(idebug, "1" + " " * 20 + "OLD POINTS" + " " * 54
                              + "NEW POINTS\n"
                              + 2 * "   N     S        P       T       R"
                                    "        L       X    Z   O16 ")
            # Question (see the page): the "old points" columns X and Z are
            # HKO and HFO, X and Z at the new points.
            for i in range(1, minmin + 1):
                iounits_mod.write(idebug, _fmt_row(
                    i, hs[i], hp[i], ht[i], hr[i], hl[i], hko[i], hfo[i],
                    hcomp[9, i]) + _fmt_row(
                    i, oldmod_mod.hss[i], oldmod_mod.hpo[i],
                    oldmod_mod.hto[i], oldmod_mod.hro[i], oldmod_mod.hlo[i],
                    hcomp[1, i],
                    hcomp[3, i], hcomp[9, i]))
            if m > minmin:
                minmin = minmin + 1
                for i in range(minmin, m + 1):
                    iounits_mod.write(idebug, _fmt_row(
                        i, hs[i], hp[i], ht[i], hr[i], hl[i], hko[i],
                        hfo[i], hcomp[9, i]))
            elif mnew > minmin:
                minmin = minmin + 1
                for i in range(minmin, mnew + 1):
                    iounits_mod.write(idebug, " " * 64 + _fmt_row(
                        i, oldmod_mod.hss[i], oldmod_mod.hpo[i],
                        oldmod_mod.hto[i], oldmod_mod.hro[i],
                        oldmod_mod.hlo[i],
                        hcomp[1, i], hcomp[3, i], hcomp[9, i]))

    # Transfer the new points (loop 1000)
    for j in range(1, mnew + 1):
        hs[j] = oldmod_mod.hss[j]
        hp[j] = oldmod_mod.hpo[j]
        ht[j] = oldmod_mod.hto[j]
        hr[j] = oldmod_mod.hro[j]
        hl[j] = oldmod_mod.hlo[j]
        hd[j] = oldmod_mod.hdo[j]
    if nmlpar_mod.lrot:
        for j in range(1, mnew + 1):  # loop 1005
            hjm[j] = oldrot_mod.hjx[j]
            omega[j] = oldrot_mod.wold[j]
            fp[j] = fpx[j]
            ft[j] = ftx[j]
            eta2[j] = oldrot_mod.eta2x[j]
            r0[j] = oldrot_mod.r0x[j]
    # The energy generation at the start of the step (MHP 6/00). Note that HS
    # already holds the new points here (loop 1000), while SESUM and SEG are on
    # the old ones.
    if nmlpar_mod.lrot or (nmlpar_mod.lexcom and nmlpar_mod.lovste):
        osplin(oldmod_mod.hss, oldphy_mod.esumo, hs, scrtch_mod.sesum, ntab,
               ntot)
        for im in range(1, m + 1):
            ytab[im] = (scrtch_mod.sesum[im] + scrtch_mod.seg[6, im]
                        + scrtch_mod.seg[7, im])
        osplin(oldmod_mod.hss, oldphy_mod.epso, hs, ytab, ntab, ntot)
    iounits_mod.write(iounits_mod.ishort, f" POINTS  OLD{m:5d}   NEW{mnew:5d}")
    m = mnew

    # Masses: HS1 = m at each point, HS2 = the mass of each shell
    clndp = constants_mod.clndp
    ds3 = math.exp(clndp * hs[1])
    ds2 = -ds3
    for i in range(2, m + 1):  # loop 1030
        ds1 = ds2
        ds2 = ds3
        ds3 = math.exp(clndp * hs[i])
        hs1[i - 1] = ds2
        hs2[i - 1] = 0.5 * (ds3 - ds1)
    hs1[m] = ds3
    hs2[m] = math.exp(constants_mod.cln * hstot) - 0.5 * (ds2 + ds3)

    if nmlpar_mod.lrot:
        # First guess of the moment of inertia
        for i in range(1, m + 1):
            hi[i] = constants_mod.cc23 * hs2[i] * math.exp(
                constants_mod.cln * 2.0 * hr[i])
        # Overshoot
        nrzone, nzone = ovrot(hcomp, hd, hp, hr, hs, ht, lc, m, lcz, mrzone,
                              mxzone, 0, 0)
        # The moment of inertia again, on the new points, before the rotation (JNT 2025/09/03).
        # Question (see the page): is this extra call still needed?
        momi(eta2, hr, hs, hs2, 1, m, omega, r0, hi, qiw)
        # Omega from the specific angular momentum HJM
        getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
        # fp, ft, r0 and eta2 from omega
        fpft(hd, hr, hs, m, omega, eta2, fp, ft, hg, r0)
        # The moment of inertia
        momi(eta2, hr, hs, hs2, 1, m, omega, r0, hi, qiw)
        # The new total angular momentum and rotational kinetic energy
        sumj = 0.0
        sumke = 0.0
        for i in range(1, m + 1):  # loop 1110
            hj = hjm[i] * hs2[i]
            hkerot[i] = 0.5 * omega[i] * hj
            sumj = sumj + hj
            sumke = sumke + hkerot[i]
        iounits_mod.write(iounits_mod.ishort,
                          f" TOTAL J OF STAR - PREVIOUS {sjtot:21.13E} NEW "
                          f"{sumj:21.13E}\n TOTAL ROTATIONAL K.E. OF STAR-"
                          f"PREVIOUS {skerot:21.13E} NEW {sumke:21.13E}")
        sjtot = sumj
        skerot = sumke
        # Store the old model structure for use in diffusion (loop 1130)
        for i in range(1, m + 1):
            oldrot_mod.wold[i] = omega[i]
            oldrot_mod.hgo[i] = hg[i]
            oldrot_mod.hio[i] = hi[i]
            oldrot_mod.eta2x[i] = eta2[i]
            oldrot_mod.r0x[i] = r0[i]
            oldmod_mod.lco[i] = lc[i]
            oldmod_mod.lczo[i] = lcz[i]
            oldrot_mod.hjx[i] = hjm[i]

    # The physics at the new points: for rotation, or for light elements
    # with envelope overshoot (to locate the base of the overshoot region)
    if nmlpar_mod.lrot or (nmlpar_mod.lexcom and nmlpar_mod.lovste):
        physic(fp, ft, hcomp, hd, hg, hl, hp, hr, hs, ht, lc, m, teffl, cmixl)
        # Store what diffusion needs at the start of the step (loop 1040)
        for im in range(1, m + 1):
            oldphy_mod.delro[im] = scrtch_mod.sdel[1, im]
            oldphy_mod.delo[im] = scrtch_mod.sdel[2, im]
            oldphy_mod.delao[im] = scrtch_mod.sdel[3, im]
            oldphy_mod.amuo[im] = temp_mod.hamu[im]
            oldphy_mod.oo[im] = scrtch_mod.so[im]
            oldphy_mod.cpo[im] = temp_mod.cp[im]
            oldphy_mod.qdto[im] = temp_mod.sqdt[im]
            oldphy_mod.velo[im] = scrtch_mod.svel[im]
            oldphy_mod.visco[im] = temp_mod.visc[im]
            oldphy_mod.thdifo[im] = temp_mod.thdif[im]
            oldrot2_mod.ddelm[im] = oldphy_mod.delao[im] - oldphy_mod.delo[im]
            dwmax_mod.qwrmx0[im] = dwmax_mod.qwrmax[im]
        # The time evolution of the angular velocity (MHP 06/02)
        cln = constants_mod.cln
        for i in range(2, m + 1):
            dr = math.exp(cln * hr[i]) - math.exp(cln * hr[i - 1])
            dw = omega[i] - omega[i - 1]
            wmid = 0.5 * (omega[i] + omega[i - 1])
            fac = (2.0 * (hr[i] + hr[i - 1]) - 0.5 * (hs[i] + hs[i - 1])
                   - constants_mod.cgl)
            oldrot2_mod.tho[i] = math.exp(cln * fac) * wmid * dw / dr
            oldrot2_mod.qwrst[i] = dw / dr

    # A new surface opacity table, if needed
    if abs(comp_mod.xnew - hcomp[1, m]) > 1.0e-8:
        comp_mod.xnew = hcomp[1, m]
        comp_mod.znew = hcomp[3, m]
        surfopac(comp_mod.xnew)

    return m, lreset, jcore, jenv, sjtot, skerot
