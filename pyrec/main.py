"""Python version of main.f (YREC a8a5950): YREC's main program, which reads
the run's settings, sets up the tables, and evolves the starting model of
each kind card model by model."""

import math

import numpy as np

from pyrec import acdpth_mod
from pyrec import cent_mod
from pyrec import cenv_mod
from pyrec import constants_mod
from pyrec import ctol_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod
from pyrec import disk_mod
from pyrec import entrop_mod
from pyrec import fluxes_mod
from pyrec import iounits_mod
from pyrec import label_mod
from pyrec import liov_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import ovrtrn_mod
from pyrec import prevcz_mod
from pyrec import runloop_mod
from pyrec import sound_mod
from pyrec.chkcal import chkcal
from pyrec.chkscal import chkscal
from pyrec.convec import convec
from pyrec.crrect import crrect
from pyrec.engeb import engeb
from pyrec.findsh import findsh
from pyrec.fpft import fpft
from pyrec.getnewenv import getnewenv
from pyrec.getw import getw
from pyrec.hpoint import hpoint
from pyrec.htimer import htimer
from pyrec.liburn import liburn
from pyrec.lirate88 import lirate88
from pyrec.massloss import massloss
from pyrec.mix import mix
from pyrec.mixcz import mixcz
from pyrec.parmin import parmin
from pyrec.pdist import pdist
from pyrec.putstore import putstore
from pyrec.setcal import setcal
from pyrec.setscal import setscal
from pyrec.setups import setups
from pyrec.starin import starin
from pyrec.version import setversion
from pyrec.wrthead import wrthead
from pyrec.wrtlst import wrtlst
from pyrec.wrtout import wrtout

JSON = dimensions_mod.json


def main():
    """YREC's main program: settings and tables, then for each kind card
    the starting model and the evolution loop (mass loss, mixing and
    burning, rezoning, four levels of CRRECT, rotation, light elements,
    time step, output), with the solar and stellar calibrations around it.

    Translates: main.f#(file level) (a8a5950)

    No arguments; it ends with STOP (SystemExit).
    """
    c = constants_mod
    nml = nmlpar_mod
    io = iounits_mod
    rl_ = runloop_mod
    ent = entrop_mod
    old = oldmod_mod

    # The model and its work arrays (index 0 unused).
    hs, hl, hr, hp, ht, hd = (np.zeros(JSON + 1) for _ in range(6))
    lc = np.zeros(JSON + 1, dtype=bool)
    hcomp = np.zeros((16, JSON + 1))
    hs1, hs2, hhc = (np.zeros(JSON + 1) for _ in range(3))
    hha = np.zeros((5, 3, JSON + 1))
    hhb = np.zeros((5, JSON + 1))
    hpp, htt = np.zeros(JSON + 1), np.zeros(JSON + 1)
    fp, ft, omega, hi, hjm = (np.zeros(JSON + 1) for _ in range(5))
    hkerot, r0, eta2, qiw, hg = (np.zeros(JSON + 1) for _ in range(5))
    hkeroto = np.zeros(JSON + 1)
    tlumx, tril, trit = np.zeros(9), np.zeros(4), np.zeros(4)
    ps, ts, rs, cfenv = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(10)
    estore, hmax = np.zeros(5), np.zeros(5)
    nhmax = np.zeros(5, dtype=int)
    hsrf, v = np.zeros(7), np.zeros(13)
    # Start-of-step j and composition when iterating structure and rotation
    # (formerly COMMON/ORIGSTART/, now local arrays).
    hjm0 = np.zeros(JSON + 1)
    hcompp0 = np.zeros((16, JSON + 1))
    # The neutrino table (off).
    be7m = np.zeros(JSON + 1)
    fluxi0 = np.zeros((11, JSON + 1))
    hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12, hr13 = (
        np.zeros(JSON + 1) for _ in range(13))
    hf1, hf2 = np.zeros(JSON + 1), np.zeros(JSON + 1)
    mxzone = np.zeros((13, 3), dtype=int)
    mxzon0 = np.zeros((13, 3), dtype=int)
    mrzone = np.zeros((14, 3), dtype=int)
    lcz = np.zeros(JSON + 1, dtype=bool)

    # Scalars passed before any routine sets them, or set by the callees
    # (GL by WRTOUT; LATMO, LDERIV, LMIX, LOCOND by CRRECT; QDT, QDP by
    # STARIN). Fortran leaves them undefined; 0 or False here.
    bl = dage = ddage = delts = deltsh = ftri = hstot = smass = teffl = 0.0
    dwmax = gl = qdt = qdp = sjtot = skerot = senvold = vel = 0.0
    ikut = istore = jcore = jenv = jxbeg = jxend = jxmid = m = model = 0
    nrzone = nzone = nzone0 = itdone = jenv0 = 0
    latmo = lderiv = lmix = locond = lnew = lshell = lnewfit = False
    lprt0 = False
    compmix = ''
    # Question (see the page): IU is never set, but passed to ENGEB.
    iu = 0
    # DATA DLDX,DRDX,DLDA,DRDA,LSOL/4*0.0D0,.FALSE./ (explicit, since
    # SETCAL may never run while LSOL is still tested at the end).
    dldx = drdx = dlda = drda = 0.0
    lsol = False
    # COMMON/CALSTAR/ is now these three locals, set by SETSCAL, CHKSCAL.
    alri = 0.0
    lstar = lpassr = False
    # DATA NAO/1/: the next of the ages AGEOUT.
    nao = 1

    def _crrect(niter, lsbc, itlvl):
        """The four calls to CRRECT (below) differ only in NITER,
        LSBC and ITLVL; the other arguments are main's current values."""
        return crrect(
            delts, m, niter, lfini, large, lnew, lreset, lsbc,
            trit, tril, cfenv, ps, ts, rs, ftri, istore, estore,
            hstot, bl, teffl, hd, hha, hhb, hhc, hl, hmax, hp,
            hpp, hr, hs, hs1, hs2, hsrf, ht, htt, hcomp, lc,
            nhmax, tlumx, latmo, lderiv, lmix, locond, qdt, qdp,
            fp, ft, eta2, omega, r0, itdone, hg, hi, hjm, itlvl,
            lcz, mxzone, qiw, hkerot, hkeroto, cmixl)

    # Start.
    yrecver, githash = setversion('', '')
    io.iowr = 9
    # LPUNCH is true once the first model is calculated.
    lpunch = False
    # Read the run's parameters; it returns the table files' names.
    (falex06, fallard, fatm, ffermi, fliv95, fpatm, fpenv, fpmod, fpurez,
     fscvh, fscvhe, fscvz) = parmin('', '', '', '', '', '', '', '', '', '',
                                    '', '', yrecver, githash)
    # Constants and tables.
    setups(falex06, fallard, fatm, ffermi, fliv95, fpurez, fscvh, fscvhe,
           fscvz)
    # The Monte Carlo option is removed: one pass, the age factor 1.
    facage = 1.0
    lsavpu = nml.lpulse
    # Automatic solar calibration: starting estimates of alpha and X.
    if nml.lcals:
        lsol = setcal(facage, lsol)
        iconv = 1
        xguess = rl_.rescal[2, 1]
        aguess = nml.cmixla[1]
        lptimesv = nml.lptime   # LPTIME, KTTAU restored at each cycle
        kttausv = nml.kttau
    else:
        iconv = 0
    # Calibration of a star of given Teff and L.
    if nml.lcalst:
        alri, lstar, lpassr = setscal(alri, lstar, lpassr)

    # The kind-card loop (DO 200).
    for nk in range(1, nml.numrun + 1):
        sound_mod.lsound = False
        label_mod.xenv0 = nml.xenv0a[nk]
        label_mod.zenv0 = nml.zenv0a[nk]
        cmixl = nml.cmixla[nk]
        cenv_mod.lsenv0 = nml.lsenv0a[nk]
        cenv_mod.senv0 = nml.senv0a[nk]
        lreset = False
        large = False
        sjtot = 0.0
        skerot = 0.0
        # Read in (or keep) the starting model.
        (bl, dage, ddage, delts, deltsh, ftri, hstot, ikut, istore, lnew,
         m, model, qdp, qdt, sjtot, skerot, smass, teffl, vel,
         compmix) = starin(
            bl, cfenv, dage, ddage, delts, deltsh, eta2, fp, ft, ftri,
            hcomp, hd, hi, hjm, hkerot, hl, hp, hr, hs, hs1, hs2, hstot, ht,
            ikut, istore, jenv, large, lc, lnew, m, model, nk, omega, ps,
            qdp, qdt, qiw, r0, rs, sjtot, skerot, smass, teffl, tlumx, tril,
            trit, ts, vel, hg, v, cmixl, compmix)
        if omega[1] == 0 and nml.lrot:
            print('LROT set to TRUE, but OMEGA(1) = 0. Stopping.'
                  ' Initialize rotation rates or set LROT to FALSE.')
            raise SystemExit
        # Disable stop conditions the starting model already meets.
        if rl_.lendag[nk]:
            if nml.end_dcen[nk] > 0.0 and hcomp[12, 1] < nml.end_dcen[nk]:
                nml.end_dcen[nk] = -nml.end_dcen[nk]
                msg = (f"STARTING D {hcomp[12, 1]:12.4E} BELOW STOP VALUE "
                       f"{nml.end_dcen[nk]:12.4E} STOP DISABLED.")  # 101
                print(msg)
                io.write(io.ishort, msg)
            if nml.end_xcen[nk] > 0.0 and hcomp[1, 1] < nml.end_xcen[nk]:
                # Question (see the page): END_DCEN set, and D printed.
                nml.end_dcen[nk] = -nml.end_xcen[nk]
                msg = (f"STARTING X {hcomp[12, 1]:12.4E} BELOW STOP VALUE "
                       f"{nml.end_dcen[nk]:12.4E} STOP DISABLED.")  # 102
                print(msg)
                io.write(io.ishort, msg)
            if nml.end_ycen[nk] > 0.0 and hcomp[2, 1] < nml.end_ycen[nk]:
                nml.end_ycen[nk] = -nml.end_ycen[nk]
                # Question (see the page): prints D and END_DCEN.
                msg = (f"STARTING Y {hcomp[12, 1]:12.4E} BELOW STOP VALUE "
                       f"{nml.end_dcen[nk]:12.4E} STOP DISABLED.")  # 103
                print(msg)
                io.write(io.ishort, msg)

        # The neutrino table, switched off (MHP 2/04).
        # Question (see the page): remove the block, or repair it?
        lnutab = False
        if lnutab:
            # Weights and masses: HS1, the shell centres (g, unlogged);
            # HS2, the shells' masses (g).
            ds3 = math.exp(c.cln * hs[1])
            ds2 = -ds3
            for i in range(2, m + 1):
                ds1 = ds2
                ds2 = ds3
                ds3 = math.exp(c.cln * hs[i])
                hs1[i - 1] = ds2
                hs2[i - 1] = 0.5 * (ds3 - ds1)
            hs1[m] = ds3
            hs2[m] = math.exp(c.cln * hstot) - 0.5 * (ds2 + ds3)
            pdt = -1.0
            pdp = 1.0
            for j in range(1, 11):
                fluxes_mod.fluxtot[j] = 0.0
                for k in range(1, m + 1):
                    fluxi0[j, k] = 0.0
            for i in range(1, m + 1):
                dl = hd[i]
                tl = ht[i]
                if tl < 6.0:
                    break                                 # go to 666
                x = hcomp[1, i]
                y = hcomp[2, i]
                z = hcomp[3, i]
                xhe3 = hcomp[4, i]
                xc12 = hcomp[5, i]
                xc13 = hcomp[6, i]
                xn14 = hcomp[7, i]
                xn15 = hcomp[8, i]
                xo16 = hcomp[9, i]
                xo17 = hcomp[10, i]
                xo18 = hcomp[11, i]
                xh2 = 0.0  # Fortran: undefined here
                (epp1, epp2, epp3, ecn, e3al, pep, pet, sum1,
                 be7massfraction, ealpca, enu) = engeb(
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dl, tl, x, y,
                    xhe3, xc12, xc13, xn14, xo16, xo18, xh2, iu, hr1, hr2,
                    hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12,
                    hr13, hf1, hf2, 0.0, 0.0, 0.0)
                be7m[i] = be7massfraction
                # From erg/g/s to erg/s for each shell.
                for j in range(1, 11):
                    fluxi0[j, i] = fluxes_mod.flux[j] * hs2[i]
                    fluxes_mod.fluxtot[j] = (fluxes_mod.fluxtot[j]
                                             + fluxi0[j, i])
                print(f"{i:5d}" + f"{hs2[i]:10.3E}"
                      + "".join(f"{fluxi0[j, i]:10.3E}"
                                for j in range(1, 11)))  # format 911
            # label 666
            io.write(76, "".join(f"{fluxes_mod.fluxtot[i]:10.3E}"
                                 for i in range(1, 11)))  # format 222
            for j in range(1, 11):
                for i in range(1, m + 1):
                    fluxi0[j, i] = fluxi0[j, i] / fluxes_mod.fluxtot[j]
            for i in range(1, m + 1):
                t6 = math.exp(c.cln * (ht[i] - 6.0))
                if t6 < 5.0:
                    break                                 # go to 141
                rhoe = hd[i] + math.log10((1.0 + hcomp[1, i]) / 2.0)
                fm = hs2[i] / 1.9891e33
                fr = math.exp(c.cln * hr[i]) / c.crsun
                io.write(76, f"{fr:9.5f}{t6:7.3f}{rhoe:6.3f}"
                         + "".join(f"{q:10.3E}" for q in (
                             fm, be7m[i], fluxi0[1, i], fluxi0[5, i],
                             fluxi0[6, i], fluxi0[7, i], fluxi0[8, i],
                             fluxi0[4, i], fluxi0[2, i],
                             fluxi0[3, i])))  # format 145
            # label 141

        # Mass in solar units, kept for the output (formerly COMMON/PULSE/).
        xmsol = smass
        # Store the starting convection-zone properties (MHP 08/02).
        deuter_mod.jcz = jenv
        ovrtrn_mod.taucz = 0.0
        # Headers of the output files.
        wrthead(smass, cmixl, nk)
        # For deciding when to dump pulsation output.
        pol1 = bl
        pot1 = teffl
        poa1 = dage
        polen = 0.0
        if nml.lkuthe:
            # Time-step cutting needs a model in unit ILAST, or it will
            # crash: copy the initial model there.
            # Question (see the page): only if LPUNCH, and before FINDSH.
            if lpunch:
                iwrite = io.ilast
                wrtlst(iwrite, hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril,
                       ps, ts, rs, cfenv, ftri, tlumx, jcore, jenv, model, m,
                       smass, teffl, bl, hstot, dage, ddage, omega, cmixl,
                       compmix)
        # The H-burning shell and the convection-zone edges.
        jcore, jenv, jxbeg, jxend, jxmid, lshell = findsh(
            hcomp, hl, lc, m, jcore, jenv, jxbeg, jxend, jxmid, lshell)
        # The first time step.
        delts, deltsh, ddage, dwmax = htimer(
            delts, deltsh, m, hd, hl, hs1, hs2, ht, hcomp, jcore, jxmid,
            tlumx, dage, ddage, nk, hp, hr, omega, dwmax, jxbeg, teffl)
        delts0 = delts
        # Zero the rates of change for the first guess.
        for i in range(1, m + 1):
            ent.htt0[i] = 0.0
            ent.hpp0[i] = 0.0
            ent.hll0[i] = 0.0
            ent.hrr0[i] = 0.0
        # Light-element rates of the surface zone.
        if nml.lexcom:
            prevcz_mod.bli6old = 0.0
            prevcz_mod.bli7old = 0.0
            prevcz_mod.bbe9old = 0.0

        # The model loop (DO 100): one evolutionary step each.
        goto200 = False
        for modeln in range(1, nml.nmodls[nk] + 1):
            if nml.lrwsh:
                io.rewind(io.ishort)
            # Pulsation output on the last model of the run.
            if modeln == nml.nmodls[nk] and nk == nml.numrun:
                nml.lpulse = lsavpu
            # Pulsation output at chosen ages (AGEOUT), for the steps before
            # and after each (JVS 02/11).
            a = acdpth_mod
            if a.lacout:
                a.ladon = True
            else:
                a.ladon = False
            if a.lacout:
                if a.ljwrt:
                    print('LJWRT on')
                    nml.lpulse = lsavpu
                    nao = nao + 1
                    a.lclcd = True
                    a.ljlast = False
                    a.ljwrt = False
                elif not a.ljwrt:
                    a.lclcd = False
                if nao < 6:
                    if (dage + ddage / 1.0e9 - a.ageout[nao] <= 0.0
                            and dage + 2.0 * ddage / 1.0e9
                            - a.ageout[nao] >= 0.0
                            and not a.ljwrt):
                        print('AGEOUT reached')
                        nml.lpulse = lsavpu
                        a.lclcd = True
                        a.ljlast = True
                        a.ljwrt = True
            # At the final age, pulsation output and sound speed (MHP 10/24).
            if (rl_.lendag[nk] and nml.endage[nk] > 0.0
                    and abs(nml.endage[nk] - dage * 1.0e9 - ddage) <= 1.0):
                nml.lpulse = lsavpu
                sound_mod.lsound = True
            if nml.lpout:
                pol1, pot1, poa1, polen = pdist(
                    pol1, pot1, poa1, polen, bl, teffl, modeln, nk, fpatm,
                    fpenv, fpmod, dage)

            # Label 15: the step, from STARIN (after a failure) to the end
            # of the rotation iterations. GOTO 15 is `continue` below.
            while True:
                # "STARIN called here for timestep cutting": with LARGE
                # set, STARIN prints RUN STOPPED and stops.
                if large:
                    (bl, dage, ddage, delts, deltsh, ftri, hstot, ikut,
                     istore, lnew, m, model, qdp, qdt, sjtot, skerot,
                     smass, teffl, vel, compmix) = starin(
                        bl, cfenv, dage, ddage, delts, deltsh, eta2, fp, ft,
                        ftri, hcomp, hd, hi, hjm, hkerot, hl, hp, hr, hs,
                        hs1, hs2, hstot, ht, ikut, istore, jenv, large, lc,
                        lnew, m, model, nk, omega, ps, qdp, qdt, qiw, r0,
                        rs, sjtot, skerot, smass, teffl, tlumx, tril, trit,
                        ts, vel, hg, v, cmixl, compmix)
                    if omega[1] == 0 and nml.lrot:
                        print('LROT set to TRUE, but OMEGA(1) = 0. Stopping.'
                              ' Initialize rotation rates or set LROT to'
                              ' FALSE.')  # format 18
                        raise SystemExit
                lpunch = True
                # Age the model, unless the kind card only relaxes a
                # rescaled model (pre-main-sequence models, log Tc < 6.6,
                # still age as a failsafe, unless L_ZERODT_PREMS is set).
                if rl_.iresca[nk] != 2 and model >= 0:
                    lage = True
                elif (not nml.l_zerodt_prems and model >= 0
                      and ht[1] < 6.6):
                    lage = True
                else:
                    lage = False
                lnewfit = False
                if lage:
                    # Mass loss.
                    delts, hstot, lnew, smass, senvold, lnewfit = massloss(
                        bl, dage, delts, hcomp, hd, hjm, hp, hr, hs, hs1,
                        hs2, hstot, ht, jenv, lnew, m, omega, smass, teffl,
                        senvold, lnewfit)
                    # The composition at the start of the step.
                    nspec = 11
                    if nml.lexcom:
                        nspec = 15
                    for i in range(1, m + 1):
                        for j in range(1, nspec + 1):
                            old.hcompp[j, i] = hcomp[j, i]
                    # Mixing, burning and settling.
                    itlvl = 1
                    ddage, jcore, jenv = mix(
                        delts, hcomp, hd, hl, hp, hr, hs, hs1, hs2, hstot, ht,
                        itlvl, lc, m, ddage, jcore, jenv, mxzone, mxzon0,
                        teffl, cmixl)
                    ddage = delts / c.csecyr
                    dage = dage + 1.0e-9 * ddage
                # Rezone, except for He-flash runs, which only keep the old
                # model for PTIME.
                if not nml.lkuthe:
                    m, lreset, jcore, jenv, sjtot, skerot = hpoint(
                        m, hstot, hs, hs1, hs2, ht, hp, hr, lc, hl, hd,
                        hcomp, ps, rs, estore, istore, lreset, jxbeg, model,
                        lshell, jcore, jenv, omega, eta2, r0, hi, hjm,
                        hkerot, sjtot, skerot, fp, ft, hg, qiw, teffl, cmixl)
                    deuter_mod.jcz = jenv
                else:
                    for i in range(1, m + 1):
                        old.hpo[i] = hp[i]
                        old.hto[i] = ht[i]
                        old.hro[i] = hr[i]
                        old.hlo[i] = hl[i]
                        old.hdo[i] = hd[i]
                    old.tefflo = teffl
                    old.mo = m
                if nml.lrot:
                    for i in range(1, m + 1):
                        hkeroto[i] = hkerot[i]
                # Light elements at the start of the step.
                if nml.lexcom:
                    prevcz_mod.rbase0 = 0.0
                    jenv0 = jenv
                    if nml.lovste:
                        liov_mod.pscal0 = nml.alphae * math.exp(
                            c.clndp * (hp[jenv] + 2.0 * hr[jenv] - hd[jenv]
                                       - c.cgl - hs[jenv]))
                    else:
                        liov_mod.pscal0 = 0.0
                    lirate88(hcomp, hd, ht, m, 1)
                # Flags for CRRECT: LNEW0, a new envelope triangle on the
                # first iteration; LFINI, converged; LARGE, diverged.
                if nml.lnew0:
                    lnew = True
                if not lage:
                    delts = -abs(delts)
                ctol_mod.fcorr = abs(nml.fcorr0) - nml.fcorri
                itdone = 0
                large = False
                lfini = False
                # First guess of the structure.
                if not nml.lnews or delts <= 0.0:
                    for i in range(1, m + 1):
                        htt[i] = 0.0
                        hpp[i] = 0.0
                        ent.htt0[i] = 0.0
                        ent.hpp0[i] = 0.0
                        ent.hll0[i] = 0.0
                        ent.hrr0[i] = 0.0
                        hhc[i] = 0.0
                else:
                    # Extrapolate from the previous model's rates.
                    for i in range(1, m + 1):
                        dt = ent.htt0[i] * delts
                        dp = ent.hpp0[i] * delts
                        dl = hl[i] * ent.hll0[i] * delts
                        dr = ent.hrr0[i] * delts
                        htt[i] = dt
                        hpp[i] = dp
                        ht[i] = ht[i] + dt
                        hp[i] = hp[i] + dp
                        hl[i] = hl[i] + dl
                        hr[i] = hr[i] + dr
                        hhc[i] = 0.0

                # Level 1: the previous envelope triangle.
                niter = nml.niter1
                lsbc = False
                itlvl = 1
                (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
                 latmo, lderiv, lmix, locond, qdt, qdp, itdone) = _crrect(
                    niter, lsbc, itlvl)
                if large:
                    continue                              # go to 15
                # Level 2: check the surface boundary condition.
                lsbc = True
                niter = nml.niter2
                itlvl = 2
                (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
                 latmo, lderiv, lmix, locond, qdt, qdp, itdone) = _crrect(
                    niter, lsbc, itlvl)
                if large:
                    continue                              # go to 15
                # Store the rates of change for the next step's guess.
                if delts > 0.0:
                    for ii in range(1, m + 1):
                        ent.htt0[ii] = htt[ii] / delts
                        ent.hpp0[ii] = hpp[ii] / delts
                        ent.hll0[ii] = (2.0 * (hl[ii] - old.hlo[ii])
                                        / (hl[ii] + old.hlo[ii]) / delts)
                        ent.hrr0[ii] = (hr[ii] - old.hro[ii]) / delts
                # Level 3.
                lsbc = False
                niter = nml.niter3
                itlvl = 3
                (lfini, large, lnew, lreset, lsbc, ftri, istore, bl, teffl,
                 latmo, lderiv, lmix, locond, qdt, qdp, itdone) = _crrect(
                    niter, lsbc, itlvl)
                if large:
                    continue                              # go to 15
                if not nml.lrot:
                    nml.itdif1 = 1
                # With rotation, iterate between structure and rotation
                # ITDIF1 times from the start-of-step j and composition.
                if nml.itdif1 > 1:
                    for i in range(1, m + 1):
                        hjm0[i] = hjm[i]
                        for j in range(1, 16):
                            hcompp0[j, i] = old.hcompp[j, i]
                goto15 = False
                for itrot in range(1, nml.itdif1 + 1):
                    if itrot > 1:
                        for i in range(1, m + 1):
                            for j in range(1, 16):
                                old.hcompp[j, i] = hcompp0[j, i]
                    # Level 4: repeats the iteration between mixing and
                    # structure (not for models where semiconvection matters,
                    # where it oscillates; for high-precision work such as
                    # solar models) (7/91).
                    niter = nml.niter4
                    lsbc = False
                    itlvl = 4
                    (lfini, large, lnew, lreset, lsbc, ftri, istore, bl,
                     teffl, latmo, lderiv, lmix, locond, qdt, qdp,
                     itdone) = _crrect(niter, lsbc, itlvl)
                    if not lfini:
                        # Not converged: treat as diverged.
                        large = True
                        goto15 = True
                        break
                    # Converged: fully mix the convection zones (LIMIX, so that
                    # the correct gradients are used) (G. Somers 6/14).
                    limix = True
                    mixcz(hcomp, hs2, lc, m)
                    # The total age, for disk locking (MHP 9/94).
                    disk_mod.sage = dage
                    if nml.lrot:
                        if itrot > 1:
                            for i in range(1, m + 1):
                                hjm[i] = hjm0[i]
                        # Print the rotation summary at the final age.
                        if (rl_.lendag[nk] and nml.endage[nk] > 0.0
                                and abs(nml.endage[nk] - dage * 1.0e9
                                        - ddage) <= 1.0):
                            lprt0 = True
                        else:
                            lprt0 = False
                        # The new Omega: transport and wind loss.
                        ljdot = nml.ljdot0
                        bl, dwmax, ljdot, jenv0 = getw(
                            bl, delts, dwmax, fp, ft, hcomp, hd, hjm, hl, hp,
                            hr, hs, hs1, hs2, hstot, ht, lc, ljdot, m, smass,
                            teffl, eta2, hg, hi, omega, qiw, r0, jenv0,
                            cmixl, lprt0)
                        # The distortion terms.
                        fpft(hd, hr, hs, m, omega, eta2, fp, ft, hg, r0)
                if goto15:
                    continue                              # go to 15
                break

            # The H-burning shell and convection-zone edges again.
            jcore, jenv, jxbeg, jxend, jxmid, lshell = findsh(
                hcomp, hl, lc, m, jcore, jenv, jxbeg, jxend, jxmid, lshell)
            # Light-element burning (with rotational mixing, GETW does it).
            if nml.lexcom and model >= 0 and delts > 0.0:
                if not nml.lrot or not nml.linstb:
                    jcore, jenv, nrzone, nzone, nzone0 = convec(
                        hcomp, hd, hp, hr, hs, ht, lc, m, mrzone, mxzone,
                        mxzon0, jcore, jenv, nrzone, nzone, nzone0)
                    jenv1 = jenv
                    if nml.lovste:
                        liov_mod.pscal1 = nml.alphae * math.exp(
                            c.clndp * (hp[jenv] + 2.0 * hr[jenv] - hd[jenv]
                                       - c.cgl - hs[jenv]))
                    else:
                        liov_mod.pscal1 = 0.0
                    lirate88(hcomp, hd, ht, m, 2)
                    liburn(delts, hcomp, hr, hs1, hs2, ht, jenv1, jenv0, m)
            # Restore the earlier fitting point after accretion (MHP 07/02).
            if lnewfit:
                senvold, sjtot, skerot, m, lnew = getnewenv(
                    senvold, hcomp, hd, hl, hp, hr, hs, hs1, hs2, hstot, ht,
                    lc, eta2, hi, hjm, qiw, r0, hkerot, bl, sjtot, skerot,
                    teffl, m, lnew, cmixl)
                fpft(hd, hr, hs, m, omega, eta2, fp, ft, hg, r0)
                lnewfit = False
            # The time step for the next model (HTIMER also locates the
            # H-burning shell).
            delts = abs(delts)
            delts0 = delts
            delts, deltsh, ddage, dwmax = htimer(
                delts, deltsh, m, hd, hl, hs1, hs2, ht, hcomp, jcore, jxmid,
                tlumx, dage, ddage, nk, hp, hr, omega, dwmax, jxbeg, teffl)
            # At a fixed final age, write the step just taken rather than a
            # zero one, so that the sequence can be continued.
            if rl_.lendag[nk] and nml.endage[nk] > 0.0:
                if nml.endage[nk] - dage * 1.0e9 <= 1.0:
                    # Question (see the page): the floor never acts.
                    delts = max(delts0, 1.0e-3 * dage * c.csecyr)
                    ddage = delts / c.csecyr
                else:
                    delts0 = delts
            else:
                delts0 = delts
            if rl_.iresca[nk] != 2:
                model = model + 1
            # The output driver.
            gl, lpunch = wrtout(
                hcomp, hd, hl, hp, hr, hs, hs1, ht, lc, m, model, dage,
                ddage, smass, teffl, bl, gl, lshell, jxbeg, jxmid, jxend,
                jcore, jenv, tlumx, trit, tril, ps, ts, rs, cfenv, ftri,
                hstot, omega, lpunch, fp, ft, eta2, r0, hjm, hi, sjtot,
                skerot, hs2, cmixl, xmsol, compmix)
            # Stop conditions (MHP 10/24).
            if (rl_.lendag[nk] and nml.endage[nk] > 0.0
                    and (nml.endage[nk] - dage * 1.0e9) <= 1.0):
                break                                     # go to 110
            lend_kind = False
            if (rl_.lendag[nk] and nml.end_dcen[nk] > 0.0
                    and hcomp[12, 1] < nml.end_dcen[nk]):
                print(f"CENTRAL D {hcomp[12, 1]:12.4E} BELOW STOP VALUE "
                      f"{nml.end_dcen[nk]:12.4E}")        # format 104
                lend_kind = True
            if (rl_.lendag[nk] and nml.end_xcen[nk] > 0.0
                    and hcomp[1, 1] < nml.end_xcen[nk]):
                print(f"CENTRAL X {hcomp[1, 1]:12.4E} BELOW STOP VALUE "
                      f"{nml.end_xcen[nk]:12.4E}")        # format 105
                lend_kind = True
            if (rl_.lendag[nk] and nml.end_ycen[nk] > 0.0
                    and hcomp[2, 1] < nml.end_ycen[nk]):
                print(f"CENTRAL Y {hcomp[2, 1]:12.4E} BELOW STOP VALUE "
                      f"{nml.end_ycen[nk]:12.4E}")        # format 106
                lend_kind = True
            if lend_kind:
                nml.lpulse = lsavpu
                sound_mod.lsound = True
                lprt0 = True
                break                                     # go to 110
            # Is the model near the desired Teff and L? If not, rescale and
            # try again.
            if nml.lcalst and not lstar:
                if nk % 2 == 0:
                    if modeln == 1:
                        teffi = 10.0**teffl  # unused
                    else:
                        alri, lstar, lpassr = chkscal(bl, teffl, dage, nk,
                                                      alri, lstar, lpassr)
                        if lpassr:
                            goto200 = True
                            break                         # go to 200
        # label 100 (end of the model loop)
        if goto200:
            continue                                      # go to 200

        # Label 110: store the last model in ISTOR if asked.
        if nml.lstore and nml.lstpch and lpunch:
            putstore(hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps, ts,
                     rs, cfenv, ftri, tlumx, jcore, jenv, model, m, smass,
                     teffl, bl, hstot, dage, ddage, omega, hs1, eta2, r0, fp,
                     ft, hjm, hi, cmixl, xmsol, compmix)
            lpunch = False
        # The solar calibration, in three kind cards (MHP 1/93, 5/96).
        if nml.lcals:
            acdpth_mod.ladon = False   # calcad off, for speed (JVS)
            if nk % 3 == 0:
                rll = (0.5 * (bl + c.clsunl - c.c4pil - c.csigl
                              - 4.0 * teffl) - c.crsunl)
                # The solar Z/X as an observable (MHP 06/13).
                solzx = hcomp[3, m] / hcomp[1, m]
                dldx, drdx, dlda, drda, lsol = chkcal(
                    bl, rll, nk, solzx, dldx, drdx, dlda, drda, lsol)
                nml.lptime = lptimesv
                nml.kttau = kttausv
                if lsol:
                    break                                 # go to 250
                else:
                    iconv = iconv + 1
                    # Stop after 15 attempts (the Fortran comment still says 10) (MHP 6/97).
                    if iconv >= 15:
                        break                             # go to 250
                    if nml.lpulse:
                        # The pulsation output is not final yet: rewind.
                        io.rewind(io.iopmod)
                        io.rewind(io.iopenv)
                        io.rewind(io.iopatm)
        # The star of given Teff and L is matched.
        if nml.lcalst and lstar and nk % 2 == 0:
            break                                         # go to 250
        # label 200
    # Label 250: the solar calibration's summary.
    if nml.lcals and nml.lsnu and lsol:
        for unit in (io.ilast, io.ifirst, io.idebug, io.itrack, io.ishort,
                     io.imodpt, io.istor):
            io.rewind(unit)
        ztox = hcomp[3, m] / hcomp[1, m]
        # Runs needed, initial X and alpha, final dL/dX, dR/dX, dL/dalpha,
        # dR/dalpha.
        io.write(io.isnu, f" {iconv:2d}{xguess:10.6f}{aguess:10.6f}"
                 f"{dldx:11.4E}{drdx:11.4E}{dlda:11.4E}"
                 f"{drda:11.4E}")                          # format 1518
        # Neutrino fluxes (see ENGEB).
        io.write(io.isnu, f" {fluxes_mod.clsnu:8.3f}{fluxes_mod.gasnu:8.3f}"
                 + "".join(f"{fluxes_mod.fluxtot[i]:10.3E}"
                           for i in range(1, 11)))         # format 1516
        # Structure: Tc, rho_c, Pc, Xc, Zc.
        tcen = 10.0**(cent_mod.tcenter - 6.0)
        pcen = 10.0**(cent_mod.pcenter - 17.0)
        dcen = 10.0**cent_mod.dcenter
        io.write(io.isnu, f" {tcen:7.3f}{dcen:7.2f}{pcen:6.3f}"
                 f"{hcomp[1, 1]:8.5f}{hcomp[3, 1]:8.5f}")  # format 1517
        # Initial alpha, Y, Z; final L, R.
        yini = 1.0 - rl_.rescal[2, nk - 2] - rl_.rescal[3, nk - 2]
        zini = rl_.rescal[3, nk - 2]
        io.write(io.isnu, f"{nml.cmixla[nk]:7.4f}{yini:8.5f}{zini:8.5f}"
                 f"{bl:10.3E}{rll:10.3E}")                 # format 1521
        # The convection zone's base (R, M), surface Y, Z and Z/X.
        io.write(io.isnu, f"{cent_mod.czsurf_r:8.5f}{cent_mod.czsurf_m:9.6f}"
                 f"{hcomp[2, m]:8.5f}{hcomp[3, m]:8.5f}"
                 f"{ztox:9.6f}")                           # format 1522
        # Energy generation fractions: pp I, II, III, CNO, gravitational.
        io.write(io.isnu, "".join(f"{tlumx[j]:10.3E}"
                                  for j in (1, 2, 3, 4, 7)))  # format 1523
    raise SystemExit  # STOP


if __name__ == "__main__":
    main()
