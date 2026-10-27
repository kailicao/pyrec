"""Python version of checkj.f (YREC 09c3b9b): checks one iteration of
angular-momentum transport."""

import numpy as np

from pyrec import difad2_mod
from pyrec import difad_mod
from pyrec import egrid_mod
from pyrec import errmom_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import oldrot_mod
from pyrec import quadd_mod
from pyrec import temp2_mod
from pyrec.getrot import getrot
from pyrec.solid import solid

# Locals kept between calls (SAVE): the largest change of j over each
# iteration, and where it occurred.
# Question (see the page): indexed by IT, so ITDIF2 > 16 overflows them.
delji = np.zeros(17)
imaxi = np.zeros(17, dtype=int)


def checkj(hd, hjmp, hjmsav, hr, hs, hs2, hv, imin, imax, it, lcz, m, dt,
           eta2, hi, hjm, iredo, lok, lredo, omega, qiw, r0, wsav, id_,
           idm, lokad):
    """Check one iteration of angular-momentum diffusion: stop on j <= 0,
    test convergence against nmlpar_mod.djok, get Omega (and the moment
    of inertia) from j with getrot, and remove reversals of the Omega
    gradient by local solid-body rotation (solid).

    Translates: checkj.f#CHECKJ (09c3b9b)

    Arguments: the model hd, hr, hs, hs2 (logs, and the shell masses), hv
    (the diffusion velocities), lcz (convective for angular momentum);
    hjmp, j after the previous iteration (set to hjm here); hjmsav, j at
    the start of the step; imin, imax, the range of hv; it, the iteration;
    m, the number of points; dt, the diffusion step; the arrays eta2, hi,
    omega, qiw, r0 (set by getrot and solid); hjm, j after this iteration
    (changed by solid); iredo; wsav, Omega at the start of the step (for
    the print); id_, the printed shells (ID); idm; lokad, convergence
    assumed on entry.

    Returns (dt, iredo, lok, lredo, idm): the step (halved in a branch
    that never runs), the count of cuts, converged or not, a cut needed,
    and the number of printed shells (set only in the switched-off print).
    """
    n = nmlpar_mod

    # Check for negative j.
    lok = False
    lredo = False
    if lokad:
        lok = True
        lredo = False
    for i in range(1, m + 1):
        if hjm[i] <= 0.0:
            iredo = iredo + 1
            # Question (see the page): IREDO.GT.0, not .GT.3, so the
            # time-step cut below never runs.
            if iredo > 0:
                msg = (" " + ">" * 39 + "<" * 39 + "\n"
                       + "     ERROR IN SR CHECKJ\n"
                       + f"     NEGATIVE J/M ENCOUNTERED IN ZONE{i:5d}"
                       + " AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED\n"
                       + "RUN STOPPED")
                print(msg)
                iounits_mod.write(iounits_mod.ishort, msg)
                raise SystemExit(msg)
            else:
                lredo = True
                dt = 0.5 * dt
                msg = ("     ERROR IN SR CHECKJ\n"
                       + f"     TIMESTEP CUT, NUMBER{iredo:5d}"
                       + f" DUE TO NEGATIVE J/M IN ZONE{i:5d}")
                print(msg)
                iounits_mod.write(iounits_mod.ishort, msg)
                return dt, iredo, lok, lredo, idm  # GOTO 240

    # Convergence: the largest fractional change of j over the iteration,
    # and where it occurs; hjmp becomes hjm.
    delji[it] = (hjm[1] - hjmp[1]) / hjmp[1]
    hjmp[1] = hjm[1]
    imaxi[it] = 1
    for i in range(2, m + 1):
        dj = (hjm[i] - hjmp[i]) / hjmp[i]
        if abs(dj) > abs(delji[it]):
            delji[it] = dj
            imaxi[it] = i
        hjmp[i] = hjm[i]
    if abs(delji[it]) <= n.djok:
        lok = True

    # Omega from the new j. The error in the diffusion is of order DJOK,
    # so the tolerances of the moment-of-inertia iteration are relaxed
    # except at the last iteration, and restored after.
    toler0 = errmom_mod.toleri
    acc = n.acfpft
    if it < n.itdif2 and not lok:
        errmom_mod.toleri = max(n.djok * 1.0e-2, toler0)
        n.acfpft = max(n.djok * 1.0e-2, acc)
    getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0)
    errmom_mod.toleri = toler0
    n.acfpft = acc

    # Reversals of the Omega gradient, from the surface inward: enforce
    # solid-body rotation in the offending pair of shells.
    wold = oldrot_mod.wold
    i = m
    ibot = m
    while True:  # label 20
        # Question (see the page): a threshold of 1.0, not 1.0d-15.
        if omega[i] - omega[i - 1] > 1.0:
            if wold[i] - wold[i - 1] > 1.0e-15:
                # The previous gradient was positive too: leave it.
                i = ibot - 1
                ibot = i
                if i > 1:
                    continue  # GOTO 20
                # else GOTO 130
            else:
                # ITOP, the uppermost unstable shell, and IBOT, the
                # lowest, each widened over an adjacent convection zone.
                itop = i
                if lcz[itop] and itop < m:
                    j = itop + 1
                    while j <= m and lcz[j]:
                        j += 1
                    itop = j - 1  # label 40
                ibot = i - 1
                if lcz[ibot] and ibot > 1:
                    j = ibot - 1
                    while j >= 1 and lcz[j]:
                        j -= 1
                    ibot = j + 1  # label 60
                solid(hd, hjm, hr, hs, hs2, ibot, itop, eta2, hi, omega,
                      qiw, r0, m)
                # Widen the range while its boundaries show new reversals.
                while True:  # label 70
                    lredo = False
                    if ibot > 1:
                        if omega[ibot] - omega[ibot - 1] > 1.0e-15:
                            if wold[ibot] - wold[ibot - 1] < 1.0e-15:
                                lredo = True
                                ibot = ibot - 1
                                if lcz[ibot] and ibot > 1:
                                    j = ibot - 1
                                    while j >= 1 and lcz[j]:
                                        j -= 1
                                    ibot = j + 1  # label 90
                    if itop < m:
                        if omega[itop + 1] - omega[itop] > 1.0e-15:
                            if wold[itop + 1] - wold[itop] < 1.0e-15:
                                lredo = True
                                itop = itop + 1
                                if lcz[itop] and itop < m:
                                    j = itop + 1
                                    while j <= m and lcz[j]:
                                        j += 1
                                    itop = j - 1  # label 110
                    if lredo:
                        solid(hd, hjm, hr, hs, hs2, ibot, itop, eta2, hi,
                              omega, qiw, r0, m)
                        continue  # GOTO 70
                    break
                if it == n.itdif2 or lok:
                    print(f"     OMEGA GRADIENT REVERSAL BETWEEN ZONES "
                          f"{ibot:5d} AND {itop:5d} ITERATION {it:5d}")
        # label 130: the next shell.
        i = ibot - 1
        ibot = i
        if i <= 1:
            break

    # Output at the end of the diffusion step.
    if it == n.itdif2 or lok:
        # The largest fractional change of j over the time step.
        delj = 0.0
        ijmax = 0
        ii = 1
        for i in range(ii, m + 1):
            dj = (hjm[i] - hjmsav[i]) / hjmsav[i]
            if abs(dj) > abs(delj):
                delj = dj
                ijmax = i
        if not lokad:
            line = f" MAX D(J/M)/(J/M){delj:12.3E} AT PT.{ijmax:5d} BY ITERATION"
            for j in range(1, it + 1):
                if (j - 1) % 5 == 0:
                    line += "\n"
                line += f" {delji[j]:11.3E}{imaxi[j]:4d}"
            print(line)
        # The detailed output to IMODPT is switched off (IF(.TRUE.)GOTO
        # 240, G. Somers 11/14); kept here as in the Fortran.
        # Question (see the page): can it be removed?
        if False:
            idm = _detailed_print(hjm, hjmsav, hv, imin, imax, m, omega,
                                  wsav, id_, idm)
    # label 240
    return dt, iredo, lok, lredo, idm


def _detailed_print(hjm, hjmsav, hv, imin, imax, m, omega, wsav, id_,
                    idm):
    """The switched-off detailed print of CHECKJ (after IF(.TRUE.)GOTO
    240): Omega, j and the diffusion velocities at every NPRTPT-th shell.
    Returns idm."""
    n = nmlpar_mod
    t2 = temp2_mod
    imodpt = iounits_mod.imodpt
    nprtpt = n.nprtpt
    if nprtpt > m:
        return idm
    iounits_mod.write(imodpt, " SHELL   OMEGA     DEL OMEGA      J/M"
                              "       DEL J/M")
    # The first point is always printed.
    id_[1] = 1
    idm = 2
    ibeg = max(imin, nprtpt)
    iend = min(imax, int(imax / nprtpt) * nprtpt)
    for j in range(ibeg, iend + 1, nprtpt):
        id_[idm] = j
        idm = idm + 1
    # The outermost point (or the base of a surface CZ) is always printed.
    if id_[idm - 1] != imax:
        id_[idm] = imax
    else:
        idm = idm - 1
    for i in range(1, idm + 1):
        k = id_[i]
        iounits_mod.write(imodpt, f" {k:5d}{omega[k]:12.3E}"
                          f"{omega[k] - wsav[k]:12.3E}{hjm[k]:12.3E}"
                          f"{hjm[k] - hjmsav[k]:12.3E}")
    iounits_mod.write(imodpt, " SHELL    VES0         VES       VGSF0"
                              "        VGSF         VSS         RAT"
                              "        VTOT       LENGTH        VMU")
    for i in range(1, idm + 1):
        k = id_[i]
        vals = (t2.ves0[k], t2.ves[k], t2.vgsf0[k], t2.vgsf[k], t2.vss[k],
                quadd_mod.rat[k], hv[k], t2.hle[k], t2.vmu[k])
        iounits_mod.write(imodpt, f" {k:5d}"
                          + "".join(f"{v:12.3E}" for v in vals))
    if n.ldifad:
        eg = egrid_mod
        ntot = eg.ntot

        def row7(i):
            vals = (eg.chi[i], t2.ves[i], difad2_mod.vesa[i],
                    difad2_mod.vesd[i], eg.echi[i], difad_mod.ecod3[i],
                    difad_mod.ecod4[i])
            return f" {i:5d}" + "".join(f"{v:12.3E}" for v in vals)

        if idm == ntot:
            for i in range(1, idm + 1):
                iounits_mod.write(imodpt, row7(i))
        elif idm < ntot:
            for i in range(1, idm + 1):
                iounits_mod.write(imodpt, row7(i))
            for i in range(idm + 1, ntot + 1):
                iounits_mod.write(imodpt, f" {i:5d}" + " " * 48
                                  + f"{eg.echi[i]:12.3E}"
                                  f"{difad_mod.ecod3[i]:12.3E}"
                                  f"{difad_mod.ecod4[i]:12.3E}")
        else:
            for i in range(1, ntot + 1):
                iounits_mod.write(imodpt, row7(i))
            for i in range(ntot + 1, idm + 1):
                vals = (eg.chi[i], t2.ves[i], difad2_mod.vesa[i],
                        difad2_mod.vesd[i])
                iounits_mod.write(imodpt, f" {i:5d}"
                                  + "".join(f"{v:12.3E}" for v in vals))
    return idm
