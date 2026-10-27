"""Python version of checkc.f (YREC 09c3b9b): the composition checked after
rotational mixing, and the mean molecular weight recomputed."""

from pyrec import comp2_mod
from pyrec import comp_mod
from pyrec import iounits_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# DATA ATOMWT: H, 4He, Z (as carbon), 3He; index 0 unused.
ATOMWT = (None, 1.007825, 4.002603, 12.0, 3.01603)


def checkc(hcomp, it, lprt, m, dt, iredo, lok, lredo):
    """Check abundances after a step of rotational mixing: repair tiny
    anomalies, halve the time step otherwise (stop after 3 cuts), and
    recompute the mean molecular weight AMUM of a fully ionized gas.

    Translates: checkc.f#CHECKC (09c3b9b)

    Arguments: hcomp[species, shell] (changed in place where repaired); it,
    the iteration on the diffusion coefficients; lprt, print the largest
    change at the last iteration; m, the number of shells; dt, the time
    step; iredo, the number of cuts so far; lok, lredo.

    Returns (dt, iredo, lok, lredo): dt halved, iredo increased and
    lok = False, lredo = True after an anomaly that is not repaired.
    """
    hcompp = oldmod_mod.hcompp
    # Before the last iteration only X, Y and 3He.
    jend = 11 if it == nmlpar_mod.itdif2 else 4
    lredo = False
    for j in range(1, jend + 1):  # label 20
        if j == 3:  # Z is not diffused as a unit.
            continue
        for i in range(2, m):  # label 10
            if hcomp[j, i] < 0.0 or hcomp[j, i] > 1.0:
                # Question (see the page): the threshold is relative to the
                # surface value, which may be zero.
                if abs(hcomp[j, i]) < 1.0e-5 * hcomp[j, m]:
                    # Question (see the page): never true in this loop.
                    if i == 1 or i == m:
                        hcomp[j, i] = 0.0
                    else:
                        hcomp[j, i] = max(0.0, 0.5 * (hcomp[j, i + 1]
                                                      + hcomp[j, i - 1]))
                    continue
                iredo = iredo + 1
                if iredo > 3:
                    msg = (" " + ">" * 78 + "\n ERROR IN SR CHECKC\n"
                           f" ANOMALOUS COMP NUMBER{j:2d} IN ZONE{i:5d}"
                           f" ABUNDANCE {hcomp[j, i]:12.3E}\n"
                           " 3 ATTEMPTS AT TIMESTEP CUTTING FAILED\n"
                           "RUN STOPPED")  # format 1010
                    print(msg)
                    iounits_mod.write(iounits_mod.ishort, msg)
                    raise SystemExit
                lredo = True
                lok = False
                dt = 0.5 * dt
                msg = (" ERROR IN SR CHECKC\n"
                       f" TIMESTEP CUT NUMBER {iredo:2d}"
                       f" DUE TO ANOMALOUS COMP NUMBER{j:2d} IN ZONE{i:5d}"
                       f" ABUNDANCE {hcomp[j, i]:12.3E}")  # format 1015
                print(msg)
                iounits_mod.write(iounits_mod.ishort, msg)
                # GOTO 100: the end, without recomputing mu.
                return dt, iredo, lok, lredo

    # The largest fractional change over the step.
    if it == nmlpar_mod.itdif2 and lprt:
        delc = 0.0
        icmax = 0
        jcmax = 0
        for j in range(1, jend + 1):  # label 40
            if j == 3:
                continue
            cminchk = max(1.0e-6 * hcomp[j, m], 1.0e-20)
            for i in range(1, m + 1):  # label 30
                if hcompp[j, i] < cminchk:
                    continue
                dcomp = (hcomp[j, i] - hcompp[j, i]) / hcompp[j, i]
                if abs(dcomp) > abs(delc):
                    delc = dcomp
                    icmax = i
                    jcmax = j
        print(f" MAX FRAC.COMP.CHANGE{delc:12.3E} SPECIES{jcmax:2d}"
              f" AT PT.{icmax:5d}")  # format 50
        if nmlpar_mod.lexcom:
            print(f"     NEW SURFACE LI{hcomp[14, m]:14.4E}"
                  f"OLD VALUE{hcompp[14, m]:14.4E}")  # format 60

    # Mean molecular weight of a fully ionized gas.
    if it > 1:
        for i in range(1, m + 1):  # label 90
            dfx1 = hcomp[1, i] - comp_mod.xenv
            dfx2 = hcomp[2, i] - comp2_mod.yenv
            dfx3 = hcomp[3, i] - comp_mod.zenv
            dfx4 = hcomp[4, i] - comp2_mod.y3env
            temp = (comp_mod.amuenv + dfx1 / ATOMWT[1] + dfx2 / ATOMWT[2]
                    + dfx3 / ATOMWT[3] + dfx4 / ATOMWT[4])
            amu2 = 1.0 / temp
            temp = (hcomp[1, i] / ATOMWT[1]
                    + 2.0 * (hcomp[4, i] / ATOMWT[4]
                             + hcomp[2, i] / ATOMWT[2])
                    + 0.5 * hcomp[3, i])
            emu2 = 1.0 / temp
            mdphy_mod.amum[i] = amu2 * emu2 / (amu2 + emu2)
    return dt, iredo, lok, lredo
