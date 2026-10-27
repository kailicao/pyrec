"""Python version of chkcal.f (YREC 09c3b9b): the solar calibration's
check and Newton step."""

import math

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod

# CHKCAL has SAVE: BLP and RLP are stored on one call and differenced on
# the next; the others keep their values too.
blp = 0.0
rlp = 0.0
dx = 0.0
da = 0.0
dz = 0.0
zl = 0.0
zlp = 0.0


def chkcal(bl, rl, n, solzx, dldx, drdx, dlda, drda, lsol):
    """Check a solar model against L_sun, R_sun (and the solar Z/X); if it
    is not calibrated, set new X, alpha (and Z) for kind cards n+1 to n+3
    by a Newton step with fixed derivatives.

    Translates: chkcal.f#CHKCAL (09c3b9b)

    Arguments: bl, log L/L_sun; rl, log R/R_sun; n, the current kind card;
    solzx, the model's surface Z/X; the derivatives dldx, drdx, dlda, drda
    (set here); lsol, whether the calibration has converged.

    Returns (dldx, drdx, dlda, drda, lsol). The new X, alpha and Z go to
    runloop_mod.rescal and nmlpar_mod.xenv0a, cmixla, zenv0a.
    """
    global blp, rlp, dx, da, dz, zl, zlp
    nml = nmlpar_mod
    rescal = runloop_mod.rescal
    iowr = iounits_mod.iowr

    # Converged? Question (see the page): Z/X not converged.
    if abs(bl) < nml.toll and abs(rl) < nml.tolr:
        if nml.lcalsolzx and nml.calsolzx > 0.0:
            zl = math.log10(solzx) - math.log10(nml.calsolzx)
            if abs(zl) < nml.tolz:
                lsol = True
                return dldx, drdx, dlda, drda, lsol
        else:
            lsol = True
            return dldx, drdx, dlda, drda, lsol

    # Fixed, empirical derivatives (LLP 6/28/09).
    dldx = -3.78   # empirical result:  -3.783    RMS error .070
    drdx = -0.89   # empirical result:  -0.890    RMS error .048
    # Question (see the page): a derivative off by ten?
    dlda = 0.0139  # empirical result:  +0.139    RMS error .0022
    drda = -0.050  # empirical result:  -0.0504   RMS error .0059

    # GOTO 1234 bypasses the older method, which measured the derivatives
    # with extra runs at kind cards 3, 6 and 9.
    # Question (see the page): dead code.
    if False:
        if n == 3:
            rescal[2, n + 1] = rescal[2, 1] + dx
            nml.xenv0a[n + 1] = rescal[2, n + 1]
            nml.xenv0a[n + 2] = nml.xenv0a[n + 1]
            nml.xenv0a[n + 3] = nml.xenv0a[n + 1]
            blp = bl
            rlp = rl
            return dldx, drdx, dlda, drda, lsol
        elif n == 6:
            nml.cmixla[n + 1] = nml.cmixla[1] + da
            nml.cmixla[n + 2] = nml.cmixla[n + 1]
            nml.cmixla[n + 3] = nml.cmixla[n + 1]
            rescal[2, n + 1] = rescal[2, 1]
            nml.xenv0a[n + 1] = nml.xenv0a[1]
            nml.xenv0a[n + 2] = nml.xenv0a[n + 1]
            nml.xenv0a[n + 3] = nml.xenv0a[n + 1]
            blp = bl
            rlp = rl
            return dldx, drdx, dlda, drda, lsol
        elif n == 9:
            dlda = (bl - blp) / da
            drda = (rl - rlp) / da
            iounits_mod.write(iowr, f" DA,DLDA,BL,BLP,DRDA,RL,RLP:  {da} "
                              f"{dlda} {bl} {blp} {drda} {rl} {rlp}")
            da = ((bl * drdx / dldx - rl) / (drda - dlda * drdx / dldx))
            dx = -(bl + dlda * da) / dldx
            nml.cmixla[n + 1] = nml.cmixla[1] + da
            nml.cmixla[n + 2] = nml.cmixla[n + 1]
            nml.cmixla[n + 3] = nml.cmixla[n + 1]
            rescal[2, n + 1] = rescal[2, 1] + dx
            nml.xenv0a[n + 1] = rescal[2, n + 1]
            nml.xenv0a[n + 2] = rescal[2, n + 1]
            nml.xenv0a[n + 3] = rescal[2, n + 1]
            blp = bl
            rlp = rl
            iounits_mod.write(iowr, f" New X, Old X, Calc DX:  "
                              f"{rescal[2, n + 1]} {rescal[2, 1]} {dx}")
            iounits_mod.write(iowr, f" New A, Old A, Calc DA:  "
                              f"{nml.cmixla[n + 1]} {nml.cmixla[1]} {da}")
            return dldx, drdx, dlda, drda, lsol

    # label 1234: the Newton step, from the values of card n-2.
    da = ((bl * drdx / dldx - rl) / (drda - dlda * drdx / dldx))
    dx = -(bl + dlda * da) / dldx
    nml.cmixla[n + 1] = nml.cmixla[n - 2] + da
    nml.cmixla[n + 2] = nml.cmixla[n + 1]
    nml.cmixla[n + 3] = nml.cmixla[n + 1]
    if nml.lcalsolzx:
        rescal[3, n + 1] = rescal[3, n - 2] * nml.calsolzx / solzx
        dz = rescal[3, n + 1] - rescal[3, n - 2]
        nml.zenv0a[n + 1] = rescal[3, n + 1]
        nml.zenv0a[n + 2] = rescal[3, n + 1]
        nml.zenv0a[n + 3] = rescal[3, n + 1]
    rescal[2, n + 1] = rescal[2, n - 2] + dx
    nml.xenv0a[n + 1] = rescal[2, n + 1]
    nml.xenv0a[n + 2] = rescal[2, n + 1]
    nml.xenv0a[n + 3] = rescal[2, n + 1]
    iounits_mod.write(iowr, f" New BL, Old BL, Delta BL:  {bl} {blp} "
                      f"{bl - blp}")
    iounits_mod.write(iowr, f" New RL, Old RL, Delta RL:  {rl} {rlp} "
                      f"{rl - rlp}")
    iounits_mod.write(iowr, f" New X, Old X, DX:  {rescal[2, n + 1]} "
                      f"{rescal[2, n - 2]} {dx}")
    iounits_mod.write(iowr, f" New A, Old A, DA:  {nml.cmixla[n + 1]} "
                      f"{nml.cmixla[n - 2]} {da}")
    if nml.lcalsolzx:
        iounits_mod.write(iowr, f" New Z, Old Z, DZ:  {rescal[3, n + 1]} "
                          f"{rescal[3, n - 2]} {dz}")
    blp = bl
    rlp = rl
    zlp = zl
    return dldx, drdx, dlda, drda, lsol
