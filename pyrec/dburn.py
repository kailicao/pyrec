"""Python version of dburn.f (YREC 09c3b9b): deuterium burning over a
step in one shell or one mixed zone."""

import math

from pyrec import deuter_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# SAVE: only SUMM is read before it is set (the accretion branch for a
# single-shell surface zone), so it is the one kept between calls.
summ = 0.0


def dburn(ibegin, iend, m, hs2, hcomp, dt):
    """Burn deuterium over the step in shells ibegin..iend (one shell, or a
    zone mixed to its mass-weighted means), with the rate rising linearly
    from its start-of-step to its end-of-step value in substeps; update H
    and 3He, and at the surface with accretion burn the accreted D too.

    Translates: dburn.f#DBURN (09c3b9b)

    Arguments: ibegin, iend, the shells; m, the number of shells; hs2, the
    shell masses; hcomp[species, shell], changed in place (species 1, 4,
    12: H, 3He, D); dt, the step (Gyr).

    Returns nothing.
    """
    global summ
    hcompp = oldmod_mod.hcompp
    drate = deuter_mod.drate
    drate0 = deuter_mod.drate0
    if ibegin == iend:
        # Question (see the page): summ is not set here, but the accretion
        # branch below uses it.
        x = hcompp[1, ibegin]
        xh2 = hcompp[12, ibegin]
        xhe3 = hcompp[4, ibegin]
        rate0 = drate0[ibegin]
        rate1 = drate[ibegin]
    else:
        summ = 0.0
        sumr0 = 0.0
        sumr1 = 0.0
        x = 0.0
        xh2 = 0.0
        xhe3 = 0.0
        for i in range(ibegin, iend + 1):
            summ = summ + hs2[i]
            sumr0 = sumr0 + hs2[i] * drate0[i]
            sumr1 = sumr1 + hs2[i] * drate[i]
            x = x + hcompp[1, i] * hs2[i]
            xh2 = xh2 + hcompp[12, i] * hs2[i]
            xhe3 = xhe3 + hcompp[4, i] * hs2[i]
        rate0 = sumr0 / summ
        rate1 = sumr1 / summ
        x = x / summ
        xh2 = xh2 / summ
        xhe3 = xhe3 / summ

    fmassacc = deuter_mod.fmassacc
    compacc = nmlpar_mod.compacc
    if nmlpar_mod.lmdot and iend == m and nmlpar_mod.dmdt0 > 0.0:
        xh2test = (xh2 * summ + compacc[12] * fmassacc) / (summ + fmassacc)
    else:
        xh2test = xh2
    if xh2test < 1.0e-11:
        for i in range(ibegin, iend + 1):
            hcomp[12, i] = 0.0
        return  # label 9999

    # Burn in substeps, the rate rising linearly from rate0 to rate1.
    # Question (see the page): the comment says 0.5 in the natural log per
    # substep; the code limits rate*tstep to 0.1.
    rate = max(rate0, rate1)
    istep = int(rate * dt / 0.1) + 1
    fstep = float(istep)
    tstep = dt / fstep
    r0 = rate0
    delrate = (rate1 - rate0) / fstep
    xh2no = xh2
    for i in range(1, istep + 1):
        rate = r0 + 0.5 * delrate
        xh2no = xh2no * math.exp(-2.0 * rate * tstep)
        r0 = r0 + delrate

    # Question (see the page): unlike the test above, this branch does not
    # check dmdt0 > 0, so under mass loss it uses a stale fmassacc.
    if nmlpar_mod.lmdot and iend == m:
        # Accreted matter: earlier accretion burns fully, new accretion is
        # exposed to half a substep.
        r0 = rate0
        xh2nacc = compacc[12]
        xh2na = 0.0
        for i in range(1, istep + 1):
            rate = r0 + 0.5 * delrate
            xh2na = xh2na * math.exp(-2.0 * rate * tstep)
            x2new = xh2nacc / fstep * math.exp(-rate * tstep)
            xh2na = xh2na + x2new
            r0 = r0 + delrate
        # Mass-weighted means of the original and the accreted matter.
        dxh2o = xh2no - xh2
        dxh2a = xh2na - compacc[12]
        xh2n = (xh2no * summ + xh2na * fmassacc) / (summ + fmassacc)
        xn = x + 0.5 * (dxh2o * summ + dxh2a * fmassacc) / (summ + fmassacc)
        xhe3n = xhe3 - 1.5 * (dxh2o * summ + dxh2a * fmassacc) / (summ + fmassacc)
    else:
        dxh2 = xh2no - xh2
        xh2n = xh2no
        xn = x + 0.5 * dxh2
        xhe3n = xhe3 - 1.5 * dxh2

    for i in range(ibegin, iend + 1):
        hcomp[1, i] = xn
        hcomp[4, i] = xhe3n
        hcomp[12, i] = xh2n
