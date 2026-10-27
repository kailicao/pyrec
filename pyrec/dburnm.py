"""Python version of dburnm.f (YREC 09c3b9b): deuterium burning over a
substep of a model in the middle of a time step (MIDMOD)."""

import math

from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod


def dburnm(ibegin, iend, m, hs2, hcomp, dt, dratem, dratem0, fac2):
    """As DBURN, for a substep: burn deuterium in shells ibegin..iend with
    rates interpolated to the substep, and the substep's share fac2 of the
    accreted matter at the surface.

    Translates: dburnm.f#DBURNM (09c3b9b)

    Arguments: ibegin, iend, the shells; m, the number of shells; hs2, the
    shell masses; hcomp[species, shell], changed in place (H, 3He, D); dt,
    the substep (s); dratem, dratem0, the rates at the end and start of the
    substep; fac2, the substep's fraction of the step.

    Returns nothing.
    """
    hcompp = oldmod_mod.hcompp
    dtt = dt * 1.0e-9 / constants_mod.csecyr
    if ibegin == iend:
        x = hcompp[1, ibegin]
        xh2 = hcompp[12, ibegin]
        xhe3 = hcompp[4, ibegin]
        rate0 = dratem0[ibegin]
        rate1 = dratem[ibegin]
    else:
        summ = 0.0
        sumr0 = 0.0
        sumr1 = 0.0
        x = 0.0
        xh2 = 0.0
        xhe3 = 0.0
        for i in range(ibegin, iend + 1):
            summ = summ + hs2[i]
            sumr0 = sumr0 + hs2[i] * dratem0[i]
            sumr1 = sumr1 + hs2[i] * dratem[i]
            x = x + hcompp[1, i] * hs2[i]
            xh2 = xh2 + hcompp[12, i] * hs2[i]
            xhe3 = xhe3 + hcompp[4, i] * hs2[i]
        rate0 = sumr0 / summ
        rate1 = sumr1 / summ
        x = x / summ
        xh2 = xh2 / summ
        xhe3 = xhe3 / summ

    if xh2 < 1.0e-14:
        for i in range(ibegin, iend + 1):
            hcomp[12, i] = 0.0
        return  # label 9999

    rate = max(rate0, rate1)
    istep = int(rate * dtt / 0.1) + 1
    fstep = float(istep)
    tstep = dtt / fstep
    r0 = rate0
    delrate = (rate1 - rate0) / fstep
    xh2no = xh2
    for i in range(1, istep + 1):
        rate = r0 + 0.5 * delrate
        xh2no = xh2no * math.exp(-2.0 * rate * tstep)
        r0 = r0 + delrate

    compacc = nmlpar_mod.compacc
    if nmlpar_mod.lmdot and iend == m:
        r0 = rate0
        xh2nacc = compacc[12]
        xh2na = 0.0
        for i in range(1, istep + 1):
            rate = r0 + 0.5 * delrate
            xh2na = xh2na * math.exp(-2.0 * rate * tstep)
            x2new = xh2nacc / fstep * math.exp(-rate * tstep)
            xh2na = xh2na + x2new
            r0 = r0 + delrate
        fmassa = fac2 * deuter_mod.fmassacc
        dxh2o = xh2no - xh2
        dxh2a = xh2na - compacc[12]
        # Question (see the page): fmassa is a mass in grams, weighted
        # here as if it were a fraction of the zone's mass.
        xh2n = (xh2no + xh2na * fmassa) / (1.0 + fmassa)
        xn = x + 0.5 * (dxh2o + dxh2a * fmassa) / (1.0 + fmassa)
        xhe3n = xhe3 - 1.5 * (dxh2o + dxh2a * fmassa) / (1.0 + fmassa)
    else:
        dxh2 = xh2no - xh2
        xh2n = xh2no
        xn = x + 0.5 * dxh2
        xhe3n = xhe3 - 1.5 * dxh2

    for i in range(ibegin, iend + 1):
        hcomp[1, i] = xn
        hcomp[4, i] = xhe3n
        hcomp[12, i] = xh2n
