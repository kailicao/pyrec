"""Python version of eqburn.f (YREC 09c3b9b): explicit hydrogen and helium
burning rates for a shell or a convection zone, with 3He and the CN cycle
in equilibrium."""

import math

import numpy as np

from pyrec import nmlpar_mod
from pyrec import oldmod_mod

# SAVE: the work array XSUM(11), set before use on every call.
xsum = np.zeros(11 + 1)


def eqburn(hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr10, hr12, hs2, ht,
           ibegin, iend, dcdt, dodt, dxdt, dydt, xc12, xo16, x, z):
    """Explicit burning rates for shells ibegin..iend (mass-weighted means
    over a zone), with 3He and the CN cycle in equilibrium: an initial
    guess at the rates at the start of the step.

    Translates: eqburn.f#EQBURN (09c3b9b)

    Arguments: hr1..hr12, the reaction-rate arrays from RATES; hs2, the
    shell masses; ht, log T; ibegin, iend, the shell or zone; dcdt, dodt,
    dxdt, dydt, xc12, xo16, x, z (set only, except below tcut[1], where
    xc12, xo16, x, z keep their values).

    Returns (dcdt, dodt, dxdt, dydt, xc12, xo16, x, z): the rates of change
    of 12C, 16O, H and He (mass fraction per Gyr) and the zone's mean
    12C (its CN-equilibrium value if the CN step ran), 16O, X and Z.
    """
    tcut = nmlpar_mod.tcut
    hcompp = oldmod_mod.hcompp
    dxdt = 0.0
    dydt = 0.0
    dcdt = 0.0
    dodt = 0.0
    if ht[ibegin] < tcut[1]:
        return dcdt, dodt, dxdt, dydt, xc12, xo16, x, z

    # Abundances: mass-weighted means over a zone.
    if ibegin != iend:
        wtsum = 0.0
        for j in range(1, 12):  # label 1
            xsum[j] = 0.0
        for i in range(ibegin, iend + 1):  # label 5
            wtsum = wtsum + hs2[i]
            for j in range(1, 12):  # label 3
                xsum[j] = xsum[j] + hcompp[j, i] * hs2[i]
        for j in range(1, 12):  # label 7
            xsum[j] = xsum[j] / wtsum
    else:
        for j in range(1, 12):  # label 9
            xsum[j] = hcompp[j, ibegin]

    # Rates: likewise.
    if ibegin == iend:
        gr1 = hr1[ibegin]
        gr2 = hr2[ibegin]
        gr3 = hr3[ibegin]
        gr4 = hr4[ibegin]
        gr5 = hr5[ibegin]
        gr6 = hr6[ibegin]
        gr7 = hr7[ibegin]
        gr10 = hr10[ibegin]
        gr12 = hr12[ibegin]
    else:
        gr1 = gr2 = gr3 = gr4 = gr5 = gr6 = gr7 = gr10 = gr12 = 0.0
        for i in range(ibegin, iend + 1):  # label 15
            gr1 = gr1 + hs2[i] * hr1[i]
            gr2 = gr2 + hs2[i] * hr2[i]
            gr3 = gr3 + hs2[i] * hr3[i]
            gr4 = gr4 + hs2[i] * hr4[i]
            gr5 = gr5 + hs2[i] * hr5[i]
            gr6 = gr6 + hs2[i] * hr6[i]
            gr7 = gr7 + hs2[i] * hr7[i]
            gr10 = gr10 + hs2[i] * hr10[i]
            gr12 = gr12 + hs2[i] * hr12[i]
        gr1 = gr1 / wtsum
        gr2 = gr2 / wtsum
        gr3 = gr3 / wtsum
        gr4 = gr4 / wtsum
        gr5 = gr5 / wtsum
        gr6 = gr6 / wtsum
        gr7 = gr7 / wtsum
        gr10 = gr10 / wtsum
        gr12 = gr12 / wtsum

    x = xsum[1]
    xhe3 = xsum[4]
    y = xsum[2]
    z = xsum[3]
    xc12 = xsum[5]
    xc13 = xsum[6]
    xn14 = xsum[7]
    xo16 = xsum[9]

    # pp chains, 3He in equilibrium (capped at its actual abundance).
    # Question (see the page): divides by gr2 (and gr4, gr5 below).
    if x > 1.0e-10 and ht[ibegin] > tcut[2]:
        r11 = gr1 * x**2
        r34 = gr3 * y
        xhe3e = 0.25 * (math.sqrt(r34**2 + 8.0 * gr2 * r11) - r34) / gr2
        xhe3 = min(xhe3, xhe3e)
        dxdt = -3.0 * r11 + 2.0 * gr2 * xhe3**2 - r34 * xhe3
        # CN cycle in equilibrium.
        if ht[ibegin] > tcut[3]:
            sumcno = xc12 / 12.0 + xc13 / 13.0 + xn14 / 14.0
            r1214 = gr6 / gr4
            r1314 = gr6 / gr5
            xn14 = sumcno / (r1214 / 12.0 + r1314 / 13.0 + 1.0 / 14.0)
            # Question (see the page): xc12 becomes its equilibrium value,
            # returned and used for 12C + alpha below.
            xc12 = xn14 * r1214
            xc13 = xn14 * r1314
            # Question (see the page): two protons for 16O + p.
            dxdt = (dxdt - gr4 * x * xc12 - gr5 * x * xc13
                    - 2.0 * gr6 * x * xn14 - 2.0 * gr7 * x * xo16)
    else:
        dxdt = 0.0  # GOTO 100: the CN step is skipped

    # label 100: helium burning.
    if y > 1.0e-10 and ht[ibegin] > tcut[4]:
        r3a = gr12 * y**3
        r412 = gr10 * y * xc12
        dydt = -4.0 * (r412 + 3.0 * r3a)
        dcdt = 12.0 * (r3a - r412)
        dodt = 16.0 * r412
    else:
        dydt = 0.0
        dcdt = 0.0
        dodt = 0.0

    return dcdt, dodt, dxdt, dydt, xc12, xo16, x, z
