"""Python version of kemcom.f (YREC development, 09c3b9b): the implicit
solution of the CNO and pp burning network for one shell or one convection
zone."""

import numpy as np

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec.simeqc import simeqc

# KEMCOM has SAVE: its work arrays are module-level. A(1..49) holds the
# Jacobian column by column, A(50..56) the residuals (the B of rp_1's
# EQUIVALENCE), which SIMEQC replaces by the Newton correction.
a = np.zeros(57)
w = np.zeros(8)
xsum = np.zeros(12)


def kemcom(ht, ibegin, iend, hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8,
           hr10, hr11, hr12, hf1, hs2, hcomp, ddage):
    """New abundances of H, 3He, 4He, 12C, 13C, 14N and 16O after one step.

    Solves the network implicitly by Newton iteration (VandenBerg's
    notes), for one radiative shell (IBEGIN = IEND) or one convection zone
    (mass-weighted means of abundances and rates); 18O is updated
    explicitly afterwards.

    Translates: kemcom.f#KEMCOM (09c3b9b)

    Parameters
    ----------
    ht : array
        log T.
    ibegin, iend : int
        First and last shell of the zone.
    hr1 ... hr12, hf1 : arrays
        Reaction rates (per Gyr per amu, without composition factors) and
        the 15N + p branching ratio.
    hs2 : array
        Shell masses.
    hcomp : array
        Abundances; elements 1-7, 9 and 11 set for shells IBEGIN..IEND.
    ddage : float
        Time step (years).
    """
    nml = nmlpar_mod
    hcompp = oldmod_mod.hcompp
    # Abundances at the start of the step, mixed over a zone.
    if ibegin != iend:
        wtsum = 0.0
        for j in range(1, 12):                           # label 1
            xsum[j] = 0.0
        for i in range(ibegin, iend + 1):                # labels 5, 3
            wtsum = wtsum + hs2[i]
            for j in range(1, 12):
                xsum[j] = xsum[j] + hcompp[j, i] * hs2[i]
        for j in range(1, 12):                           # label 7
            xsum[j] = xsum[j] / wtsum
    else:
        for j in range(1, 12):                           # label 9
            xsum[j] = hcompp[j, ibegin]
    # Too cool to burn: homogenize only.
    if ht[ibegin] < nml.tcut[1]:
        for i in range(ibegin, iend + 1):                # labels 13, 11
            for j in range(1, 12):
                hcomp[j, i] = xsum[j]
        return                                           # GOTO 200

    delt = ddage * 1.0e-9
    delt3 = 3.0 * delt
    delt4 = 4.0 * delt
    delt12 = 12.0 * delt
    delt13 = 13.0 * delt
    delt14 = 14.0 * delt
    delt16 = 16.0 * delt
    dcmin = nml.cmin
    icon = 0
    gr9 = 0.0
    gr13 = 0.0
    # The rates: the shell's, or the zone's mass-weighted means.
    if ibegin == iend:
        k = ibegin
        gr1, gr2, gr3, gr4 = hr1[k], hr2[k], hr3[k], hr4[k]
        gr5, gr6, gr7, gr8 = hr5[k], hr6[k], hr7[k], hr8[k]
        gr10, gr11, gr12 = hr10[k], hr11[k], hr12[k]
        f3 = hf1[ibegin]
        f4 = 1.0 - f3
    else:
        sl = slice(ibegin, iend + 1)
        wm = hs2[sl]
        gr1 = np.sum(wm * hr1[sl]) / wtsum               # label 15
        gr2 = np.sum(wm * hr2[sl]) / wtsum
        gr3 = np.sum(wm * hr3[sl]) / wtsum
        gr4 = np.sum(wm * hr4[sl]) / wtsum
        gr5 = np.sum(wm * hr5[sl]) / wtsum
        gr6 = np.sum(wm * hr6[sl]) / wtsum
        gr7 = np.sum(wm * hr7[sl]) / wtsum
        gr8 = np.sum(wm * hr8[sl]) / wtsum
        gr10 = np.sum(wm * hr10[sl]) / wtsum
        gr11 = np.sum(wm * hr11[sl]) / wtsum
        gr12 = np.sum(wm * hr12[sl]) / wtsum
        f3 = np.sum(wm * hf1[sl]) / wtsum
        f4 = 1.0 - f3
    # Starting abundances (note the order of w: X, 3He, Y, 12C, 13C, 14N,
    # 16O = HCOMP elements 1, 4, 2, 5, 6, 7, 9).
    xpr = xsum[1]
    y3r = xsum[4]
    ypr = xsum[2]
    x2r = xsum[5]
    x3r = xsum[6]
    x4r = xsum[7]
    x6r = xsum[9]
    # First guess: the shell below, or the shell itself for the first one.
    # Question (see the page): the comment's ITLVL no longer exists.
    k = ibegin - 1 if ibegin != 1 else ibegin
    w[1] = hcomp[1, k]
    w[2] = hcomp[4, k]
    w[3] = hcomp[2, k]
    w[4] = hcomp[5, k]
    w[5] = hcomp[6, k]
    w[6] = hcomp[7, k]
    w[7] = hcomp[9, k]

    while True:                                          # label 10
        # Residuals G_k.
        a[50] = w[1] - xpr - delt * (-3.0 * gr1 * w[1] ** 2 + 2.0 * gr2 * w[2] ** 2
                                     - gr3 * w[2] * w[3] - gr4 * w[1] * w[4]
                                     - gr5 * w[1] * w[5]
                                     - 2.0 * w[1] * (gr6 * w[6] + gr7 * w[7]))
        a[51] = w[2] - y3r - delt3 * (gr1 * w[1] ** 2 - 2.0 * gr2 * w[2] ** 2
                                      - gr3 * w[2] * w[3])
        a[52] = w[3] - ypr - delt4 * (gr2 * w[2] ** 2 + gr3 * w[2] * w[3]
                                      + f3 * gr6 * w[6] * w[1] + gr7 * w[1] * w[7]
                                      - gr10 * w[3] * w[4] - gr8 * w[3] * w[5]
                                      - gr11 * w[3] * w[6]
                                      - 3.0 * gr12 * w[3] ** 3)
        a[53] = w[4] - x2r - delt12 * (gr12 * w[3] ** 3 - gr4 * w[1] * w[4]
                                       - gr10 * w[3] * w[4]
                                       + f3 * gr6 * w[1] * w[6])
        a[54] = w[5] - x3r - delt13 * (gr4 * w[1] * w[4] - gr5 * w[1] * w[5]
                                       - gr8 * w[3] * w[5])
        a[55] = w[6] - x4r - delt14 * (gr5 * w[1] * w[5] + gr7 * w[1] * w[7]
                                       - gr6 * w[1] * w[6] - gr11 * w[3] * w[6])
        a[56] = w[7] - x6r - delt16 * (f4 * gr6 * w[1] * w[6] + gr10 * w[3] * w[4]
                                       - gr7 * w[1] * w[7] + gr8 * w[3] * w[5])
        # The Jacobian, column by column (A(1..7) is d G / d X, ...).
        a[1] = 1.0 + delt * (6.0 * gr1 * w[1] + gr4 * w[4] + gr5 * w[5]
                             + 2.0 * (gr6 * w[6] + gr7 * w[7]))
        a[2] = -2.0 * delt3 * gr1 * w[1]
        a[3] = -delt4 * (f3 * gr6 * w[6] + gr7 * w[7])
        a[4] = -delt12 * (f3 * gr6 * w[6] - gr4 * w[4])
        a[5] = -delt13 * (gr4 * w[4] - gr5 * w[5])
        a[6] = -delt14 * (gr5 * w[5] + gr7 * w[7] - gr6 * w[6])
        a[7] = -delt16 * (f4 * gr6 * w[6] - gr7 * w[7])

        a[8] = -delt * (4.0 * gr2 * w[2] - gr3 * w[3])
        a[9] = 1.0 + delt3 * (4.0 * gr2 * w[2] + gr3 * w[3])
        a[10] = -delt4 * (2.0 * gr2 * w[2] + gr3 * w[3])
        a[11] = 0.0
        a[12] = 0.0
        a[13] = 0.0
        a[14] = 0.0

        a[15] = delt * gr3 * w[2]
        a[16] = delt3 * gr3 * w[2]
        a[17] = 1.0 - delt4 * (gr3 * w[2] - gr10 * w[4] - gr8 * w[5]
                               - gr11 * w[6] - 9.0 * gr12 * w[3] ** 2)
        a[18] = -delt12 * (3.0 * gr12 * w[3] ** 2 - gr10 * w[4])
        a[19] = delt13 * gr8 * w[5]
        a[20] = delt14 * gr11 * w[6]
        a[21] = -delt16 * (gr10 * w[4] + gr8 * w[5])

        a[22] = delt * gr4 * w[1]
        a[23] = 0.0
        a[24] = delt4 * gr10 * w[3]
        a[25] = 1.0 + delt12 * (gr4 * w[1] + gr10 * w[3])
        a[26] = -delt13 * gr4 * w[1]
        a[27] = 0.0
        a[28] = -delt16 * gr10 * w[3]

        a[29] = delt * gr5 * w[1]
        a[30] = 0.0
        a[31] = delt4 * gr8 * w[3]
        a[32] = 0.0
        a[33] = 1.0 + delt13 * (gr5 * w[1] + gr8 * w[3])
        a[34] = -delt14 * gr5 * w[1]
        a[35] = -delt16 * gr8 * w[3]

        a[36] = 2.0 * delt * gr6 * w[1]
        a[37] = 0.0
        a[38] = -delt4 * (f3 * gr6 * w[1] - gr11 * w[3])
        a[39] = -delt12 * f3 * gr6 * w[1]
        a[40] = 0.0
        a[41] = 1.0 + delt14 * (gr6 * w[1] + gr11 * w[3])
        a[42] = -f4 * delt16 * gr6 * w[1]

        a[43] = 2.0 * delt * gr7 * w[1]
        a[44] = 0.0
        a[45] = -delt4 * gr7 * w[1]
        a[46] = 0.0
        a[47] = 0.0
        a[48] = -delt14 * gr7 * w[1]
        a[49] = 1.0 + delt16 * gr7 * w[1]

        # Zero small off-diagonal entries (the diagonal is every 8th, from
        # 1); solve only if any off-diagonal entry is left.
        j1j = 0
        j2j = 1
        for i in range(1, 50):                           # label 20
            if i == j2j:
                j2j = j2j + 8
            elif abs(a[i]) < 1.0e-8:
                a[i] = 0.0
            else:
                j1j = 1
        # Question (see the page): without SIMEQC the correction is the
        # residual itself, not divided by the diagonal.
        if j1j != 0:
            simeqc(a, 8, 7)
        # Newton update and the convergence test.
        bigb = 0.0
        bigd = 0.0
        for i in range(1, 8):                            # label 30
            w[i] = w[i] - a[49 + i]
            if abs(a[49 + i]) >= bigd:
                bigd = abs(a[49 + i])
            if w[i] < dcmin:
                w[i] = 0.0
            else:
                xb = abs(a[49 + i] / w[i])
                if xb >= bigb:
                    bigb = xb
        if bigd >= nml.abstol or bigb >= nml.reltol:
            icon = icon + 1
            if icon >= nml.kemmax:
                # Question (see the page): the message says 50, the limit
                # is KEMMAX.
                iounits_mod.write(iounits_mod.ishort,
                                  " " + ">" * 39 + "<" * 40 + "\n"
                                  " ERROR IN SUBROUTINE KEMCOM\n"
                                  f" UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL{ibegin:4d}\n"
                                  " RUN STOPPED AFTER 50 ATTEMPTS")  # format 1000
                raise SystemExit
            continue                                     # GOTO 10
        break

    # 18O explicitly. Question (see the page): production only.
    xo18 = xsum[11] + 18.0 * delt * gr11 * w[3] * w[6]
    # Question (see the page): the threshold 5e-7 is hard-coded.
    if hcompp[1, iend] < 5.0e-7:
        z = 1.0 - w[1] - w[2] - w[3]
    else:
        z = xsum[3]
    for iu in range(ibegin, iend + 1):                   # label 40
        hcomp[1, iu] = w[1]
        hcomp[4, iu] = w[2]
        hcomp[2, iu] = w[3]
        hcomp[5, iu] = w[4]
        hcomp[6, iu] = w[5]
        hcomp[7, iu] = w[6]
        hcomp[9, iu] = w[7]
        hcomp[3, iu] = z
        hcomp[11, iu] = xo18
    # label 200
