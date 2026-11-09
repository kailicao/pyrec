"""Python version of kemcom.f (YREC development, 7be4133): the implicit
solution of the CNO and pp burning network for one shell or one convection
zone."""

import numpy as np

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec.simeqc import simeqc

# KEMCOM has SAVE: its work arrays are module-level. A(1..49) holds the
# Jacobian column by column, A(50..56) the residuals (once a separate
# array B, EQUIVALENCEd to it), which SIMEQC replaces by the Newton correction.
a = np.zeros(57)
w = np.zeros(8)
xsum = np.zeros(12)


def kemcom(ht: np.ndarray, ibegin: int, iend: int, hr1: np.ndarray, hr2: np.ndarray, hr3: np.ndarray, hr4: np.ndarray,
           hr5: np.ndarray, hr6: np.ndarray, hr7: np.ndarray, hr8: np.ndarray, hr10: np.ndarray, hr11: np.ndarray,
           hr12: np.ndarray, hf1: np.ndarray, hs2: np.ndarray, hcomp: np.ndarray, ddage: float) -> None:
    """New abundances of H, 3He, 4He, 12C, 13C, 14N and 16O after one step.

    Solves the network implicitly by Newton iteration (VandenBerg's
    notes), for one radiative shell (IBEGIN = IEND) or one convection zone
    (mass-weighted means of abundances and rates); 18O is updated
    explicitly afterwards.

    Translates: kemcom.f#KEMCOM (7be4133)

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

    # Abundances at the start of the step, mixed over a zone (Fortran lines
    # 97–119).
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

    # Too cool to burn: homogenize only (Fortran lines 121–128).
    if ht[ibegin] < nml.tcut[1]:
        for i in range(ibegin, iend + 1):                # labels 13, 11
            for j in range(1, 12):
                hcomp[j, i] = xsum[j]
        return                                           # GOTO 200

    # The step in Gyr, times the species' masses (Fortran lines 134–143).
    delt = ddage * 1.0e-9
    delt3 = 3.0 * delt
    delt4 = 4.0 * delt
    delt12 = 12.0 * delt
    delt13 = 13.0 * delt
    delt14 = 14.0 * delt
    delt16 = 16.0 * delt
    dcmin = nml.cmin
    icon = 0

    # Fortran lines 148–149.
    gr9 = 0.0  # noqa: F841  (unused in the Fortran too)
    gr13 = 0.0  # noqa: F841  (unused in the Fortran too)

    # The rates: the shell's, or the zone's mass-weighted means, summed in
    # the Fortran's order (Fortran lines 150–223).
    if ibegin == iend:
        k = ibegin
        gr1, gr2, gr3, gr4 = hr1[k], hr2[k], hr3[k], hr4[k]
        gr5, gr6, gr7, gr8 = hr5[k], hr6[k], hr7[k], hr8[k]
        gr10, gr11, gr12 = hr10[k], hr11[k], hr12[k]
        f3 = hf1[ibegin]
        f4 = 1.0 - f3
    else:
        gr1 = gr2 = gr3 = gr4 = gr5 = gr6 = gr7 = gr8 = 0.0
        gr10 = gr11 = gr12 = 0.0
        f3 = 0.0
        for i in range(ibegin, iend + 1):                # label 15
            gr1 = gr1 + hs2[i] * hr1[i]
            gr2 = gr2 + hs2[i] * hr2[i]
            gr3 = gr3 + hs2[i] * hr3[i]
            gr4 = gr4 + hs2[i] * hr4[i]
            gr5 = gr5 + hs2[i] * hr5[i]
            gr6 = gr6 + hs2[i] * hr6[i]
            gr7 = gr7 + hs2[i] * hr7[i]
            gr8 = gr8 + hs2[i] * hr8[i]
            gr10 = gr10 + hs2[i] * hr10[i]
            gr11 = gr11 + hs2[i] * hr11[i]
            gr12 = gr12 + hs2[i] * hr12[i]
            f3 = f3 + hs2[i] * hf1[i]
        gr1 = gr1 / wtsum
        gr2 = gr2 / wtsum
        gr3 = gr3 / wtsum
        gr4 = gr4 / wtsum
        gr5 = gr5 / wtsum
        gr6 = gr6 / wtsum
        gr7 = gr7 / wtsum
        gr8 = gr8 / wtsum
        gr10 = gr10 / wtsum
        gr11 = gr11 / wtsum
        gr12 = gr12 / wtsum
        f3 = f3 / wtsum
        f4 = 1.0 - f3

    # Starting abundances (note the order of w: X, 3He, Y, 12C, 13C, 14N,
    # 16O = HCOMP elements 1, 4, 2, 5, 6, 7, 9) (Fortran lines 227–233).
    xpr = xsum[1]
    y3r = xsum[4]
    ypr = xsum[2]
    x2r = xsum[5]
    x3r = xsum[6]
    x4r = xsum[7]
    x6r = xsum[9]

    # First guess: the shell below, or the shell itself for the first one
    # (Fortran lines 238–254).
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
        # Residuals G_k (Fortran lines 264–284).
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

        # The Jacobian, column by column (A(1..7) is d G / d X, ...)
        # (Fortran lines 290–346).
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
        # 1); solve only if any off-diagonal entry is left (Fortran lines
        # 349–360).
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

        # Newton update (Fortran lines 362–373).
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

        # The convergence test (Fortran lines 374–389).
        if bigd >= nml.abstol or bigb >= nml.reltol:
            icon = icon + 1
            if icon >= nml.kemmax:
                # Suspected bug (see the page): the message says 50, the
                # limit is KEMMAX.
                iounits_mod.write(iounits_mod.ishort,
                                  " " + ">" * 39 + "<" * 40 + "\n"
                                  " ERROR IN SUBROUTINE KEMCOM\n"
                                  f" UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL{ibegin:4d}\n"
                                  " RUN STOPPED AFTER 50 ATTEMPTS")  # format 1000
                raise SystemExit(1)    # STOP 1, the run-time error code
            continue                                     # GOTO 10
        break

    # 18O explicitly (Fortran line 393). Question (see the page): production
    # only.
    xo18 = xsum[11] + 18.0 * delt * gr11 * w[3] * w[6]

    # Z: reset, or the zone's mean (Fortran lines 395–399).
    # Question (see the page): the threshold 5e-7 is hard-coded.
    if hcompp[1, iend] < 5.0e-7:
        z = 1.0 - w[1] - w[2] - w[3]
    else:
        z = xsum[3]

    # The new abundances for every shell of the zone (Fortran lines
    # 400–410). Question (see the page): elements 8 and 10 are not mixed.
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
