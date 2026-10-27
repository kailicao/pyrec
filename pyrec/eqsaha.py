"""Python version of eqsaha.f (YREC development, 09c3b9b): the Saha
equation of state for partially ionized gas."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import iounits_mod

# DATA: ionization temperatures chi/k (K) of the 11 species, in the order
# Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne (index 0 unused), and of He, He+.
SAHAT = np.array([0.0, 59630., 69450., 88710., 91330., 94570., 130630.,
                  157800., 158000., 168630., 183080., 250200.])
SAHET0 = 285270.
SAHET1 = 631370.
# DATA SAHAU: the constant u_i of each Saha equation (13 values: the 11
# species, then He and He+). Double precision since 09c3b9b.
SAHAU = np.array([0.0, -.01, -.47, .60, .49, .08, .11, .0, -.04, .64, .99,
                  1.03, .6, .0])
NZ = 11
MAXS = 20
STOL = 4.0e1
STOLEE = 1.0e-09

# DATA E,SAHEX1,SAHEX2 and SAVE: the electrons per atom and the He+ and
# He++ fractions, kept between calls as the next starting guess.
# Question (see the page): the result can depend on the order of calls.
e = 0.1
sahex1 = 0.0
sahex2 = 0.0


def eqsaha(fx, tl, t, pl, p, dl, d, beta, betai, beta14, fxion, rmu,
           amu, emu, lderiv, latmo, qdt, qdp, qcp, dela, qdtt, qdtp, qat,
           qap, qcpt, qcpp, ksaha):
    """Solve the Saha equations for hydrogen, helium (two stages) and ten
    metals by Newton iteration on the electrons per atom, and return the
    density, the ionization fractions and the thermodynamic derivatives.

    Translates: eqsaha.f#EQSAHA (09c3b9b)

    Arguments: fx, the number fractions of the 12 species (index 0
    unused); tl, t, pl, p (log T, T, log P, P); beta and its relatives
    betai (1/beta) and beta14; amu; lderiv (also the second derivatives);
    latmo (skip the derivatives); ksaha, the caller's iteration counter.
    Array changed in place: fxion (H+, He+, He++ fractions, index 0
    unused).

    Returns (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
    qcpt, qcpp, ksaha). With latmo, qdt ... qcpp are returned as passed in;
    without lderiv, qdtt ... qcpp are.
    """
    global e, sahex1, sahex2
    c = constants_mod
    sahatt = np.zeros(NZ + 1)
    sahak = np.zeros(NZ + 1)
    sahax = np.zeros(NZ + 1)
    fxs = np.zeros(NZ + 1)
    fxt = np.zeros(NZ + 1)

    fxhe = fx[12]
    emax = 1.0 + fxhe
    ti = 1.0 / t
    # The Saha K's, in the order of increasing ionization temperature. The
    # first species with an exponent below -STOL ends the list: it and all
    # later ones count as neutral. Above +STOL, the species is fully
    # ionized (K = 1e16) and left out of the derivative sums (nz1).
    stemp = 2.50 * tl - pl + c.cmkh
    nz1 = 1
    for i in range(1, NZ + 1):
        sahatt[i] = SAHAT[i] * ti
        sahak[i] = c.cln * (SAHAU[i] + stemp) - sahatt[i]
        if sahak[i] < -STOL:
            break                                  # go to 12
        if sahak[i] > +STOL:
            sahak[i] = 1.0e16                      # label 10
            nz1 = i + 1
            continue
        sahak[i] = betai * math.exp(sahak[i])
    else:
        i = NZ + 1
    nz0 = i - 1                                    # label 12
    nz1 = min(NZ, nz1)
    # Helium: neutral (lnohe), or no He++ (lnohe2).
    lnohe = True
    lnohe2 = True
    sahek0 = sahek1 = 0.0
    if fxhe < 1.0e-10:
        sahex1 = 0.0                               # label 13
        sahex2 = 0.0                               # label 14
    else:
        sahek0 = c.cln * (SAHAU[12] + stemp) - SAHET0 / t
        if sahek0 < -STOL:
            sahex1 = 0.0
            sahex2 = 0.0
        else:
            lnohe = False
            sahek0 = betai * math.exp(sahek0)
            sahek1 = c.cln * (SAHAU[13] + stemp) - SAHET1 / t
            if sahek1 < -STOL:
                sahex2 = 0.0
            else:
                lnohe2 = False
                sahek1 = betai * math.exp(sahek1)

    # Newton iterations for e, sahex1, sahex2 (the 3x3 system reduced by
    # hand).
    cr12 = cr13 = cr23 = 0.0
    c11 = c12 = c13 = c22 = c23 = 0.0
    r1 = r2 = 0.0
    for isaha in range(1, MAXS + 1):
        ksaha = ksaha + 1
        ep1 = 1.0 / (1.0 + e)
        ep1e = ep1 / e
        ep12 = ep1 * ep1
        eep1 = e * ep1
        c33 = -e * (1.0 + e)
        r3 = e
        for i in range(1, nz0 + 1):
            div = 1.0 / (sahak[i] + eep1)
            temp8 = sahak[i] * div
            sahax[i] = temp8
            temp8 = temp8 * fx[i]
            r3 = r3 - temp8
            temp8 = temp8 * eep1 * div
            c33 = c33 - temp8
            fxs[i] = temp8
        c33 = c33 * ep1e                           # label 22
        r3 = r3 - fxhe * (sahex1 + sahex2 + sahex2)
        lfin = abs(r3) < STOLEE
        if not lnohe:
            c22 = sahek0 + eep1
            c23 = sahex1 * ep12
            c32 = fxhe
            r2 = sahek0 * (1.0 - sahex1 - sahex2) - sahex1 * eep1
            if not lnohe2:
                c11 = eep1
                c12 = -sahek1
                c13 = sahex2 * ep12
                c21 = sahek0
                c31 = 2.0 * fxhe
                r1 = sahex1 * sahek1 - sahex2 * eep1
                # Reduction of the Saha matrix to diagonal form.
                cr12 = -c21 / c11
                cr13 = -c31 / c11
                c22 = c22 + cr12 * c12
                c32 = c32 + cr13 * c12
                c23 = c23 + cr12 * c13
                c33 = c33 + cr13 * c13
                # Products below 1e-38 are skipped (an old underflow guard).
                if r1 != 0.0:
                    r1l = math.log10(abs(r1))
                    if cr12 != 0.0:
                        fact = math.log10(abs(cr12)) + r1l
                        if not fact < -38.0:
                            r2 = r2 + cr12 * r1
                    if cr13 != 0.0:                # label 200
                        fact = math.log10(abs(cr13)) + r1l
                        if not fact < -38.0:
                            r3 = r3 + cr13 * r1
            # label 23: entry for no fully ionized helium
            cr23 = -c32 / c22
            c33 = c33 + cr23 * c23
            r3 = r3 + cr23 * r2
        # label 25: entry for neutral helium
        delte = r3 / c33
        if not lnohe:
            deltx1 = (r2 - c23 * delte) / c22
            sahex1 = sahex1 + deltx1
            sahex1 = max(0.0, min(1.0, sahex1))
            if not lnohe2:
                # Products below 1e-38 are replaced by +-1e-38 ("statement
                # recalculated for DEC-20 system").
                fact1 = c12 * deltx1               # label 100
                if c12 != 0.0 and deltx1 != 0.0:
                    if math.log10(abs(c12)) + math.log10(abs(deltx1)) < -38.0:
                        fact1 = (math.exp(c.cln * -38.0)
                                 * math.copysign(1.0, c12)
                                 * math.copysign(1.0, deltx1))
                fact2 = c13 * delte                # label 110
                if c13 != 0.0 and delte != 0.0:
                    if math.log10(abs(c13)) + math.log10(abs(delte)) < -38.0:
                        fact2 = (math.exp(c.cln * -38.0)
                                 * math.copysign(1.0, c13)
                                 * math.copysign(1.0, delte))
                deltx2 = (r1 - fact1 - fact2) / c11
                sahex2 = sahex2 + deltx2
                sahex2 = max(0.0, min(1.0, sahex2))
        e = e + delte                              # label 26
        e = max(1.0e-11, min(emax, e))
        lfin = lfin and abs(delte) < STOLEE
        if lfin:
            break                                  # go to 29
    else:
        # Question (see the page): the run goes on with the unconverged e.
        iounits_mod.write(iounits_mod.ishort,
                          f" -----SAHA FAILURE (TL,PL)={tl:10.6f}{pl:10.6f}"
                          f"  (E,DE)={e:20.12f}{delte:20.12f}")  # format 99

    # label 29: the system has been solved for e.
    emu = e * amu
    rmu = c.cgas * (emu + amu)
    d = beta * p / (rmu * t)
    dl = math.log10(d)
    fxion[1] = 0.0
    if nz0 >= 7:
        fxion[1] = sahax[7]                        # hydrogen is species 7
    fxion[2] = sahex1
    fxion[3] = sahex2
    if latmo:
        return (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
                qcpt, qcpp, ksaha)                 # go to 60

    # First derivatives: d ln E / d ln T (sqet) and d ln P (sqep).
    rmub = rmu * betai                             # set, never used
    r3t = 0.0
    r3p = 0.0
    if nz0 > 0:
        for i in range(nz1, nz0 + 1):
            r3p = r3p - fxs[i]
            r3t = r3t - fxs[i] * (2.5 + beta14 + sahatt[i])
        r3p = -betai * r3p
    r1t = r1p = r2t = r2p = 0.0
    sahex0 = sk0qt = sk1qt = 0.0
    if not lnohe:                                  # label 31
        sahex0 = 1.0 - sahex1 - sahex2
        sk0qt = 2.5 + beta14 + SAHET0 * ti
        temp = sahek0 * sahex0
        r2t = temp * sk0qt
        r2p = -temp * betai
        if not lnohe2:
            sk1qt = 2.5 + beta14 + SAHET1 * ti
            temp = sahek1 * sahex1
            r1t = temp * sk1qt
            r1p = -temp * betai
            r2t = r2t + cr12 * r1t
            r2p = r2p + cr12 * r1p
            r3t = r3t + cr13 * r1t
            r3p = r3p + cr13 * r1p
        r3t = r3t + cr23 * r2t                     # label 32
        r3p = r3p + cr23 * r2p
    sqet = r3t / c33                               # label 33
    sqep = r3p / c33
    qdt = -1.0 - beta14 - ep1 * sqet
    qdp = betai - ep1 * sqep
    # Internal energy temperature derivative (qut).
    usum = 0.0
    sx0qt = sx0qp = sx1qt = sx1qp = sx2qt = sx2qp = 0.0
    if not lnohe:
        sx1qt = (r2t - c23 * sqet) / c22
        sx1qp = (r2p - c23 * sqep) / c22
        sx2qt = 0.0
        sx2qp = 0.0
        if not lnohe2:
            sx2qt = (r1t - c12 * sx1qt - c13 * sqet) / c11
            sx2qp = (r1p - c12 * sx1qp - c13 * sqep) / c11
        sx0qt = -(sx1qt + sx2qt)                   # label 34
        sx0qp = -(sx1qp + sx2qp)
        usum = fxhe * (sx1qt * SAHET0 + sx2qt * (SAHET0 + SAHET1))
    if nz0 > 0:                                    # label 35
        stemp = 2.50 + beta14 - ep1e * sqet
        for i in range(nz1, nz0 + 1):
            fxt[i] = fxs[i] * (stemp + sahatt[i])
            usum = usum + fxt[i] * SAHAT[i]
    betaut = 0.75 * (2.0 + beta14)                 # label 37
    ramu = c.cgas * amu
    qut = ramu * ((1.0 + e) * (betaut + 3.0 * betai * beta14)
                  + betaut * sqet + usum * ti)
    # Specific heat (qcp) and adiabatic gradient (dela).
    pdt = p / (d * t)
    pdtq = pdt * qdt
    qcp = qut - pdtq
    qcpi = 1.0 / qcp
    dela = -pdtq * qcpi
    if not lderiv:
        return (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
                qcpt, qcpp, ksaha)                 # go to 60

    # Derivatives of qdt, qcp and dela.
    beta16 = 4.0 * betai * beta14
    r3t = 0.0
    r3p = 0.0
    utsum = 0.0
    upsum = 0.0
    if nz0 > 0:
        etemp = (1.0 + e + e) / (e * (1.0 + e))**2
        stemp1 = 2.5 + beta14 - sqet * ep1e
        stemp2 = -betai - sqep * ep1e
        stemp3 = etemp * sqep * sqet - betai * beta14
        stemp4 = etemp * sqet**2 + beta16
        for i in range(nz1, nz0 + 1):
            stemp5 = (stemp1 + sahatt[i]) * (1.0 - 2.0 * sahax[i])
            stemp6 = fxs[i] * (stemp4 - sahatt[i]) + fxt[i] * stemp5
            r3t = r3t - stemp6
            utsum = utsum + SAHAT[i] * (stemp6 - fxt[i])
            stemp7 = fxs[i] * (stemp3 + stemp2 * stemp5)
            r3p = r3p - stemp7
            upsum = upsum + SAHAT[i] * stemp7
    if not lnohe:                                  # label 42
        stemp = 2.0 * sahex1 * ep1 * sqet
        r2t = (sahek0 * (sahex0 * (beta16 - SAHET0 * ti + sk0qt**2)
                         + 2.0 * sx0qt * sk0qt)
               + ep12 * (stemp * sqet - 2.0 * sx1qt * sqet))
        r2p = (sahek0 * (sahex0 * betai * (-beta14 - sk0qt) + sx0qp * sk0qt
                         - betai * sx0qt)
               + ep12 * (stemp * sqep - sx1qp * sqet - sx1qt * sqep))
        if not lnohe2:
            stemp = 2.0 * sahex2 * ep1 * sqet
            r1t = (sahek1 * (sahex1 * (beta16 - SAHET1 * ti + sk1qt**2)
                             + 2.0 * sx1qt * sk1qt)
                   + ep12 * (stemp * sqet - 2.0 * sx2qt * sqet))
            r1p = (sahek1 * (sahex1 * betai * (-beta14 - sk1qt)
                             + sx1qp * sk1qt - betai * sx1qt)
                   + ep12 * (stemp * sqep - sx2qp * sqet - sx2qt * sqep))
            r2t = r2t + cr12 * r1t
            r2p = r2p + cr12 * r1p
            r3p = r3p + cr13 * r1p
            r3t = r3t + cr13 * r1t
        r3p = r3p + cr23 * r2p                     # label 43
        r3t = r3t + cr23 * r2t
    sqett = r3t / c33                              # label 44
    sqetp = r3p / c33
    qdtt = -beta16 + (ep1 * sqet)**2 - ep1 * sqett
    qdtp = betai * beta14 + ep12 * sqet * sqep - ep1 * sqetp
    if nz0 > 0:
        stemp = ep1e * sqett
        stemp1 = ep1e * sqetp
        for i in range(nz1, nz0 + 1):
            utsum = utsum - fxs[i] * SAHAT[i] * stemp
            upsum = upsum - fxs[i] * SAHAT[i] * stemp1
    if not lnohe:                                  # label 46
        sx1qtt = (r2t - c23 * sqett) / c22
        sx1qtp = (r2p - c23 * sqetp) / c22
        if not lnohe2:
            sx2qtt = (r1t - c12 * sx1qtt - c13 * sqett) / c11
            sx2qtp = (r1p - c12 * sx1qtp - c13 * sqetp) / c11
            utsum = utsum + fxhe * (SAHET0 + SAHET1) * (sx2qtt - sx2qt)
            upsum = upsum + fxhe * (SAHET0 + SAHET1) * sx2qtp
        utsum = utsum + fxhe * SAHET0 * (sx1qtt - sx1qt)   # label 47
        upsum = upsum + fxhe * SAHET0 * sx1qtp
    btemp = (1.0 + e) * beta14 * betai * (24.0 * betai - 9.0)   # label 48
    qqutt = ramu * (btemp + 0.75 * (2.0 + beta14 * (8.0 * betai + 1.0)) * sqet
                    + betaut * sqett + utsum * ti)
    qqutp = ramu * (-0.25 * btemp
                    + 0.75 * ((2.0 + beta14 * (4.0 * betai + 1.0)) * sqep
                              - beta14 * betai * sqet)
                    + betaut * sqetp + upsum * ti)
    qcpt = qcpi * (qqutt + pdtq * (1.0 + qdt) - pdt * qdtt)
    qcpp = qcpi * (qqutp - pdtq * (1.0 - qdp) - pdt * qdtp)
    qat = -1.0 - qdt - qcpt + qdtt / qdt
    qap = 1.0 - qdp - qcpp + qdtp / qdt
    return (dl, d, rmu, emu, qdt, qdp, qcp, dela, qdtt, qdtp, qat, qap,
            qcpt, qcpp, ksaha)                     # label 60
