"""Python version of rates.f (YREC 09c3b9b): the nuclear reaction rates at
one shell, for the abundance updates (a stripped-down ENGEB)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import cross_mod
from pyrec import nmlpar_mod

NELEM = 13
NRXNS = 13

# Thresholds still written as single-precision constants in the Fortran
# (1.E-2, 1.E-30, 1.E-5), kept at their single-precision values.
_SP_1EM2 = float(np.float32(1.0e-2))
_SP_1EM30 = float(np.float32(1.0e-30))
_SP_1EM5 = float(np.float32(1.0e-5))

# DATA (index 0 unused). Isotopes: n, H1, D, H3, He3, He4, C12, C13, N14,
# O16, O18, Ne20, Mg24. ANUC: bare nuclear masses (amu); ZNUC: charges.
ANUC = np.array([0.0, 1.008665, 1.007276, 2.013553, 3.015501, 3.014933,
                 4.001506, 11.996709, 13.000064, 13.999233, 15.990526,
                 17.994772, 19.986954, 23.978458])
ZNUC = np.array([0.0, 0., 1., 1., 1., 2., 2., 6., 6., 7., 8., 8., 10., 12.])
# The coefficients of eq. 3.14 of Neutrino Astrophysics, reactions 1-7 and
# Be7+p (8).
Q1 = np.array([0.0, 0.12317, .03392, .0325, .0304, .03035, .0273, .02494,
               .040572])
Q2 = np.array([0.0, 1.08749, -.273, -.2085, .7630, -0.4044, -1.60, -1.224,
               -0.2095])
Q3 = np.array([0.0, .93833, -.0648, -.0474, .1626, -.08598, -.3064, -.2139,
               -0.0595])
Q4 = np.array([0.0, 0., 0., 0., 4.79, 7.456, 0.0, .69703, 0.16762])
Q5 = np.array([0.0, 0., 0., 0., 2.595, 4.032, 0.0, 0.3097, 0.12114])
Q6 = np.array([0.0, -3.3804, -12.2757, -12.826, -13.6899, -13.7173,
               -15.2281, -16.6925])
Q7 = np.array([0.0, 20.8964, 76.6003, 67.8036, 69.130, 70.3809, 69.8517,
               70.8012])
Q8 = np.array([0.0, 0., 0., 0., 0.0, 0.0, 0.0, 0.])
# Screening: charge products and the coefficients of Graboske et al. 1973.
ZPRD = np.array([0.0, 1., 4., 4., 6., 6., 7., 8., 12., 16., 12., 14., 12.,
                 36.])
Z53 = np.array([0.0, 1.175, 3.73, 3.73, 4.804, 4.804, 5.385, 5.941, 9.014,
                11.24, 9.014, 10.15, 9.104, 23.28])
Z43 = np.array([0.0, 0.52, 1.31, 1.31, 1.488, 1.488, 1.61, 1.721, 2.577,
                3.025, 2.577, 2.81, 2.577, 5.668])
Z23 = np.array([0.0, -0.413, -0.655, -0.655, -0.643, -0.643, -0.659,
                -0.673, -0.889, -0.946, -0.889, -0.92, -0.889, -1.36])
Z86 = np.array([0.0, 1.630, 5.917, 5.917, 8.302, 8.302, 9.520, 10.716,
                16.192, 20.978, 16.192, 18.606, 16.192, 45.6635])
# (10^9 years / 1 second) * (1 amu / 1 gram).
C21 = 5.240358e-8

# Locals kept between calls (SAVE): the branching fractions, which reach
# the output unchanged below TCUT(1) (question 7).
f1 = 0.0
f3 = 0.0


def rates(dl, tl, x, y, xhe3, xc12, xc13, xn14, xo16, xo18, iu, hr1, hr2,
          hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12, hr13, hf1,
          hf2):
    """Nuclear reaction rates at shell iu: the 13 rates (pp chain, CNO,
    alpha captures) with weak, intermediate or strong screening, times
    C21, and the branching fractions of N15+p (hf1) and Be7 (hf2).

    Translates: rates.f#RATES (09c3b9b)

    Arguments: log10 rho (dl), log10 T (tl), the mass fractions x ... xo18,
    the shell index iu. Arrays changed in place: hr1 ... hr13, hf1, hf2
    (element iu). Nothing is returned.
    """
    global f1, f3
    c = constants_mod
    cr = cross_mod
    nml = nmlpar_mod
    xfrac = np.zeros(NELEM + 1)
    rate = np.zeros(NRXNS + 1)
    utot = np.zeros(NRXNS + 1)

    xfrac[1] = 0.0
    xfrac[2] = x
    xfrac[3] = 0.0
    xfrac[4] = 0.0
    xfrac[5] = xhe3
    xfrac[6] = y
    xfrac[7] = xc12
    xfrac[8] = xc13
    xfrac[9] = xn14
    xfrac[10] = xo16
    xfrac[11] = xo18
    xfrac[12] = 0.0
    xfrac[13] = 0.0
    # Composition sums: 1/mu_I, 1/mu_e, and the sums for screening.
    amu = 0.0
    emu = 0.0
    xtr = 0.0
    zet = 0.0
    for i in range(1, NELEM + 1):
        trm = xfrac[i] / ANUC[i]
        amu = amu + trm
        emu = emu + trm * ZNUC[i]
        xtr = xtr + trm * ZNUC[i]**1.58
        zet = zet + trm * ZNUC[i]**2
    rwe = math.exp(c.cln * dl) * emu
    dd = dl

    if tl <= nml.tcut[1]:
        for i in range(1, NRXNS + 1):
            rate[i] = 0.0
        # go to 200. Question (see the page): f1, f3 are left from the
        # previous call.
    else:
        rho = math.exp(c.cln * dd)
        t9 = math.exp(c.cln * (tl - 9.0))
        t9p13 = t9**c.cc13
        t9p23 = t9p13**2
        t9m13 = 1.0 / t9p13
        t9m23 = t9m13**2
        t9m1 = 1.0 / t9
        t9m2 = t9m1**2
        t9m12 = 1.0 / math.sqrt(t9)
        t9m32 = t9m1 * t9m12
        # Electron degeneracy, for Salpeter's screening.
        pfmc2 = 1.017677e-4 * rwe**c.cc23
        efmkt = 5.92986 * t9m1 * (math.sqrt(1.0 + pfmc2) - 1.0)
        if efmkt <= _SP_1EM2:
            fprf = 1.0
        else:
            degd = math.log10(efmkt)
            if degd >= 1.5:
                fprf = 0.0
            else:
                fprf = (0.75793 - 0.54621 * degd - 0.30964 * degd**2
                        + 0.12535 * degd**3 + 0.1203 * degd**4
                        - 0.012857 * degd**5 - 0.014768 * degd**6)
        # Screening.
        xxl = 5.9426e-6 * t9m32 * math.sqrt(rho * amu)
        xxl6 = xxl**c.cc23
        xxl8 = xxl**0.86
        zcurl = math.sqrt((zet + fprf * emu) / amu)
        zbar = emu / amu
        z58 = zcurl**0.58
        z28 = zbar**0.28
        z33 = zbar**c.cc13
        tm1 = xxl * zcurl
        for i in range(1, NRXNS + 1):
            uwk = tm1 * ZPRD[i]
            if uwk <= nml.weakscreening:
                utot[i] = uwk
            else:
                uint = 0.38 * xxl8 * xtr * Z86[i] / (amu * z58 * z28)
                if uwk <= 2.0:
                    utot[i] = uint
                else:
                    ustr = 0.624 * z33 * xxl6 * (
                        Z53[i] + 0.316 * z33 * Z43[i]
                        + 0.737 * Z23[i] / (zbar * xxl6))
                    if ustr < uint or uwk >= 5.0:
                        utot[i] = ustr
                    else:
                        utot[i] = uint

        nz = 1
        if x == 0.0:
            f1 = 0.0
            f3 = 0.0
        else:
            nz = 8
            # Reactions 1-7.
            for i in range(1, 8):
                r1 = (t9m23 + Q1[i] * t9m13
                      + cr.qs0e[i] * (Q2[i] + Q3[i] * t9p13)
                      + cr.qqs0ee[i] * (Q4[i] * t9p23 + Q5[i] * t9))
                rate[i] = rho * r1 * math.exp(Q6[i] * t9m13 + Q7[i]
                                              + (Q8[i] * t9)**2 + utot[i])
                rate[i] = rate[i] * nml.sstandard[i]
                if rate[i] < _SP_1EM30:
                    rate[i] = 0.0
            # Be7: electron and proton capture.
            be7electron = 1.752e-10 * t9m12 * (1.0 + 0.004 * (1000.0 * t9 - 16.0))
            be7electron = be7electron * emu * nml.sstandard[15]
            temp3 = -10.2625 * t9m13
            be7proton = 3.128813e+5 * x * nml.sstandard[16] * math.exp(temp3)
            qrbe7 = (t9m23 + Q1[8] * t9m13
                     + cr.qs0e[8] * (Q2[8] + Q3[8] * t9p13)
                     + cr.qqs0ee[8] * (Q4[8] * t9p23 + Q5[8] * t9))
            be7proton = be7proton * qrbe7
            # Question (see the page): Be7+p's own screening constants,
            # with no strong branch.
            zprdbe7p = 4.0
            z86be7p = 5.7790
            uwk = tm1 * zprdbe7p
            if uwk <= nml.weakscreening:
                utotbe7p = uwk
            else:
                uint = 0.38 * xxl8 * xtr * z86be7p / (amu * z58 * z28)
                utotbe7p = uint
            be7proton = be7proton * math.exp(utotbe7p)
            camube7 = rho * 8.584981e+22
            be7proton = camube7 * be7proton
            be7electron = camube7 * be7electron
            f1 = be7electron / (be7electron + be7proton)
            # N15+p: the (p,gamma) and (p,alpha) branches.
            o16gamma = (t9m23 + 0.0273016 * t9m13 + 0.14374
                        + 0.027490 * t9p13 + 6.14685 * t9p23 + 2.98940 * t9)
            o16gamma = o16gamma * 64 * cr.fo16
            # Question (see the page): the comment says 675000.
            c12alpha = (t9m23 + 0.0273016 * t9m13 + 2.01186
                        + 0.384763 * t9p13 + 17.0579 * t9p23 + 8.29580 * t9)
            c12alpha = c12alpha * 67500 * cr.fc12
            f3 = c12alpha / (c12alpha + o16gamma)
            f4 = 1.0 - f3
        # label 50
        for i in range(nz, NRXNS + 1):
            rate[i] = 0.0
        if tl >= nml.tcut[4]:
            # C13(alpha,n)O16
            r1 = t9m23 + 0.0129 * t9m13 + 2.04 + 0.184 * t9p13
            a1 = 6.77e15 * math.exp(-32.329 * t9m13 - (t9 / 1.284)**2)
            a2 = 3.82e5 * math.exp(-9.373 * t9m1)
            a3 = 1.41e6 * math.exp(-11.873 * t9m1)
            a4 = 2.00e9 * math.exp(-20.409 * t9m1)
            a5 = 2.92e9 * math.exp(-29.283 * t9m1)
            rate[8] = 1.157126e22 * rho * math.exp(utot[8]) * (
                a1 * r1 + t9m32 * (a2 + a3 + a4 + a5))
            # C12(alpha,gamma)O16
            r1 = 1.0 / (1.0 + 0.0489 * t9m23)
            r2 = 1.0 / (1.0 + 0.2654 * t9m23)
            a1 = t9m2 * math.exp(-32.120 * t9m13)
            a2 = 1.04e8 * r1**2 * math.exp(-(t9 / 3.496)**2)
            a3 = 1.76e8 * r2**2
            a4 = 1.25e3 * t9m32 * math.exp(-27.499 * t9m1)
            a5 = 1.43e-2 * t9**5 * math.exp(-15.541 * t9m1)
            rate[10] = 1.25388e22 * rho * math.exp(utot[10]) * (
                a1 * (a2 + a3) + a4 + a5)
            # N14(alpha,gamma)F18, F18 => O18 + e+ + nu
            r1 = (t9m23 + 0.012 * t9m13 + 1.45 + 0.177 * t9p13
                  + 1.97 * t9p23 + 0.406 * t9)
            a1 = 7.78e9 * math.exp(-36.031 * t9m13 - (t9 / 0.881)**2)
            a2 = t9m32 * 2.36e-10 * math.exp(-2.798 * t9m1)
            a3 = t9m32 * 2.03 * math.exp(-5.054 * t9m1)
            a4 = t9m23 * 1.15e4 * math.exp(-12.310 * t9m1)
            rate[11] = 1.07452e22 * rho * math.exp(utot[11]) * (
                a1 * r1 + a2 + a3 + a4)
            # Triple alpha
            rate[12] = (1.565315e21 * rho**2 * t9m1 * t9m2 * 2.79e-8
                        * math.exp(-4.4027 * t9m1 + utot[12]))
        # label 100
        rate[9] = 0.0
        rate[13] = 0.0
        # Question (see the page): a second, much higher cut to zero.
        for i in range(1, NRXNS + 1):
            if rate[i] <= _SP_1EM5:
                rate[i] = 0.0

    # label 200: the rates in output units.
    hr1[iu] = rate[1] * C21
    hr2[iu] = rate[2] * C21
    hr3[iu] = rate[3] * C21
    hr4[iu] = rate[4] * C21
    hr5[iu] = rate[5] * C21
    hr6[iu] = rate[6] * C21
    hr7[iu] = rate[7] * C21
    hr8[iu] = rate[8] * C21
    hr9[iu] = rate[9] * C21
    hr10[iu] = rate[10] * C21
    hr11[iu] = rate[11] * C21
    hr12[iu] = rate[12] * C21
    hr13[iu] = rate[13] * C21
    hf1[iu] = f3
    hf2[iu] = f1
