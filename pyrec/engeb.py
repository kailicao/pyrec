"""Python version of engeb.f (YREC 09c3b9b): nuclear energy generation at
one point, its derivatives, the reaction rates, the solar neutrino fluxes
and the neutrino losses."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import cross_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod  # JSON, the array sizes (USEd by the Fortran)
from pyrec import fluxes_mod
from pyrec import grab_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.neutrino import neutrino

# ---- DATA (all arrays with index 0 unused) ----

# Old neutrino-loss fits (pair, photo, plasma): DATA V1,V2,V3/21 values/.
V1 = np.array([0.0, 6.002e19, 2.084e20, 1.872e21, 9.383e-1, -4.141e-1,
               5.829e-2, 5.5924])
V2 = np.array([0.0, 4.886e10, 7.580e10, 6.023e10, 6.290e-3, 7.483e-3,
               3.061e-4, 1.5654])
V3 = np.array([0.0, 2.320e-7, 8.449e-8, 1.787e-8, 2.581e-2, 1.734e-2,
               6.990e-4, 0.56457])
# Bare nuclear masses (amu) and charges of the 13 isotopes:
# n, p, 2H, 3H, 3He, 4He, 12C, 13C, 14N, 16O, 18O, 20Ne, 24Mg.
ANUC = np.array([0.0, 1.008665, 1.007276, 2.013553, 3.015501, 3.014933,
                 4.001506, 11.996709, 13.000064, 13.999233, 15.990526,
                 17.994772, 19.986954, 23.978458])
ZNUC = np.array([0.0, 0., 1., 1., 1., 2., 2., 6., 6., 7., 8., 8., 10., 12.])
NELEM = 13
# Rate coefficients of reactions 1-7 (eq. 3.14 form); Q1-Q5 element 8 is
# Be7 + p.
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
Q8 = np.zeros(8)
NRXNS = 13
# Screening: charge products and the strong/intermediate coefficients.
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
# (s per 1e9 yr) * (amu in g), per the comment.
C21 = 5.240358e-8

# ---- SAVE: locals kept between calls ----
# Below TCUT(1) the routine jumps to label 200, and these keep the previous
# call's values there (see the page, question 6).
xfrac = np.zeros(NELEM + 1)
rate = np.zeros(NRXNS + 1)
dratro = np.zeros(NRXNS + 1)
dratt = np.zeros(NRXNS + 1)
utot = np.zeros(NRXNS + 1)
dscr = np.zeros(NRXNS + 1)
dsct = np.zeros(NRXNS + 1)
dg = np.zeros(NRXNS + 1)
eg = np.zeros(51)
f1 = f2 = f3 = f4 = 0.0
rho = t9 = t9m1 = t9m12 = t9m13 = t9m23 = t9p13 = 0.0
tm1 = xxl8 = z58 = z28 = 0.0
be7proton = be7electron = 0.0
convert = 0.0  # set above TCUT(1); read at label 200 on every path


def _screen_weak_or_int(uwk, z86, amu, xtr):
    """Weak screening, or Graboske et al.'s intermediate form (as used for
    Be7 + p and He3 + p): UWK if at most WEAKSCREENING, else UINT."""
    if uwk <= nmlpar_mod.weakscreening:
        return uwk
    return 0.38 * xxl8 * xtr * z86 / (amu * z58 * z28)


def engeb(epp1, epp2, epp3, ecn, e3al, pep, pet, sum1, dl,
          tl, x, y, xhe3, xc12, xc13, xn14, xo16,
          xo18, xh2, iu, hr1, hr2, hr3, hr4, hr5, hr6, hr7,
          hr8, hr9, hr10, hr11, hr12, hr13, hf1, hf2, be7massfraction,
          ealpca, enu):
    """Nuclear energy generation at one point: screening, the pp/CNO rates
    with S-factor scaling, the Be7 and N15 branchings, alpha captures,
    deuterium burning, energies net of neutrinos, the derivatives PEP and
    PET, the per-reaction rates HR1-HR13, the solar neutrino fluxes and
    the neutrino losses.

    Translates: engeb.f#ENGEB (09c3b9b)

    Arguments: dl, tl (log rho, log T); the mass fractions x, y, xhe3,
    xc12, xc13, xn14, xo16, xo18, xh2; iu, the point's index; hr1-hr13,
    hf1, hf2, arrays set at index iu (index 0 unused). The other scalars
    are only set (passed in with any value).

    Returns (epp1, epp2, epp3, ecn, e3al, pep, pet, sum1, be7massfraction,
    ealpca, enu): the energy of the pp chain's three branches, CNO and
    triple alpha; d eps/d ln rho and d eps/d ln T; the total (erg/g/s);
    the Be7 mass fraction (set only with LSNU and X > 1e-6); the energy
    from alpha captures; the neutrino losses (negative). Also sets
    fluxes_mod.flux, grab_mod.enghe3 and grab_mod.he3all.
    """
    global f1, f2, f3, f4, rho, t9, t9m1, t9m12, t9m13, t9m23, t9p13
    global tm1, xxl8, z58, z28, be7proton, be7electron, convert
    c = constants_mod
    nml = nmlpar_mod
    cr = cross_mod

    enu = 0.0
    ealpca = 0.0
    # Mass fractions of the 13 isotopes (n, 3H, Ne, Mg taken as zero).
    xfrac[1] = 0.0
    xfrac[2] = x
    xfrac[3] = xh2
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
    # 1/mu_I, 1/mu_e, the sum for intermediate screening, and zeta.
    amu = 0.0
    emu = 0.0
    xtr = 0.0
    zet = 0.0
    for i in range(1, NELEM + 1):  # loop 10
        trm = xfrac[i] / ANUC[i]
        amu = amu + trm
        emu = emu + trm * ZNUC[i]
        xtr = xtr + trm * ZNUC[i]**1.58
        zet = zet + trm * ZNUC[i]**2
    rwe = math.exp(c.cln * dl) * emu
    dd = dl

    low_t = tl <= nml.tcut[1]
    if low_t:
        # Below TCUT(1): no burning; jump to label 200.
        # Question (see the page): SUM1, F1, F3 and the Be7 terms are stale.
        dgdeut = 0.0
        qrtdeut = 0.0
        en = -20.0
        pep = 0.0
        pet = 0.0
        for i in range(1, NRXNS + 1):  # loop 20
            eg[i] = 0.0
            rate[i] = 0.0
            dg[i] = 0.0
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

        # Deuterium burning (MHP 5/02), as in DEUTRATE.
        if xh2 <= 1.0e-11:
            rdeut = 0.0
            qrtdeut = 0.0
        else:
            qdeut = 5.494
            zz = -3.72 * t9m13
            tfacdeut = 1.0 + 0.112 * t9p13 + 3.38 * t9p23 + 2.65 * t9
            # Question (see the page): T9**(+2/3); the derivative below is
            # that of the T9**(-2/3) form.
            rdeut = (rho * 2.240e3 * t9p23 * math.exp(zz) * tfacdeut
                     * 6.023e23 / ANUC[3])
            tfacdeut2 = 0.112 * t9p13 + 6.76 * t9p23 + 7.95 * t9
            qrtdeut = c.cc13 * ((tfacdeut2 / tfacdeut) - 2.0 - zz)
            # Cap in the surface convection zone by the overturn time.
            # Question (see the page): IU is never set by main.f's call.
            if iu >= deuter_mod.jcz and ovrtrn_mod.taucz > 1.0:
                rdeutmax = 6.023e23 / ANUC[3] / ovrtrn_mod.taucz
                rdeut2 = rdeut * x
                if rdeut2 > rdeutmax:
                    if x > 1.0e-6:
                        rdeut = rdeutmax / x

        # Degeneracy factor f'/f, a fit to Salpeter's figure 1.
        pfmc2 = 1.017677e-4 * rwe**c.cc23
        efmkt = 5.92986 * t9m1 * (math.sqrt(1.0 + pfmc2) - 1.0)
        if efmkt <= 1.0e-2:
            fprf = 1.0
        else:
            degd = math.log10(efmkt)
            if degd >= 1.5:
                fprf = 0.0
            else:
                fprf = (0.75793 - 0.54621 * degd - 0.30964 * degd**2
                        + 0.12535 * degd**3 + 0.1203 * degd**4
                        - 0.012857 * degd**5 - 0.014768 * degd**6)

        # Screening: weak, intermediate or strong, for each reaction.
        xxl = 5.9426e-6 * t9m32 * math.sqrt(rho * amu)
        xxl6 = xxl**c.cc23
        xxl8 = xxl**0.86
        zcurl = math.sqrt((zet + fprf * emu) / amu)
        zbar = emu / amu
        z58 = zcurl**0.58
        z28 = zbar**0.28
        z33 = zbar**c.cc13
        tm1 = xxl * zcurl
        for i in range(1, NRXNS + 1):  # loop 30
            uwk = tm1 * ZPRD[i]
            if uwk <= nml.weakscreening:
                utot[i] = uwk
                dscr[i] = 0.5 * uwk
                dsct[i] = -1.5 * uwk
            else:
                uint = 0.38 * xxl8 * xtr * Z86[i] / (amu * z58 * z28)
                if uwk <= 2.0:
                    utot[i] = uint
                    dscr[i] = 0.43 * uint
                    dsct[i] = -1.29 * uint
                else:
                    ustr = 0.624 * z33 * xxl6 * (
                        Z53[i] + 0.316 * z33 * Z43[i]
                        + 0.737 * Z23[i] / (zbar * xxl6))
                    if ustr < uint or uwk >= 5.0:
                        utot[i] = ustr
                        dscr[i] = (0.208 * z33 * (Z53[i] + 0.316 * z33
                                                  * Z43[i]) * xxl6)
                        dsct[i] = -3.0 * dscr[i]
                    else:
                        utot[i] = uint
                        dscr[i] = 0.43 * uint
                        dsct[i] = -1.29 * uint

        # pp chain and CNO, reactions 1-7 (none without hydrogen).
        if x == 0.0:
            nz = 1
            f1 = 0.0
            f2 = 0.0
            f3 = 0.0
            f4 = 0.0
        else:
            nz = 8
            for i in range(1, 8):  # loop 40
                r1 = (t9m23 + Q1[i] * t9m13
                      + cr.qs0e[i] * (Q2[i] + Q3[i] * t9p13)
                      + cr.qqs0ee[i] * (Q4[i] * t9p23 + Q5[i] * t9))
                rate[i] = rho * r1 * math.exp(Q6[i] * t9m13 + Q7[i]
                                              + (Q8[i] * t9)**2 + utot[i])
                rate[i] = rate[i] * nml.sstandard[i]
                if rate[i] < 1.0e-30:
                    rate[i] = 0.0
                    dratt[i] = 0.0
                else:
                    dratro[i] = 1.0 + dscr[i]
                    # Question (see the page): no QS0E, QQS0EE factors here.
                    dratt[i] = (dsct[i] - (Q6[i] * t9m13
                                           + (2.0 * t9m23 + Q1[i] * t9m13
                                              - Q3[i] * t9p13
                                              - 2.0 * Q4[i] * t9p23
                                              - 3.0 * Q5[i] * t9) / r1) / 3.0
                                + 2.0 * (Q8[i] * t9)**2)

            # Be7: electron capture against proton capture.
            be7electron = (1.752e-10 * t9m12
                           * (1.0 + 0.004 * (1000.0 * t9 - 16.0)))
            be7electron = be7electron * emu * nml.sstandard[15]
            temp3 = -10.2625 * t9m13
            be7proton = 3.128813e+5 * x * nml.sstandard[16] * math.exp(temp3)
            qrbe7 = (t9m23 + Q1[8] * t9m13
                     + cr.qs0e[8] * (Q2[8] + Q3[8] * t9p13)
                     + cr.qqs0ee[8] * (Q4[8] * t9p23 + Q5[8] * t9))
            be7proton = be7proton * qrbe7
            # Screening of Be7 + p (ZPRD = 4).
            utotbe7p = _screen_weak_or_int(tm1 * 4.0, 5.7790, amu, xtr)
            be7proton = be7proton * math.exp(utotbe7p)
            camube7 = rho * 8.584981e+22
            be7proton = camube7 * be7proton
            be7electron = camube7 * be7electron
            f1 = be7electron / (be7electron + be7proton)
            f2 = be7proton / (be7electron + be7proton)

            # N15 + p: (p,alpha) to C12 against (p,gamma) to O16.
            o16gamma = (t9m23 + 0.0273016 * t9m13 + 0.14374
                        + 0.027490 * t9p13 + 6.14685 * t9p23 + 2.98940 * t9)
            o16gamma = o16gamma * 64 * cr.fo16
            c12alpha = (t9m23 + 0.0273016 * t9m13 + 2.01186
                        + 0.384763 * t9p13 + 17.0579 * t9p23 + 8.29580 * t9)
            c12alpha = c12alpha * 67500 * cr.fc12
            f3 = c12alpha / (c12alpha + o16gamma)
            f4 = 1.0 - f3

        for i in range(nz, NRXNS + 1):  # label 50, loop 60
            rate[i] = 0.0
            dratro[i] = 0.0
            dratt[i] = 0.0

        # Alpha captures, above TCUT(4).
        # Question (see the page): the derivatives DRATT(8, 10, 11).
        if tl >= nml.tcut[4]:
            # C13(alpha,n)O16.
            r1 = t9m23 + 0.0129 * t9m13 + 2.04 + 0.184 * t9p13
            a1 = 6.77e15 * math.exp(-32.329 * t9m13 - (t9 / 1.284)**2)
            a2 = 3.82e5 * math.exp(-9.373 * t9m1)
            a3 = 1.41e6 * math.exp(-11.873 * t9m1)
            a4 = 2.00e9 * math.exp(-20.409 * t9m1)
            a5 = 2.92e9 * math.exp(-29.283 * t9m1)
            rate[8] = (1.157126e22 * rho * math.exp(utot[8])
                       * (a1 * r1 + t9m32 * (a2 + a3 + a4 + a5)))
            dratro[8] = 1.0 + dscr[8]
            dr1 = c.cc13 * (-2.0 * t9m23 - 0.0129 * t9m13 + 0.184 * t9p13)
            da1 = a1 * (c.cc13 * 32.329 * t9m13 - 2.0 * (t9 / 1.284)**2)
            dratt[8] = dsct[8] + rho / rate[8] * (
                dr1 * a1 + r1 * da1
                + a2 * (9.373 * t9m1 - 1.5) + a3 * (11.873 * t9m1 - 1.5)
                + a4 * (20.409 * t9m1 - 1.5) + a5 * (29.283 * t9m1 - 1.5))
            # C12(alpha,gamma)O16.
            r1 = 1.0 / (1.0 + 0.0489 * t9m23)
            r2 = 1.0 / (1.0 + 0.2654 * t9m23)
            a1 = t9m2 * math.exp(-32.120 * t9m13)
            a2 = 1.04e8 * r1**2 * math.exp(-(t9 / 3.496)**2)
            a3 = 1.76e8 * r2**2
            a4 = 1.25e3 * t9m32 * math.exp(-27.499 * t9m1)
            a5 = 1.43e-2 * t9**5 * math.exp(-15.541 * t9m1)
            rate[10] = (1.25388e22 * rho * math.exp(utot[10])
                        * (a1 * (a2 + a3) + a4 + a5))
            dratro[10] = 1.0 + dscr[10]
            dratt[10] = dsct[10] + rho / rate[10] * (
                a1 * ((c.cc13 * 32.120 * t9m13 - 2.0) * (a2 + a3)
                      + a2 * (r1 * c.cc13 * 0.1956 - 2.0 * (t9 / 3.496)**2)
                      + a3 * (c.cc13 * 1.0616 * r2))
                + a4 * (27.499 * t9m1 - 1.5) + a5 * (5.0 + 15.541 * t9m1))
            # N14(alpha,gamma)F18, F18 -> O18 + e+ + nu.
            r1 = (t9m23 + 0.012 * t9m13 + 1.45 + 0.177 * t9p13
                  + 1.97 * t9p23 + 0.406 * t9)
            a1 = 7.78e9 * math.exp(-36.031 * t9m13 - (t9 / 0.881)**2)
            a2 = t9m32 * 2.36e-10 * math.exp(-2.798 * t9m1)
            a3 = t9m32 * 2.03 * math.exp(-5.054 * t9m1)
            a4 = t9m23 * 1.15e4 * math.exp(-12.310 * t9m1)
            rate[11] = (1.07452e22 * rho * math.exp(utot[11])
                        * (a1 * r1 + a2 + a3 + a4))
            dratro[11] = 1.0 + dscr[11]
            dr1 = (c.cc13 * (-2.0 * t9m23 - 0.012 * t9m13 + 0.177 * t9p13
                             + 3.94 * t9p23) + 0.406 * t9)
            da1 = a1 * (c.cc13 * 36.031 * t9m13 - 2.0 * (t9 / 0.881)**2)
            dratt[11] = dsct[11] + rho / rate[11] * (
                dr1 * a1 + r1 * da1 + a2 * (2.798 * t9m1 - 1.5)
                + a3 * (5.054 * t9m1 - 1.5) + a4 * (12.310 * t9m1 - c.cc23))
            # Triple alpha.
            rate[12] = (1.565315e21 * rho**2 * t9m1 * t9m2 * 2.79e-8
                        * math.exp(-4.4027 * t9m1 + utot[12]))
            dratro[12] = 2.0 + dscr[12]
            dratt[12] = -3.0 + dsct[12] + 4.4027 * t9m1

        # Label 100: rates times abundances.
        eg[1] = rate[1] * x * x
        if xh2 > 1.0e-11:
            egdeut = rdeut * x * xh2
        else:
            egdeut = 0.0
        eg[2] = rate[2] * xhe3 * xhe3
        eg[3] = rate[3] * xhe3 * y
        eg[4] = rate[4] * x * xc12
        eg[5] = rate[5] * x * xc13
        eg[6] = rate[6] * x * xn14
        eg[7] = rate[7] * x * xo16
        eg[8] = rate[8] * y * xc13
        eg[10] = rate[10] * y * xc12
        eg[11] = rate[11] * y * xn14
        eg[12] = rate[12] * y**3

        # Energies, MeV net of neutrinos, converted to erg.
        convert = 1.602177e-6
        dg[1] = eg[1] * 6.664 * convert
        if xh2 > 1.0e-11:
            dgdeut = egdeut * qdeut * convert
        else:
            dgdeut = 0.0
        dg[2] = eg[2] * 12.860 * convert
        dg[3] = eg[3] * (1.586 + f1 * 17.394 + f2 * 11.499) * convert
        dg[4] = eg[4] * 3.457372 * convert
        dg[5] = eg[5] * 7.550628 * convert
        dg[6] = eg[6] * (9.054 + f3 * 4.966 + f4 * 12.128) * convert
        dg[7] = eg[7] * 3.553 * convert
        dg[8] = eg[8] * 2.216 * convert
        dg[10] = eg[10] * 7.162 * convert
        dg[11] = eg[11] * 5.815 * convert
        dg[12] = eg[12] * 7.275 * convert
        grab_mod.enghe3 = dg[2]
        grab_mod.he3all = dg[2] + dg[3]
        # Reactions 9 and 13 are zeroed.
        dratro[9] = 0.0
        dratt[9] = 0.0
        dg[9] = 0.0
        dratro[13] = 0.0
        dratt[13] = 0.0
        dg[13] = 0.0

        # Total, and its absolute derivatives d eps/d ln rho, d eps/d ln T.
        # Question (see the page): the comment says relative derivatives.
        sum1 = 0.0
        sum2 = 0.0
        sum3 = 0.0
        for i in range(1, NRXNS + 1):  # loop 110
            sum1 = sum1 + dg[i]
            sum2 = sum2 + dg[i] * dratro[i]
            sum3 = sum3 + dg[i] * dratt[i]
        sum1 = sum1 + dgdeut
        sum2 = sum2 + dgdeut
        sum3 = sum3 + dgdeut * qrtdeut
        # Question (see the page): the hard-coded thresholds.
        if sum1 <= 1.0e-12:
            en = -20.0
            pep = 0.0
            pet = 0.0
            for i in range(1, NRXNS + 1):  # loop 120
                eg[i] = 0.0
        else:
            pep = sum2
            pet = sum3
        for i in range(1, NRXNS + 1):  # loop 130
            if rate[i] <= 1.0e-5:
                rate[i] = 0.0

    # Label 200: rates per Gyr per amu, for the abundance routines.
    # Question (see the page): main.f passes an IU it never sets.
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

    # Energy by branch: pp I (with deuterium), pp II, pp III; CNO.
    epp1 = dg[1] + dg[2] + dgdeut
    epp3 = eg[3] * f2 * (1.586 + 11.499) * convert
    epp2 = dg[3] - epp3
    ecn = dg[4] + dg[5] + dg[6] + dg[7]
    e3al = dg[12]
    ealpca = dg[8] + dg[10] + dg[11]

    # Solar neutrino fluxes per gram at 1 AU, in 1e10 cm^-2 s^-1.
    flux = fluxes_mod.flux
    if nml.lsnu:
        if x <= 1.0e-6:
            for i in range(1, 11):
                flux[i] = 0.0
        else:
            fourpiau2 = 2.812295e+27
            flux[1] = eg[1] / fourpiau2                      # pp
            flux[2] = (3.4848e-6 * rwe * t9m12 * (1.0 + 20.0 * t9)
                       * eg[1])                               # pep
            flux[2] = flux[2] * nml.sstandard[14] / fourpiau2
            q6hep = -6.1399
            flux[3] = 1.71724e+11 * rho * t9m23 * math.exp(q6hep * t9m13)
            flux[3] = (1.0 + 0.067862 * t9p13) * nml.sstandard[17] * flux[3]
            # Screening of He3 + p (ZPRD = 2), for the hep flux.
            utothe3p = _screen_weak_or_int(tm1 * 2.0, 3.08687, amu, xtr)
            flux[3] = flux[3] * math.exp(utothe3p)
            flux[3] = flux[3] * x * xhe3 / fourpiau2          # hep
            # Kept by main.f's SAVE between the calls that set it.
            be7massfraction = eg[3] / (be7proton + be7electron)
            flux[4] = eg[3] * f1 / fourpiau2                  # Be7
            flux[5] = eg[3] * f2 / fourpiau2                  # B8
            flux[6] = eg[4] / fourpiau2                       # N13
            flux[7] = eg[6] / fourpiau2                       # O15
            flux[8] = eg[7] / fourpiau2                       # F17
            flux[9] = eg[2] / fourpiau2                       # "fictional"
            flux[10] = eg[3] / fourpiau2                      # "fictional"
            for k in range(1, 11):
                flux[k] = 1.0e-10 * flux[k]
                temp = flux[k]
                if temp <= 0.0:
                    flux[k] = 0.0

    # Neutrino losses, above TCUT(5).
    if tl <= nml.tcut[5]:
        return (epp1, epp2, epp3, ecn, e3al, pep, pet, sum1,
                be7massfraction, ealpca, enu)
    zc = xc12 + xc13
    zo = xo16 + xo18
    t = 10.0**tl
    d = 10.0**dl
    if nml.lnulos1:
        # Itoh et al. 1996 (Grant Newsham 9/06).
        snu, dsnudt, dsnudd = neutrino(t, d, x, y, zc, zo, 0.0, 0.0, 0.0)
        enu = -snu
        # Question (see the page): relative derivatives added to absolute.
        dsnudt = -dsnudt * t / enu
        dsnudd = -dsnudd * d / enu
        sum1 = sum1 + enu
        pep = pep + dsnudd
        pet = pep + dsnudt  # Question (see the page): PEP, not PET.
    else:
        # The old fits: pair, photo- and plasma neutrinos.
        el = t9 / 5.9302
        eli = 1.0 / el
        ez = math.exp(c.cc13 * c.cln * (dl - 9.0)) * eli * (0.7937 + 0.2063 * x)
        ez3 = ez**3
        emue = 0.5 * (1.0 + x)
        ex1 = 0.0
        if t9 >= 0.2:
            el2 = el * el
            polx10 = 1.0 + el2 * (-13.04 + el2 * (133.5 + el2 * (
                1534.0 + el2 * 918.6)))
            polx11 = V1[1] + ez * (V1[2] + ez * V1[3])
            polx12 = ez3 + eli * (V1[4] + eli * (V1[5] + eli * V1[6]))
            ex1 = (math.exp(-ez * V1[7] - eli - eli - c.cln * dl)
                   * polx10 * polx11 / polx12)
        polx21 = V2[1] + ez * (V2[2] + ez * V2[3])
        polx22 = ez3 + eli * (V2[4] + eli * (V2[5] + eli * V2[6]))
        ex2 = emue * el**5 * math.exp(-ez * V2[7]) * polx21 / polx22
        polx31 = V3[1] + ez * (V3[2] + ez * V3[3])
        polx32 = ez3 + eli * (V3[4] + eli * (V3[5] + eli * V3[6]))
        ex3 = (emue**3 * math.exp(-ez * V3[7] + c.cln * (dl + dl))
               * polx31 / polx32)
        enu = -(ex1 + ex2 + ex3)
        sum1 = sum1 + enu
        # Derivatives of the losses.
        qetnx = 0.0
        qednx = 0.0
        if t9 >= 0.2:
            qedn = (ez * (V1[2] + 2.0 * ez * V1[3]) / polx11 - ez * V1[7]
                    - 3.0 * ez3 / polx12)
            qetn = -qedn + eli * (V1[4] + eli * (2.0 * V1[5]
                                                 + 3.0 * eli * V1[6])) / polx12
            qetnx = el2 * (-26.08 + el2 * (534.0 + el2 * (
                9204.0 + el2 * 7348.8))) / polx10
            qetnx = (qetnx + qetn + eli + eli) * ex1
            qednx = (-1.0 + c.cc13 * qedn) * ex1
        qedn = (ez * (V2[2] + 2.0 * ez * V2[3]) / polx21 - ez * V2[7]
                - 3.0 * ez3 / polx22)
        qetn = -qedn + eli * (V2[4] + eli * (2.0 * V2[5]
                                             + 3.0 * eli * V2[6])) / polx22
        qetnx = qetnx + (5.0 + qetn) * ex2
        qednx = qednx + c.cc13 * qedn * ex2
        qedn = (ez * (V3[2] + 2.0 * ez * V3[3]) / polx31 - ez * V3[7]
                - 3.0 * ez3 / polx32)
        qetn = -qedn + eli * (V3[4] + eli * (2.0 * V3[5]
                                             + 3.0 * eli * V3[6])) / polx32
        qetnx = qetnx + qetn * ex3
        qednx = qednx + (2.0 + c.cc13 * qedn) * ex3
        pet = pet - qetnx
        pep = pep - qednx
    return (epp1, epp2, epp3, ecn, e3al, pep, pet, sum1, be7massfraction,
            ealpca, enu)
