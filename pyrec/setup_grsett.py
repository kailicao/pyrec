"""Python version of setup_grsett.f (YREC 7be4133): the range, units and
diffusion coefficients for the Bahcall-Loeb settling solver (GRSETT)."""

import math

import numpy as np

from pyrec import confac_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod
from pyrec.thdiff import thdiff

# Atomic weights and charges of H, He, Fe and electrons (DATA), for the full
# Thoul et al. coefficients (MHP 8/94).
NSPEC = 4
ATOMWT = np.array([0.0, 1.008, 4.004, 55.86, 5.486e-4])
ZSPEC = np.array([0.0, 1.0, 2.0, 26.0, -1.0])

# Locals kept between calls (SAVE).
tap = np.zeros(NSPEC + 1)
tat = np.zeros(NSPEC + 1)
tac = np.zeros((NSPEC + 1, NSPEC + 1))
tcl = np.zeros((NSPEC + 1, NSPEC + 1))
xfrac = np.zeros(NSPEC + 1)
cfrac = np.zeros(NSPEC + 1)       # C(4) in the Fortran: concentrations


def setup_grsett(dt: float, hqpr: np.ndarray, hr: np.ndarray, hd: np.ndarray, hs1: np.ndarray, ht: np.ndarray,
                 lc: np.ndarray, m: int, stot: float, cod1: np.ndarray, cod2: np.ndarray, hcomp: np.ndarray,
                 hru: np.ndarray, htu: np.ndarray, ibegin: int, iend: int, lallcz: bool, qcod1x: np.ndarray,
                 qcod2x: np.ndarray, cod1z: np.ndarray, cod2z: np.ndarray, qcod1z: np.ndarray,
                 qcod2z: np.ndarray) -> tuple[float, float, int, int, bool]:
    """Find the shells where settling is computed (above a convective
    core with X > XMIN, below a convective envelope), convert the model to
    Bahcall and Loeb's units (conversion factors to confac_mod), and
    compute the helium (and, with LDIFZ, metal) diffusion coefficients and
    their X-derivatives for every shell.

    Translates: setup_grsett.f#SETUP_GRSETT (7be4133)

    Arguments: dt (s), stot (g), hqpr (d ln P/dr), hs1 (g): converted in
    place or returned; hr, hd, ht (logs); lc (convective flags); m; hcomp;
    hru, htu (filled: r, T in B&L units); ibegin, iend, lallcz (only
    set); cod1, cod2, qcod1x, qcod2x, cod1z, cod2z, qcod1z, qcod2z
    (filled, index 0 unused).

    Returns (dt, stot, ibegin, iend, lallcz).

    """

    c = constants_mod
    n = nmlpar_mod
    cf = confac_mod
    sdel = scrtch_mod.sdel

    # Solar radius (cm) and seconds per year, Bahcall and Loeb's values
    # (Fortran lines 107–108).
    crsun_bah = 6.9598e10
    csecyr_bah = 3.1558e7

    # LALLCZ: fully convective (or no shell qualifies), no settling
    # (Fortran line 110).
    lallcz = False

    # Convective core (Fortran lines 112–126).
    if lc[1]:
        for i in range(2, m + 1):
            if not lc[i]:
                break                                   # go to 20
        else:
            lallcz = True
            iounits_mod.write(iounits_mod.ishort,
                              "  FULLY CONVECTIVE MODEL - NO SETTLING")
            return dt, stot, ibegin, iend, lallcz       # go to 9999
        # Overshoot (to be added).
        ibegin = i - 1
    else:
        ibegin = 1

    # Hydrogen-exhausted core (MHP 6/90) (Fortran lines 128–137).
    for i in range(ibegin, m + 1):
        if hcomp[1, i] > n.xmin:
            break                                       # go to 25
    else:
        iounits_mod.write(iounits_mod.ishort,
                          f" X BELOW {n.xmin:9.6f} IN WHOLE MODEL-NO SETTLING")
        lallcz = True
        return dt, stot, ibegin, iend, lallcz
    ibegin = i                                          # label 25

    # Convective envelope (Fortran lines 139–148).
    if lc[m]:
        i = m - 1
        for i in range(m - 1, 1, -1):
            if not lc[i]:
                break                                   # go to 40
        else:
            i = 1           # the Fortran's DO index after the last pass
        # Overshoot (to be added).
        iend = i + 1
    else:
        iend = m

    # Helium-exhausted surface: the outer point is wherever Y > YMIN
    # (Fortran lines 151–160).
    for i in range(iend, 0, -1):
        if hcomp[2, i] > n.ymin:
            break                                       # go to 47
    else:
        iounits_mod.write(iounits_mod.ishort,
                          f" Y BELOW {n.ymin:9.6f} IN WHOLE MODEL-NO SETTLING")
        lallcz = True
        return dt, stot, ibegin, iend, lallcz
    iend = i                                            # label 47

    # Conversion factors (Fortran lines 165–169).
    cf.con_rad = 1.0 / crsun_bah
    cf.con_mass = 1.0e-2 * cf.con_rad**3
    cf.con_temp = 1.0e-7
    # Includes the factor 2.2 from ln Lambda. Question (see the page): ln
    # Lambda counted twice?
    cf.con_time = 2.7e13 * csecyr_bah

    # To Bahcall and Loeb's units (Fortran lines 172–181).
    for i in range(1, m + 1):
        hru[i] = math.exp(c.cln * hr[i]) * cf.con_rad
        htu[i] = math.exp(c.cln * ht[i]) * cf.con_temp
        hs1[i] = hs1[i] * cf.con_mass
        hqpr[i] = hqpr[i] / cf.con_rad
    dt = dt / cf.con_time
    stot = stot * cf.con_mass

    # Diffusion coefficients, for every shell (BC 5/90) (Fortran lines
    # 189–327).
    for i in range(1, m + 1):                           # loop 60
        # The Coulomb logarithm (Fortran lines 191–204). Question (see the
        # page): which ILAMBDA is recommended?
        x = hcomp[1, i]
        if n.ilambda == 2:
            # Noerdlinger (1977).
            ln_lambda = (-1.97e1 - 0.5 * math.log(1.0 / (0.75 * x + 0.25))
                         - 0.5 * hd[i] * c.cln + 1.5 * ht[i] * c.cln)
        elif n.ilambda == 3:
            # Loeb (1989).
            ln_lambda = -1.9105747e1 - 0.50 * hd[i] * c.cln \
                + 1.5 * ht[i] * c.cln
        else:
            ln_lambda = 2.2

        # Fortran lines 206–210.
        fac = n.fgry * hru[i]**2 * htu[i]**2.5 / ln_lambda
        z = hcomp[3, i] + hcomp[4, i]
        zz = hcomp[3, i]
        x2 = x * x

        if not n.lthoul:
            # Bahcall and Loeb (Fortran lines 211–222).
            xz = x * z
            x3 = x2 * x  # noqa: F841  (unused in the Fortran too)
            fac2 = 5.4 + 6.3 * x - 4.5 * x2
            # Question (see the page): +0.32 here, -0.32 in the comment.
            cod1[i] = fac * hqpr[i] * (x - x2 - xz) * (
                1.25 + sdel[2, i] * 6.0 * (x + 0.32) / fac2)
            cod2[i] = fac * (x + 3.0) / (5.0 * x2 + 8.0 * x + 3.0)
            qcod1x[i] = fac * hqpr[i] * (
                (1.0 - 2.0 * x - z)
                * (1.25 + (6.0 * sdel[2, i] * (x + 0.32)) / fac2)
                + (x - x2 - xz) * 6.0 * sdel[2, i]
                * (3.384 + 2.88 * x + 4.5 * x2) / fac2**2)
            qcod2x[i] = (-fac * (5.0 * x2 + 3.0e1 * x + 2.1e1)
                         / (5.0 * x2 + 8.0 * x + 3.0)**2)
        else:
            # Thoul et al.: A_P and A_T from THDIFF or the fits (Fortran
            # lines 224–281).
            xfrac[1] = hcomp[1, i]
            xfrac[2] = hcomp[2, i]
            xfrac[3] = hcomp[3, i]
            if not n.lthoulfit:
                if n.ilambda == 4:
                    _pair_coulomb_logs(hd[i], ht[i])
                else:
                    for ii in range(1, NSPEC + 1):
                        for jj in range(1, NSPEC + 1):
                            tcl[ii, jj] = ln_lambda
                thdiff(NSPEC, ATOMWT, ZSPEC, xfrac, tcl, tap, tat, tac)
                ap = -tap[1]
                at = -sdel[2, i] * tat[1]
            else:
                ap = 1.58 - 2.42 * x + 0.844 * x2
                at = sdel[2, i] * (1.90 - 2.69 * x + 0.805 * x2)

            # A_C and the derivatives, from the fits (Fortran lines
            # 282–285). Question (see the page): derivatives from the fits
            # always.
            ac = 1.15 - 1.42 * x + 0.647 * x2
            qapx = -2.42 + 1.688 * x
            qatx = sdel[2, i] * (-2.69 + 1.61 * x)
            qacx = -1.42 + 1.294 * x

            # Mimic mixing, to reduce settling (FD 10/09) (Fortran lines
            # 288–291). Suspected bug (see the page): the derivatives lack
            # CSTMIXING.
            cod1[i] = n.cstmixing * fac * hqpr[i] * x * (ap + at)
            cod2[i] = n.cstmixing * fac * ac
            qcod1x[i] = fac * hqpr[i] * (ap + at + x * (qapx + qatx))
            qcod2x[i] = fac * qacx

        # Metals: Thoul et al.'s coefficients for fully ionized iron
        # (Fortran lines 295–319); PARMIN sets LTHOUL with LDIFZ.
        if n.ldifz:
            fac = n.fgrz * hru[i]**2 * htu[i]**2.5 / ln_lambda
            if n.lthoul:
                if n.lthoulfit:
                    ap = -0.157 - 0.511 * x + 0.389 * x2
                    at = sdel[2, i] * (-1.36 - 1.42 * x + 0.549 * x2)
                else:
                    ap = -tap[3]
                    at = -sdel[2, i] * tat[3]
                ah = -0.0375 - 0.193 * x + 0.107 * x2
                cod1z[i] = (n.cstdiffmix * n.cstmixing * fac * hqpr[i]
                            * zz * (ap + at))
                # Positive diffusion coefficients needed.
                cod2z[i] = n.cstmixing * abs(fac * ah)
                # Suspected bug (see the page): QCOD1Z lacks CSTDIFFMIX.
                qcod1z[i] = n.cstmixing * fac * hqpr[i] * (ap + at)
                qcod2z[i] = 0.0

        # Noerdlinger's ln Lambda adds a term to d(D1)/dX (Fortran lines
        # 323–326).
        if n.ilambda == 2:
            qcod1x[i] = qcod1x[i] + cod1[i] * 1.5 / (
                ln_lambda * (3.0 * x + 1.0))

    return dt, stot, ibegin, iend, lallcz              # label 9999


def _pair_coulomb_logs(hdi: float, hti: float) -> None:
    """ILAMBDA = 4: a Coulomb logarithm for each pair of species, from the
    larger of the Debye length and the interionic distance (into tcl)
    (Fortran lines 229–267).

    """

    c = constants_mod
    rho = math.exp(c.cln * hdi)
    t = math.exp(c.cln * hti)

    # Concentrations from mass fractions (Fortran lines 232–239).
    zxa = 0.0
    for ii in range(1, NSPEC):
        zxa = zxa + ZSPEC[ii] * xfrac[ii] / ATOMWT[ii]
    for ii in range(1, NSPEC):
        cfrac[ii] = xfrac[ii] / (ATOMWT[ii] * zxa)
    cfrac[NSPEC] = 1.0

    # Electron density from the mass density (Fortran lines 241–245).
    ac = 0.0
    for ii in range(1, NSPEC + 1):
        ac = ac + ATOMWT[ii] * cfrac[ii]
    ne = rho / (1.6726e-24 * ac)

    # Interionic distance (Fortran lines 247–251).
    ni = 0.0
    for ii in range(1, NSPEC):
        ni = ni + cfrac[ii] * ne
    ao = (0.23873 / ni)**c.cc13

    # Debye length (Fortran lines 253–257).
    cz = 0.0
    for ii in range(1, NSPEC + 1):
        cz = cz + cfrac[ii] * ZSPEC[ii]**2
    lambdad = 6.9010 * math.sqrt(t / (ne * cz))

    # The Coulomb logarithms (Fortran lines 259–267).
    lam = max(lambdad, ao)
    for ii in range(1, NSPEC + 1):
        for jj in range(1, NSPEC + 1):
            xij = 2.3939e3 * t * lam / abs(ZSPEC[ii] * ZSPEC[jj])
            tcl[ii, jj] = 0.81245 * math.log(1.0 + 0.18769 * xij**1.2)
