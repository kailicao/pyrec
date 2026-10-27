"""Python version of liburn.f (YREC development, 09c3b9b)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod
from pyrec import iounits_mod
from pyrec import liov_mod
from pyrec import mdphy_mod
from pyrec import newrat_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldrat_mod
from pyrec import prevcz_mod
from pyrec import scrtch_mod
from pyrec.ratext import ratext

# DATA: the numbers of substeps tried, and the relative tolerances.
NSTEPE = [0, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96]
YTOL = [0.0, 1.0e-6, 1.0e-6, 1.0e-6]


def _converged(yest, yerr):
    """Loops 30 and 110: relative errors of the extrapolation within YTOL.

    Abundances below 1e-24 are not tested. ``yerr`` is made relative in
    place, as in the Fortran.
    """
    ldone = True
    for k2 in range(1, 4):
        if yest[k2] < 1.0e-24:
            continue
        yerr[k2] = abs(yerr[k2] / yest[k2])
        if yerr[k2] > YTOL[k2]:
            ldone = False
    return ldone


def _write_911(i, yerr):
    """Format 911: convergence failure of the extrapolation."""
    iounits_mod.write(iounits_mod.ishort,
                      f" ***LIBURN CONVERGENCE FAILURE IN SHELL {i:4d}ERRORS ")
    iounits_mod.write(iounits_mod.ishort,
                      "".join(f"{yerr[k2]:10.3E}" for k2 in range(1, 4)))


def liburn(delts, hcomp, hr, hs1, hs2, ht, jnenv, jnenv0, m):
    """Burn Li6, Li7 and Be9 over one time step (non-rotating mixing).

    The radiative interior is burned shell by shell and the surface
    convection zone as one mixed reservoir, the rates interpolated in their
    logs between the start (``oldrat_mod``) and end (``newrat_mod``) of the
    step, with substeps and rational extrapolation (RATEXT) to zero step
    size. Points left behind by a retreating zone are burned once with a
    blended rate.

    Translates: liburn.f#LIBURN (09c3b9b)

    Parameters
    ----------
    delts : float
        The time step (main.f passes seconds; see the page, question 1).
    hcomp : ndarray, shape (16, JSON + 1)
        Composition; rows 13, 14, 15 (Li6, Li7, Be9) updated in place.
    hr : ndarray
        log r at each point.
    hs1, hs2 : ndarray
        Mass inside each point, and each point's mass (g). ``hs2[jenv]``
        is changed and restored.
    ht : ndarray
        log T at each point.
    jnenv, jnenv0 : int
        The base of the surface convection zone at the end and the start
        of the step.
    m : int
        Number of points.

    Also changes ``newrat_mod`` (the bottom point's rates, question 8) and
    ``prevcz_mod`` (``rbase0``, ``bli6old``, ``bli7old``, ``bbe9old``,
    ``jenvos``). Nothing is returned.
    """
    json = dimensions_mod.json
    cln = constants_mod.cln
    rli6 = newrat_mod.rli6
    rli7 = newrat_mod.rli7
    rbe9 = newrat_mod.rbe9
    rli60 = oldrat_mod.rli60
    rli70 = oldrat_mod.rli70
    rbe90 = oldrat_mod.rbe90
    hro = oldmod_mod.hro
    sdel = scrtch_mod.sdel

    # SAVE in the Fortran: no local is read before it is set in a call.
    hsav = np.zeros((4, json + 1))
    ratli6 = np.zeros(json + 1)
    ratli7 = np.zeros(json + 1)
    ratbe9 = np.zeros(json + 1)
    yest = np.zeros(4)
    yerr = np.zeros(4)
    yext = np.zeros(4)
    fx = 0.0

    for i in range(1, m + 1):  # loop 10
        hsav[1, i] = hcomp[13, i]
        hsav[2, i] = hcomp[14, i]
        hsav[3, i] = hcomp[15, i]

    # The base of the convection zone: the offset fx (in zones) where
    # del_ad - del_rad crosses zero, between jnenv - 1 and jnenv.
    if 1 < jnenv < m:
        if nmlpar_mod.lrot and nmlpar_mod.linstb:
            ddj = mdphy_mod.delam[jnenv] - mdphy_mod.delrm[jnenv]
            ddjm1 = mdphy_mod.delam[jnenv - 1] - mdphy_mod.delrm[jnenv - 1]
        else:
            ddj = sdel[3, jnenv] - sdel[1, jnenv]
            ddjm1 = sdel[3, jnenv - 1] - sdel[1, jnenv - 1]
        fx = max(-0.5, 0.5 - ddjm1 / (ddjm1 - ddj))
        fx = min(0.5, fx)
        if not nmlpar_mod.lovste:
            jenv = jnenv
            jenv0 = jnenv0
        else:
            # Envelope overshoot: move both bases down by PSCAL0, PSCAL1.
            if prevcz_mod.rbase0 == 0.0:
                prevcz_mod.rbase0 = 0.5 * (math.exp(cln * hro[jnenv0])
                                           + math.exp(cln * hro[jnenv0 - 1]))
                ros = prevcz_mod.rbase0 - liov_mod.pscal0
                jenv0 = 1
                for j in range(jnenv0 - 1, 0, -1):
                    r = math.exp(cln * hro[j])
                    if r < ros:
                        jenv0 = j + 1
                        break  # to label 11
            else:
                jenv0 = prevcz_mod.jenvos
            dr = math.exp(cln * hr[jnenv]) - math.exp(cln * hr[jnenv - 1])
            rbase = (0.5 * (math.exp(cln * hr[jnenv])
                            + math.exp(cln * hr[jnenv - 1])) - fx * dr)
            prevcz_mod.rbase0 = rbase
            ros = rbase - liov_mod.pscal1
            jenv = 1
            for j in range(jnenv - 1, 0, -1):
                r = math.exp(cln * hr[j])
                if r < ros:
                    jenv = j + 1
                    # Question (see the page): HRO, the old radius, here.
                    dr = math.exp(cln * hro[j + 1]) - r
                    fx = 0.5 - ((ros - r) / dr)
                    fx = max(-0.5, fx)
                    fx = min(0.5, fx)
                    break  # to label 12
    else:
        jenv = jnenv
        jenv0 = jnenv0
    prevcz_mod.jenvos = jenv

    # The radiative interior, below both bases (loop 50).
    minj = min(jenv, jenv0)
    maxj = max(jenv, jenv0)
    for i in range(1, minj):
        # Question (see the page): leaves the loop, not the point.
        if rbe9[i] <= 1.0e-32 or rbe90[i] <= 1.0e-32:
            break  # to label 60
        if (hcomp[13, i] < 1.0e-24 and hcomp[14, i] < 1.0e-24
                and hcomp[15, i] < 1.0e-24):
            continue
        if ht[i] > 7.0:
            hcomp[13, i] = 0.0
            hcomp[14, i] = 0.0
            hcomp[15, i] = 0.0
            continue
        ldone = False
        for j in range(1, 12):  # loop 40
            frac = 1.0 / float(NSTEPE[j])
            dt = frac * delts
            hcomp[13, i] = hsav[1, i]
            hcomp[14, i] = hsav[2, i]
            hcomp[15, i] = hsav[3, i]
            # log rates, started half a substep before the step.
            tli6 = (math.log(rli60[i])
                    - 0.5 * frac * (math.log(rli6[i]) - math.log(rli60[i])))
            tli7 = (math.log(rli70[i])
                    - 0.5 * frac * (math.log(rli7[i]) - math.log(rli70[i])))
            tbe9 = (math.log(rbe90[i])
                    - 0.5 * frac * (math.log(rbe9[i]) - math.log(rbe90[i])))
            for k in range(1, NSTEPE[j] + 1):  # loop 20
                tli6 = tli6 + frac * (math.log(rli6[i]) - math.log(rli60[i]))
                tli7 = tli7 + frac * (math.log(rli7[i]) - math.log(rli70[i]))
                tbe9 = tbe9 + frac * (math.log(rbe9[i]) - math.log(rbe90[i]))
                ratli6[i] = dt * math.exp(tli6)
                ratli7[i] = dt * math.exp(tli7)
                ratbe9[i] = dt * math.exp(tbe9)
                if ratli6[i] < 3.0e1:
                    hcomp[13, i] = hcomp[13, i] / math.exp(ratli6[i])
                else:
                    hcomp[13, i] = 0.0
                if ratli7[i] < 3.0e1:
                    hcomp[14, i] = hcomp[14, i] / math.exp(ratli7[i])
                else:
                    hcomp[14, i] = 0.0
                if ratbe9[i] < 3.0e1:
                    hcomp[15, i] = hcomp[15, i] / math.exp(ratbe9[i])
                else:
                    hcomp[15, i] = 0.0
            yest[1] = hcomp[13, i]
            yest[2] = hcomp[14, i]
            yest[3] = hcomp[15, i]
            xest = frac ** 2
            ratext(j, xest, yest, yext, yerr, 3, 7)
            if j > 1:
                ldone = _converged(yest, yerr)
                if ldone:
                    break  # to label 45
        if not ldone:
            _write_911(i, yerr)
        # label 45
        hcomp[13, i] = yext[1]
        if hcomp[13, i] < 1.0e-24:
            hcomp[13, i] = 0.0
        hcomp[14, i] = yext[2]
        if hcomp[14, i] < 1.0e-24:
            hcomp[14, i] = 0.0
        hcomp[15, i] = yext[3]
        if hcomp[15, i] < 1.0e-24:
            hcomp[15, i] = 0.0

    # label 60: the convection zone.
    if rbe90[jenv0] <= 1.0e-32 or rbe9[jenv] <= 1.0e-32:
        return  # label 200

    # Mean abundances over the zone at the start of the step (loop 65).
    fli60 = 0.0
    fli70 = 0.0
    fbe90 = 0.0
    fm0 = 0.0
    for i in range(jenv0, m + 1):
        fli60 = fli60 + hcomp[13, i] * hs2[i]
        fli70 = fli70 + hcomp[14, i] * hs2[i]
        fbe90 = fbe90 + hcomp[15, i] * hs2[i]
        fm0 = fm0 + hs2[i]
    fli60 = fli60 / fm0
    fli70 = fli70 / fm0
    fbe90 = fbe90 / fm0
    # log of the zone's mean start rates, or the previous step's end ones.
    if prevcz_mod.bli6old <= 0.0:
        bli60 = 0.0
        bli70 = 0.0
        bbe90 = 0.0
        for i in range(jenv0, m + 1):  # loop 68
            bli60 = bli60 + rli60[i] * hs2[i]
            bli70 = bli70 + rli70[i] * hs2[i]
            bbe90 = bbe90 + rbe90[i] * hs2[i]
        bli60 = math.log(bli60 / fm0)
        bli70 = math.log(bli70 / fm0)
        bbe90 = math.log(bbe90 / fm0)
    else:
        bli60 = prevcz_mod.bli6old
        bli70 = prevcz_mod.bli7old
        bbe90 = prevcz_mod.bbe9old

    # The bottom point's mass and rates adjusted for the offset fx.
    # Question (see the page): the rates are not restored afterwards.
    hs2sav = hs2[jenv]
    if 1 < jenv < m:
        hs2[jenv] = hs2[jenv] + fx * (hs1[jenv] - hs1[jenv - 1])
        rli6[jenv] = rli6[jenv] + 0.5 * fx * (rli6[jenv - 1] - rli6[jenv])
        rli7[jenv] = rli7[jenv] + 0.5 * fx * (rli7[jenv - 1] - rli7[jenv])
        rbe9[jenv] = rbe9[jenv] + 0.5 * fx * (rbe9[jenv - 1] - rbe9[jenv])

    # Mean abundances and rates over the zone at the end (loop 70).
    fli6 = 0.0
    fli7 = 0.0
    fbe9 = 0.0
    fm = 0.0
    bli6 = 0.0
    bli7 = 0.0
    bbe9 = 0.0
    for i in range(jenv, m + 1):
        bli6 = bli6 + rli6[i] * hs2[i]
        bli7 = bli7 + rli7[i] * hs2[i]
        bbe9 = bbe9 + rbe9[i] * hs2[i]
        fli6 = fli6 + hcomp[13, i] * hs2[i]
        fli7 = fli7 + hcomp[14, i] * hs2[i]
        fbe9 = fbe9 + hcomp[15, i] * hs2[i]
        fm = fm + hs2[i]
    fli6 = fli6 / fm
    fli7 = fli7 / fm
    fbe9 = fbe9 / fm
    bli6 = math.log(bli6 / fm)
    bli7 = math.log(bli7 / fm)
    bbe9 = math.log(bbe9 / fm)
    hs2[jenv] = hs2sav

    # A deepened zone starts from the end-of-step mean.
    if jenv < jenv0:
        fli60 = fli6
        fli70 = fli7
        fbe90 = fbe9

    ldone = False
    for j in range(1, 12):  # loop 120
        frac = 1.0 / float(NSTEPE[j])
        dt = frac * delts
        fli6 = fli60
        fli7 = fli70
        fbe9 = fbe90
        czli6 = bli60 - 0.5 * frac * (bli6 - bli60)
        czli7 = bli70 - 0.5 * frac * (bli7 - bli70)
        czbe9 = bbe90 - 0.5 * frac * (bbe9 - bbe90)
        # Accretion: material of composition COMPACC added each substep.
        if nmlpar_mod.lmdot and nmlpar_mod.dmdt0 > 0.0:
            ldoacc = True
            fli6a = 0.0
            fli7a = 0.0
            fbe9a = 0.0
        else:
            ldoacc = False
        for k in range(1, NSTEPE[j] + 1):  # loop 100
            czli6 = czli6 + frac * (bli6 - bli60)
            czli7 = czli7 + frac * (bli7 - bli70)
            czbe9 = czbe9 + frac * (bbe9 - bbe90)
            dli6 = dt * math.exp(czli6)
            dli7 = dt * math.exp(czli7)
            dbe9 = dt * math.exp(czbe9)
            if dli6 < 3.0e1 and fli6 > 1.0e-24:
                fli6 = fli6 / math.exp(dli6)
            else:
                fli6 = 0.0
            if fli6 < 1.0e-24:
                fli6 = 0.0
            if dli7 < 3.0e1 and fli7 > 1.0e-24:
                fli7 = fli7 / math.exp(dli7)
            else:
                fli7 = 0.0
            if fli7 < 1.0e-24:
                fli7 = 0.0
            if dbe9 < 3.0e1 and fbe9 > 1.0e-24:
                fbe9 = fbe9 / math.exp(dbe9)
            else:
                fbe9 = 0.0
            if fbe9 < 1.0e-24:
                fbe9 = 0.0
            if ldoacc:
                # Accreted material burns half a substep on arrival.
                fli6add = nmlpar_mod.compacc[13] * frac / math.exp(0.5 * dli6)
                fli6a = fli6a / math.exp(dli6) + fli6add
                fli7add = nmlpar_mod.compacc[14] * frac / math.exp(0.5 * dli7)
                fli7a = fli7a / math.exp(dli7) + fli7add
                fbe9add = nmlpar_mod.compacc[15] * frac / math.exp(0.5 * dbe9)
                fbe9a = fbe9a / math.exp(dbe9) + fbe9add
        if ldoacc:
            fmassacc = deuter_mod.fmassacc
            print("".join(f"{v:12.3E}" for v in (fli60, fli6, fli6a, fli70,
                                                  fli7, fli7a, fmassacc)))
            fli6 = (fli6 * fm + fli6a * fmassacc) / (fm + fmassacc)
            fli7 = (fli7 * fm + fli7a * fmassacc) / (fm + fmassacc)
            fbe9 = (fbe9 * fm + fbe9a * fmassacc) / (fm + fmassacc)
        yest[1] = fli6
        yest[2] = fli7
        yest[3] = fbe9
        xest = frac ** 2
        ratext(j, xest, yest, yext, yerr, 3, 7)
        if j > 1:
            ldone = _converged(yest, yerr)
            if ldone:
                break  # to label 125
    if not ldone:
        _write_911(jenv, yerr)
    # label 125
    fli6 = yext[1]
    if fli6 < 1.0e-24:
        fli6 = 0.0
    fli7 = yext[2]
    if fli7 < 1.0e-24:
        fli7 = 0.0
    fbe9 = yext[3]
    if fbe9 < 1.0e-24:
        fbe9 = 0.0
    # Question (see the page): from maxj, so points newly in a deepened
    # zone keep their old abundances.
    for i in range(maxj, m + 1):  # loop 130
        hcomp[13, i] = fli6
        hcomp[14, i] = fli7
        hcomp[15, i] = fbe9
    prevcz_mod.bli6old = bli6
    prevcz_mod.bli7old = bli7
    prevcz_mod.bbe9old = bbe9

    # The region left by a retreating zone: burned once, blended rate.
    if jenv <= jenv0:
        return  # label 200
    if jenv0 > 1:
        hsbeg = 0.5 * (hs1[jenv0] + hs1[jenv0 - 1])
    else:
        hsbeg = 0.0
    hsend = 0.5 * (hs1[jenv] + hs1[jenv - 1])
    # Question (see the page): dmass < 0 here, so frad lies in (-1, 0).
    dmass = hsbeg - hsend
    for i in range(jenv0, jenv):  # loop 140
        if rbe9[i] <= 1.0e-32:
            break  # to label 145
        frad = (hs1[i] - hsbeg) / dmass
        dli6 = delts * math.exp(frad * math.log(rli6[i]) + (1.0 - frad) * bli60)
        dli7 = delts * math.exp(frad * math.log(rli7[i]) + (1.0 - frad) * bli70)
        dbe9 = delts * math.exp(frad * math.log(rbe9[i]) + (1.0 - frad) * bbe90)
        hcomp[13, i] = safedivexp(hcomp[13, i], dli6)
        if hcomp[13, i] < 1.0e-24:
            hcomp[13, i] = 0.0
        hcomp[14, i] = safedivexp(hcomp[14, i], dli7)
        if hcomp[14, i] < 1.0e-24:
            hcomp[14, i] = 0.0
        # Question (see the page): Be9 is not guarded by SAFEDIVEXP.
        hcomp[15, i] = hcomp[15, i] / math.exp(dbe9)
        if hcomp[15, i] < 1.0e-24:
            hcomp[15, i] = 0.0


def safedivexp(numerator, exponent):
    """KC's guard (2025-05-31): numerator / exp(exponent), or 0 on overflow.

    Translates: liburn.f#SAFEDIVEXP (09c3b9b)

    Returns ``numerator``, the Fortran's in-out argument.
    """
    if exponent < 709.7827:
        numerator = numerator / math.exp(exponent)
    else:
        numerator = 0.0
    return numerator
