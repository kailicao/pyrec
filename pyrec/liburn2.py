"""Python version of liburn2.f (YREC development, 09c3b9b)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import liov_mod
from pyrec import mdphy_mod
from pyrec import newrat_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldrat_mod
from pyrec import prevcz_mod
from pyrec import scrtch_mod


def liburn2(delts, hcomp, hr, hs1, hs2, ht, jnenv, jnenv0, m):
    """Burn Li6, Li7 and Be9 over one step in one go (cheap LIBURN).

    The same as LIBURN without substeps, extrapolation or accretion: each
    radiative point, and the convection zone as a whole, decays with the
    geometric mean of the start and end rates. Points left behind by a
    retreating zone get LIBURN's blended rate.

    Translates: liburn2.f#LIBURN2 (09c3b9b)

    Parameters
    ----------
    delts : float
        The time step.
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

    Also changes ``newrat_mod`` (the bottom point's rates) and
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
    ratli6 = np.zeros(json + 1)
    ratli7 = np.zeros(json + 1)
    ratbe9 = np.zeros(json + 1)
    fx = 0.0

    # The base of the convection zone, as in LIBURN.
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

    # The radiative interior, below both bases: skip, not leave, points
    # with a negligible Be9 rate (label 50 is inside the loop).
    minj = min(jenv, jenv0)
    maxj = max(jenv, jenv0)
    for i in range(1, minj):
        if rbe9[i] <= 1.0e-32 or rbe90[i] <= 1.0e-32:
            continue
        if (hcomp[13, i] < 1.0e-24 and hcomp[14, i] < 1.0e-24
                and hcomp[15, i] < 1.0e-24):
            continue
        if ht[i] > 7.0:
            hcomp[13, i] = 0.0
            hcomp[14, i] = 0.0
            hcomp[15, i] = 0.0
            continue
        # The geometric mean of the start and end rates.
        tli6 = 0.5 * (math.log(rli6[i]) + math.log(rli60[i]))
        tli7 = 0.5 * (math.log(rli7[i]) + math.log(rli70[i]))
        tbe9 = 0.5 * (math.log(rbe9[i]) + math.log(rbe90[i]))
        ratli6[i] = delts * math.exp(tli6)
        ratli7[i] = delts * math.exp(tli7)
        ratbe9[i] = delts * math.exp(tbe9)
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
        if hcomp[13, i] < 1.0e-24:
            hcomp[13, i] = 0.0
        if hcomp[14, i] < 1.0e-24:
            hcomp[14, i] = 0.0
        if hcomp[15, i] < 1.0e-24:
            hcomp[15, i] = 0.0

    # The convection zone.
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

    # The bottom point's mass and rates adjusted for the offset fx
    # (the rates are not restored; LIBURN's question 8).
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
    fli6 = fli60
    fli7 = fli70
    fbe9 = fbe90
    # The geometric mean of the zone's start and end rates.
    czli6 = 0.5 * (bli6 + bli60)
    czli7 = 0.5 * (bli7 + bli70)
    czbe9 = 0.5 * (bbe9 + bbe90)
    dli6 = delts * math.exp(czli6)
    dli7 = delts * math.exp(czli7)
    dbe9 = delts * math.exp(czbe9)
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
    for i in range(maxj, m + 1):
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
    # Question (see LIBURN's page, question 7): dmass < 0, frad in (-1, 0).
    dmass = hsbeg - hsend
    for i in range(jenv0, jenv):
        # Question (see the page): ends the routine, not the point.
        if rbe9[i] <= 1.0e-32:
            return  # label 200
        frad = (hs1[i] - hsbeg) / dmass
        dli6 = delts * math.exp(frad * math.log(rli6[i]) + (1.0 - frad) * bli60)
        dli7 = delts * math.exp(frad * math.log(rli7[i]) + (1.0 - frad) * bli70)
        dbe9 = delts * math.exp(frad * math.log(rbe9[i]) + (1.0 - frad) * bbe90)
        # Question (see the page): no overflow guard here.
        hcomp[13, i] = hcomp[13, i] / math.exp(dli6)
        if hcomp[13, i] < 1.0e-24:
            hcomp[13, i] = 0.0
        hcomp[14, i] = hcomp[14, i] / math.exp(dli7)
        if hcomp[14, i] < 1.0e-24:
            hcomp[14, i] = 0.0
        hcomp[15, i] = hcomp[15, i] / math.exp(dbe9)
        if hcomp[15, i] < 1.0e-24:
            hcomp[15, i] = 0.0
