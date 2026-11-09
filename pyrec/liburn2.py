"""Python version of liburn2.f (YREC development, 7be4133)."""

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


def liburn2(delts: float, hcomp: np.ndarray, hr: np.ndarray, hs1: np.ndarray, hs2: np.ndarray, ht: np.ndarray,
            jnenv: int, jnenv0: int, m: int) -> None:
    """Burn Li6, Li7 and Be9 over one step in one go (cheap LIBURN).

    The same as LIBURN without substeps, extrapolation or accretion: each
    radiative point, and the convection zone as a whole, decays with the
    geometric mean of the start and end rates. Points left behind by a
    retreating zone get LIBURN's blended rate.

    Translates: liburn2.f#LIBURN2 (7be4133)

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
    ``prevcz_mod`` (``rbase0``, ``lbliold``, ``jenvos``). Nothing is
    returned.

    """

    # Module data (Fortran lines 12–39).
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

    # Local arrays (Fortran lines 76–80). SAVE in the Fortran: no local is
    # read before it is set in a call.
    ratli6 = np.zeros(json + 1)
    ratli7 = np.zeros(json + 1)
    ratbe9 = np.zeros(json + 1)
    fx = 0.0

    # The base of the convection zone, as in LIBURN (Fortran lines 88–153).
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
            # Envelope overshoot: the start of the step (Fortran lines
            # 112–128).
            if prevcz_mod.rbase0 == 0.0:
                # Suspected bug (see the page): reads hro[0] when jnenv0 = 1
                # (0 here, outside the array in the Fortran); LIBURN guards
                # this.
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

            # The end of the step, from the true base (Fortran lines
            # 129–147).
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
                    # Suspected bug (see the page): HRO, the old radius, here.
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
    # with a negligible Be9 rate (label 50 is inside the loop) (Fortran
    # lines 154–196).
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

        # The geometric mean of the start and end rates, and the burning
        # (Fortran lines 167–194).
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

    # The convection zone; skipped if it does not burn (Fortran lines
    # 197–200).
    if rbe90[jenv0] <= 1.0e-32 or rbe9[jenv] <= 1.0e-32:
        return  # label 200

    # Mean abundances over the zone at the start of the step (loop 65,
    # Fortran lines 201–215).
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

    # Always the zone's mean start rates: BURSMIX restarts from the
    # start-of-step model, so no end rate from an earlier call is reused
    # (Fortran lines 216–232).
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

    # The bottom point's mass and rates adjusted for the offset fx
    # (Fortran lines 233–242).
    # Suspected bug (see the page): the rates are not restored afterwards.
    hs2sav = hs2[jenv]
    if 1 < jenv < m:
        hs2[jenv] = hs2[jenv] + fx * (hs1[jenv] - hs1[jenv - 1])
        rli6[jenv] = rli6[jenv] + 0.5 * fx * (rli6[jenv - 1] - rli6[jenv])
        rli7[jenv] = rli7[jenv] + 0.5 * fx * (rli7[jenv - 1] - rli7[jenv])
        rbe9[jenv] = rbe9[jenv] + 0.5 * fx * (rbe9[jenv - 1] - rbe9[jenv])

    # Mean abundances and rates over the zone at the end, then the mass
    # restored (loop 70, Fortran lines 243–271).
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

    # A deepened zone starts from the end-of-step mean (Fortran lines
    # 272–278).
    if jenv < jenv0:
        fli60 = fli6
        fli70 = fli7
        fbe90 = fbe9

    # Burn the zone once with the geometric mean of its start and end
    # rates (Fortran lines 279–307).
    fli6 = fli60
    fli7 = fli70
    fbe9 = fbe90
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

    # The zone's new abundances, and LBLIOLD cleared (Fortran lines
    # 308–316).
    # Suspected bug (see the page): from maxj, so points newly in a
    # deepened zone keep their old abundances.
    for i in range(maxj, m + 1):
        hcomp[13, i] = fli6
        hcomp[14, i] = fli7
        hcomp[15, i] = fbe9
    prevcz_mod.lbliold = False

    # The region left by a retreating zone: burned once, blended rate
    # (Fortran lines 317–347).
    if jenv <= jenv0:
        return  # label 200
    if jenv0 > 1:
        hsbeg = 0.5 * (hs1[jenv0] + hs1[jenv0 - 1])
    else:
        hsbeg = 0.0
    hsend = 0.5 * (hs1[jenv] + hs1[jenv - 1])
    # Suspected bug (see the page): dmass < 0, so frad lies in (-1, 0),
    # and the weights are the other way round from the comment.
    dmass = hsbeg - hsend
    for i in range(jenv0, jenv):
        # Question (see the page): ends the routine, not the point.
        if rbe9[i] <= 1.0e-32:
            return  # label 200
        frad = (hs1[i] - hsbeg) / dmass
        dli6 = delts * math.exp(frad * math.log(rli6[i]) + (1.0 - frad) * bli60)
        dli7 = delts * math.exp(frad * math.log(rli7[i]) + (1.0 - frad) * bli70)
        dbe9 = delts * math.exp(frad * math.log(rbe9[i]) + (1.0 - frad) * bbe90)
        # Suspected bug (see the page): no overflow guard here (math.exp
        # raises OverflowError where the Fortran's EXP overflows).
        hcomp[13, i] = hcomp[13, i] / math.exp(dli6)
        if hcomp[13, i] < 1.0e-24:
            hcomp[13, i] = 0.0
        hcomp[14, i] = hcomp[14, i] / math.exp(dli7)
        if hcomp[14, i] < 1.0e-24:
            hcomp[14, i] = 0.0
        hcomp[15, i] = hcomp[15, i] / math.exp(dbe9)
        if hcomp[15, i] < 1.0e-24:
            hcomp[15, i] = 0.0
