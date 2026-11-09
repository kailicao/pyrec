"""Python version of htimer.f (YREC 7be4133): the time step of the next
model, the most restrictive of several criteria."""

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod
from pyrec.entime import entime
from pyrec.ptime import ptime
from pyrec.wtime import wtime
from pyrec.xtime import xtime
from pyrec.ytime import ytime

# SAVE: the work arrays passed to XTIME and YTIME, and EG (YTIME's)
# (Fortran lines 56–61).
_json = dimensions_mod.json
eg = np.zeros(6 + 1)
hr1 = np.zeros(_json + 1)
hr2 = np.zeros(_json + 1)
hr3 = np.zeros(_json + 1)
hr4 = np.zeros(_json + 1)
hr5 = np.zeros(_json + 1)
hr6 = np.zeros(_json + 1)
hr7 = np.zeros(_json + 1)
hr8 = np.zeros(_json + 1)
hr9 = np.zeros(_json + 1)
hr10 = np.zeros(_json + 1)
hr11 = np.zeros(_json + 1)
hr12 = np.zeros(_json + 1)
hr13 = np.zeros(_json + 1)
hf1 = np.zeros(_json + 1)
hf2 = np.zeros(_json + 1)


def htimer(delts: float, deltsh: float, m: int, hd: np.ndarray, hl: np.ndarray, hs1: np.ndarray, hs2: np.ndarray,
           ht: np.ndarray, hcomp: np.ndarray, jcore: int, jxmid: int, tlumx: np.ndarray, dage: float, ddage: float,
           nk: int, hp: np.ndarray, hr: np.ndarray, omega: np.ndarray, dwmax: float, jxbeg: int,
           teffl: float) -> tuple[float, float, float, float]:
    """Set the next time step: the minimum of the steps proposed by
    PTIME (structure), WTIME (rotation), XTIME (hydrogen burning), YTIME
    (helium burning) and ENTIME (the envelope triangle), and of
    atime[13] times the previous step; a fixed step for the run, or a
    negative step (relaxation), skips the criteria; the step is cut so
    that the age does not pass the run's final age.

    Translates: htimer.f#HTIMER (7be4133)

    Arguments: delts (the previous step, s), deltsh (left as passed on
    the relaxation branch), the model (m, hd, hl, hs1, hs2, ht, hcomp,
    hp, hr, omega, teffl), jcore, jxmid, jxbeg (core and shell indices),
    tlumx (luminosities by source, solar units), dage (age, Gyr), nk
    (the run); ddage is only set; dwmax is set by WTIME when rotating.

    Returns (delts, deltsh, ddage, dwmax).

    """

    nml = nmlpar_mod
    run = runloop_mod
    c = constants_mod

    # A negative step marks a model being relaxed; a fixed step for this run
    # (Fortran lines 89–95).
    relax_or_fixed = False
    # Question (see the page): is the negative branch still needed?
    if delts < 0.0:
        relax_or_fixed = True  # GOTO 310
    elif run.lsetdt[nk]:
        deltsh = nml.setdt[nk] * c.csecyr
        delts = deltsh
        relax_or_fixed = True  # GOTO 310

    if not relax_or_fixed:
        # Switch the structure criterion off near the ZAMS (Fortran lines
        # 96–104).
        # Question (see the page): it stays off for the rest of the job.
        if nml.lptime:
            if ht[1] > 7.1 and tlumx[7] < 0.0:
                print(f'LPTIME SET FALSE - TC {ht[1]:7.4f}'
                      f' EGRAV {tlumx[7]:10.2E}')
                nml.lptime = False

        # The structure criterion; a zero time step in the model will fritz
        # out the code (Fortran lines 109–116).
        if nml.lptime:
            pdelt = ptime(delts, hl, hp, hr, ht, m, 0.0)
        else:
            pdelt = 1.0e20
        if pdelt < 1.0:
            pdelt = 1.0e20

        # The rotation criterion (Fortran lines 124–131).
        if nml.lrot:
            deltsw, dwmax = wtime(delts, m, omega, 0.0, dwmax)
        else:
            deltsw = 1.0e20
        if deltsw < 1.0:
            deltsw = 1.0e20

        # Hydrogen burning (Fortran lines 134–143).
        # Question (see the page): the burning proposals are not checked
        # against values below 1 s; the thresholds are hard-coded.
        hydlum = tlumx[1] + tlumx[2] + tlumx[3] + tlumx[4]
        if hydlum > 1.0e-34:
            deltsh = xtime(hd, hcomp, hl, hs1, hs2, ht, hydlum, jcore, jxmid,
                           m, deltsh, hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8,
                           hr9, hr10, hr11, hr12, hr13, hf1, hf2)
        else:
            deltsh = 1.0e20

        # Helium burning (Fortran lines 146–153).
        if tlumx[5] > 1.0e-34:
            deltsy = ytime(eg, hcomp, hd, hl, hs1, ht, jcore, m, 0.0, hr1,
                           hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11,
                           hr12, hr13, hf1, hf2, jxbeg)
        else:
            deltsy = 1.0e20

        # The envelope triangle (Fortran lines 158–162).
        # Suspected bug (see the page): edelt is not checked against values
        # below 1 s, as pdelt and deltsw are.
        if nml.ltrist and not nml.lptime:
            edelt = entime(delts, hl, teffl, m, 0.0)
        else:
            edelt = 1.0e20

        # The minimum; the step may grow by at most atime[13] (Fortran
        # lines 171–182).
        if delts > 1.0:
            delts = nml.atime[13] * delts
            delts = min(delts, deltsh, deltsy, deltsw, pdelt, edelt)
        else:
            delts = min(deltsh, deltsy, deltsw, pdelt, edelt)

        # Fortran line 211.
        deltsh = delts

    # Label 310: do not pass the run's final age (endage in years, dage in
    # Gyr) (Fortran lines 212–226).
    delts = abs(delts)
    ddage = delts / c.csecyr
    if run.lendag[nk] and nml.endage[nk] > 0.0:
        tleft = nml.endage[nk] - dage * 1.0e9
        if tleft < ddage:
            ddage = tleft
            deltsh = ddage * c.csecyr
            delts = deltsh

    return delts, deltsh, ddage, dwmax
