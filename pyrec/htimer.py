"""Python version of htimer.f (YREC 09c3b9b): the time step of the next
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

# SAVE: the work arrays passed to XTIME and YTIME, and EG (YTIME's).
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


def htimer(delts, deltsh, m, hd, hl, hs1, hs2, ht, hcomp, jcore, jxmid,
           tlumx, dage, ddage, nk, hp, hr, omega, dwmax, jxbeg, teffl):
    """Set the next time step: the minimum of the steps proposed by
    PTIME (structure), WTIME (rotation), XTIME (hydrogen burning), YTIME
    (helium burning) and ENTIME (the envelope triangle), and of
    atime[13] times the previous step; a fixed step for the run, or a
    negative step (relaxation), skips the criteria; the step is cut so
    that the age does not pass the run's final age.

    Translates: htimer.f#HTIMER (09c3b9b)

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

    relax_or_fixed = False
    # Question (see the page): a negative step marks a model being relaxed.
    if delts < 0.0:
        relax_or_fixed = True  # GOTO 310
    elif run.lsetdt[nk]:
        # A fixed step for this run.
        deltsh = nml.setdt[nk] * c.csecyr
        delts = deltsh
        relax_or_fixed = True  # GOTO 310

    if not relax_or_fixed:
        # Question (see the page): why log Tc > 7.1 with negative
        # gravitational luminosity, and for the rest of the job.
        if nml.lptime:
            if ht[1] > 7.1 and tlumx[7] < 0.0:
                print(f'LPTIME SET FALSE - TC {ht[1]:7.4f}'
                      f' EGRAV {tlumx[7]:10.2E}')
                nml.lptime = False

        # The structure criterion.
        if nml.lptime:
            pdelt = ptime(delts, hl, hp, hr, ht, m, 0.0)
        else:
            pdelt = 1.0e20
        # A zero time step in the model will fritz out the code.
        if pdelt < 1.0:
            pdelt = 1.0e20

        # The rotation criterion.
        if nml.lrot:
            deltsw, dwmax = wtime(delts, m, omega, 0.0, dwmax)
        else:
            deltsw = 1.0e20
        if deltsw < 1.0:
            deltsw = 1.0e20

        # Question (see the page): the burning proposals are not checked
        # against values below 1 s; the thresholds are hard-coded.
        # Hydrogen burning.
        hydlum = tlumx[1] + tlumx[2] + tlumx[3] + tlumx[4]
        if hydlum > 1.0e-34:
            deltsh = xtime(hd, hcomp, hl, hs1, hs2, ht, hydlum, jcore, jxmid,
                           m, deltsh, hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8,
                           hr9, hr10, hr11, hr12, hr13, hf1, hf2)
        else:
            deltsh = 1.0e20

        # Helium burning.
        if tlumx[5] > 1.0e-34:
            deltsy = ytime(eg, hcomp, hd, hl, hs1, ht, jcore, m, 0.0, hr1,
                           hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11,
                           hr12, hr13, hf1, hf2, jxbeg)
        else:
            deltsy = 1.0e20

        # The envelope triangle.
        if nml.ltrist and not nml.lptime:
            edelt = entime(delts, hl, teffl, m, 0.0)
        else:
            edelt = 1.0e20

        # The step may grow by at most atime[13].
        if delts > 1.0:
            delts = nml.atime[13] * delts
            delts = min(delts, deltsh, deltsy, deltsw, pdelt, edelt)
        else:
            delts = min(deltsh, deltsy, deltsw, pdelt, edelt)
        deltsh = delts

    # label 310
    delts = abs(delts)
    ddage = delts / c.csecyr
    # Do not pass the run's final age (endage in years, dage in Gyr).
    if run.lendag[nk] and nml.endage[nk] > 0.0:
        tleft = nml.endage[nk] - dage * 1.0e9
        if tleft < ddage:
            ddage = tleft
            deltsh = ddage * c.csecyr
            delts = deltsh
    return delts, deltsh, ddage, dwmax
