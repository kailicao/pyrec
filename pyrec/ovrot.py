"""Python version of ovrot.f (YREC 09c3b9b): convection flags for rotation
and mixing, overshoot included."""

import numpy as np

from pyrec.convec import convec

# Locals kept between calls (SAVE).
mxzon0 = np.zeros((13, 3), dtype=int)
jcore = 0
jenv = 0
nzone0 = 0


def ovrot(hcomp, hd, hp, hr, hs, ht, lc, m, lcz, mrzone, mxzone,
          nrzone, nzone):
    """Call CONVEC for the mixed zones with overshoot, then mark every
    shell of a radiative zone not convective and every shell of a mixed
    zone convective in lcz.

    Translates: ovrot.f#OVROT (09c3b9b)

    Arguments: the model hcomp, hd, hp, hr, hs, ht; lc, the Schwarzschild
    flags; m, the number of shells; lcz, set in place; mrzone[13+1, 2+1],
    mxzone[12+1, 2+1], the radiative and mixed zones, set in place by
    CONVEC; nrzone, nzone (set). Returns (nrzone, nzone).
    """
    global jcore, jenv, nzone0
    jcore, jenv, nrzone, nzone, nzone0 = convec(
        hcomp, hd, hp, hr, hs, ht, lc, m, mrzone, mxzone, mxzon0,
        jcore, jenv, nrzone, nzone, nzone0)
    for j in range(1, nrzone + 1):
        for i in range(mrzone[j, 1], mrzone[j, 2] + 1):
            # Question (see the page): KC's guard against i <= 0; an upper
            # bound beyond m is not caught.
            if i > 0:
                lcz[i] = False
    for j in range(1, nzone + 1):
        for i in range(mxzone[j, 1], mxzone[j, 2] + 1):
            lcz[i] = True
    return nrzone, nzone
