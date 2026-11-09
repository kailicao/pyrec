"""Python version of ovrot.f (YREC 7be4133): convection flags for rotation
and mixing, overshoot included."""

import numpy as np

from pyrec.convec import convec

# Locals kept between calls (SAVE).
mxzon0 = np.zeros((13, 3), dtype=int)
jcore = 0
jenv = 0
nzone0 = 0


def ovrot(hcomp: np.ndarray, hd: np.ndarray, hp: np.ndarray, hr: np.ndarray, hs: np.ndarray, ht: np.ndarray,
          lc: np.ndarray, m: int, lcz: np.ndarray, mrzone: np.ndarray, mxzone: np.ndarray, nrzone: int,
          nzone: int) -> tuple[int, int]:
    """Call CONVEC for the mixed zones with overshoot, then mark every
    shell of a radiative zone not convective and every shell of a mixed
    zone convective in lcz.

    Translates: ovrot.f#OVROT (7be4133)

    Arguments: the model hcomp, hd, hp, hr, hs, ht; lc, the Schwarzschild
    flags; m, the number of shells; lcz, set in place; mrzone[13+1, 2+1],
    mxzone[12+1, 2+1], the radiative and mixed zones, set in place by
    CONVEC; nrzone, nzone (set). Returns (nrzone, nzone).

    """

    global jcore, jenv, nzone0

    # The convective regions with overshoot (Fortran lines 21–23).
    jcore, jenv, nrzone, nzone, nzone0 = convec(
        hcomp, hd, hp, hr, hs, ht, lc, m, mrzone, mxzone, mxzon0,
        jcore, jenv, nrzone, nzone, nzone0)

    # The radiative zones (Fortran lines 24–29). The guard i > 0 cannot
    # fail; an upper bound beyond m (CONVEC's suspected bug, see its page)
    # is not caught.
    for j in range(1, nrzone + 1):
        for i in range(mrzone[j, 1], mrzone[j, 2] + 1):
            if i > 0:
                lcz[i] = False

    # The mixed zones (Fortran lines 30–34).
    for j in range(1, nzone + 1):
        for i in range(mxzone[j, 1], mxzone[j, 2] + 1):
            lcz[i] = True

    return nrzone, nzone
