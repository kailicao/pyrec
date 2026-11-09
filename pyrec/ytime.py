"""Python version of ytime.f (YREC 7be4133): the time step allowed by
helium burning."""

from pyrec import constants_mod
import numpy as np
from pyrec import nmlpar_mod
from pyrec.engeb import engeb

# SAVE: the extended-composition abundances keep their values between
# calls, so that without LEXCOM ENGEB gets the last ones set (or 0)
# (Fortran line 61).
xh2 = 0.0
xli6 = 0.0
xli7 = 0.0
xbe9 = 0.0

# Central Y below which the core criterion is dropped in the shell branch
# (Fortran lines 31–32).
YFLOOR = 1.0e-5


def ytime(eg: np.ndarray, hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hs1: np.ndarray, ht: np.ndarray,
          jcore: int, m: int, deltsy: float, hr1: np.ndarray, hr2: np.ndarray, hr3: np.ndarray, hr4: np.ndarray,
          hr5: np.ndarray, hr6: np.ndarray, hr7: np.ndarray, hr8: np.ndarray, hr9: np.ndarray, hr10: np.ndarray,
          hr11: np.ndarray, hr12: np.ndarray, hr13: np.ndarray, hf1: np.ndarray, hf2: np.ndarray, jxbeg: int) -> float:
    """The helium-burning time step DELTSY: in the helium flash (central
    log rho >= 5), 1e15 over the triple-alpha rate at the temperature
    maximum; otherwise the time to burn a fraction of the core's helium,
    or, below ATIME(1), the smallest of the shell and core limits.

    Translates: ytime.f#YTIME (7be4133)

    Arguments: the energy rates eg(1..6) (set here, in the flash), the model
    (hcomp, hd = log rho, hl, hs1 = m, ht = log T, jcore, m, jxbeg), the
    step deltsy, and the rate work arrays hr1..hr13, hf1, hf2 passed on to
    ENGEB. Returns deltsy.

    """

    global xh2, xli6, xli7, xbe9
    n = nmlpar_mod

    if hd[1] >= 5.0:
        # Helium flash: find the temperature maximum (Fortran lines 69–78).
        htmax = 0.0
        itmax = 1
        for it in range(1, m + 1):  # label 210
            if ht[it] > htmax:
                itmax = it
                htmax = ht[it]

        # The composition there (Fortran lines 80–99).
        dl = hd[itmax]
        tl = ht[itmax]
        x = hcomp[1, itmax]
        y = hcomp[2, itmax]
        z = hcomp[3, itmax]  # noqa: F841  (unused in the Fortran too)
        xhe3 = hcomp[4, itmax]
        xc12 = hcomp[5, itmax]
        xc13 = hcomp[6, itmax]
        xn14 = hcomp[7, itmax]
        xn15 = hcomp[8, itmax]  # noqa: F841  (unused in the Fortran too)
        xo16 = hcomp[9, itmax]
        xo17 = hcomp[10, itmax]  # noqa: F841  (unused in the Fortran too)
        xo18 = hcomp[11, itmax]
        if n.lexcom:
            xh2 = hcomp[12, itmax]
            xli6 = hcomp[13, itmax]
            xli7 = hcomp[14, itmax]
            xbe9 = hcomp[15, itmax]
        iu = itmax

        # The energy rates at the temperature maximum; ENGEB's set scalars
        # are passed in with any value, and BE7MFLOC is not used here
        # (Fortran lines 104–117).
        (egg11, egg22, egg33, egg44, egg55, qed, qet, sum1, be7mfloc,
         ealpca, enu) = engeb(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dl,
                              tl, x, y, xhe3, xc12, xc13, xn14, xo16,
                              xo18, xh2, iu, hr1, hr2, hr3, hr4, hr5, hr6,
                              hr7, hr8, hr9, hr10, hr11, hr12, hr13, hf1,
                              hf2, 0.0, 0.0, 0.0)
        e = sum1  # noqa: F841  (unused in the Fortran too)
        eg[1] = egg11
        eg[2] = egg22
        eg[3] = egg33
        eg[4] = egg44
        eg[5] = egg55
        eg[6] = enu
        eca = ealpca  # noqa: F841  (unused in the Fortran too)

        # The step: 1e15 erg/g over the triple-alpha rate (Fortran lines
        # 118–119).
        # Question (see the page): the constant 1e15 and the floor 1e-22.
        if eg[5] < 1.0e-22:
            eg[5] = 1.0e-22
        deltsy = 1.0e15 / eg[5]

    else:
        # Limits on helium burning for stars outside the flash: the core
        # criterion (Fortran lines 122–127).
        # Question (see the page): 1 - Z is X + Y, not Y.
        ycore = 1.0 - hcomp[3, jcore]
        if ycore >= n.atime[1]:
            deltsy = min(n.atime[4], n.atime[5] * ycore)
            deltsy = ((5.85e17 / constants_mod.clsun) * deltsy
                      * (hs1[jcore] / hl[jcore]))

        else:
            # The shell branch: atime[14] solar masses of helium at the
            # luminosity just below the hydrogen shell (Fortran lines
            # 142–146).
            deltsy = 1.0e20
            if jxbeg > 1:
                if hl[jxbeg - 1] > 0.0:
                    deltsy = ((5.85e17 / constants_mod.clsun) * n.atime[14]
                              * constants_mod.cmsun / hl[jxbeg - 1])

            # The fraction atime[12] of Y at each helium-burning point
            # (Fortran lines 147–170).
            for j in range(1, m + 1):  # label 220
                if ht[j] < n.tcut[4] or hcomp[2, j] < n.atime[1]:
                    continue
                dl = hd[j]
                tl = ht[j]
                x = hcomp[1, j]
                y = hcomp[2, j]
                xhe3 = hcomp[4, j]
                xc12 = hcomp[5, j]
                xc13 = hcomp[6, j]
                xn14 = hcomp[7, j]
                xo16 = hcomp[9, j]
                xo18 = hcomp[11, j]
                xh2 = 0.0
                if n.lexcom:
                    xh2 = hcomp[12, j]
                iu = j
                (egg11, egg22, egg33, egg44, egg55, qed, qet, sum1,
                 be7mfloc, ealpca, enu) = engeb(
                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dl, tl, x, y,
                    xhe3, xc12, xc13, xn14, xo16, xo18, xh2, iu, hr1, hr2,
                    hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12,
                    hr13, hf1, hf2, 0.0, 0.0, 0.0)
                # Rate of Y consumption: 3-alpha and alpha captures.
                dydt = max(egg55, 0.0) / 5.85e17 + max(ealpca, 0.0) / 1.7276e18
                if dydt > 0.0:
                    deltsy = min(deltsy, n.atime[12] * y / dydt)

            # The core criterion without its cap, while helium is left at
            # the center (Fortran lines 171–173).
            if ycore > YFLOOR and hl[jcore] > 0.0:
                deltsy = min(deltsy, (5.85e17 / constants_mod.clsun)
                             * n.atime[5] * ycore * (hs1[jcore] / hl[jcore]))

    return deltsy
