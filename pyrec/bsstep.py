"""Python version of bsstep.f (YREC 7be4133): one Bulirsch-Stoer step."""

from collections.abc import Callable

import numpy as np

from pyrec import intpar_mod
from pyrec import iounits_mod
from pyrec.mmid import mmid
from pyrec.ratext import ratext

# Fortran line 26.
ONE = 1.0
SHRINK = 0.95
GROW = 1.2

# DATA NSEQ: the number of midpoint substeps of each sequence, index 0
# unused (Fortran line 30).
NSEQ = np.array([0, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96])


def bsstep(y: np.ndarray, dydx: np.ndarray, nv: int, x0: float, htry: float, eps: float, yscal: np.ndarray, hdid: float,
           hnext: float, deriv: Callable, b: float, fpl: float, ftl: float, gl: float, latmo: bool, lderiv: bool,
           locond: bool, lsave: bool, rl: float, teffl: float, x: float, z: float, kount: int, ksaha: int,
           err: np.ndarray) -> tuple[float, float, float, float, int, int]:
    """One Bulirsch-Stoer step (modified midpoint plus rational
    extrapolation) of the 3 envelope variables, with step-size control.

    Translates: bsstep.f#BSSTEP (7be4133)

    Arguments: y, dydx (arrays of 3, index 0 unused) at x0; the trial step
    htry; the tolerance eps; the scales yscal; the derivative routine deriv
    (qatm or qenv); the stellar quantities b ... ksaha passed on to deriv.
    Arrays changed in place: y (advanced), err (the scaled error of each
    variable).

    Returns (x0, hdid, hnext, rl, kount, ksaha): the advanced independent
    variable, the step taken and the suggested next step, then rl and the
    counters as deriv (through mmid) leaves them. The hdid and hnext passed
    in are not used.

    """

    # The local arrays (Fortran lines 27–28).
    yerr = np.zeros(4)
    ysav = np.zeros(4)
    dysav = np.zeros(4)
    yseq = np.zeros(4)

    # Save the starting values (Fortran lines 32–37).
    h = htry
    xsav = x0
    for i in range(1, nv + 1):
        ysav[i] = y[i]
        dysav[i] = dydx[i]

    # Up to IMAX sequences of midpoint steps, extrapolated to zero step;
    # accept the step once the error is below eps (Fortran lines 38–61).
    while True:  # label 20
        for i in range(1, intpar_mod.imax + 1):
            rl, kount, ksaha = mmid(ysav, dysav, nv, xsav, h, NSEQ[i], yseq,
                                    deriv, b, fpl, ftl, gl, latmo, lderiv,
                                    locond, lsave, rl, teffl, x, z, kount,
                                    ksaha)
            xest = (h / NSEQ[i])**2
            ratext(i, xest, yseq, y, yerr, nv, intpar_mod.nuse)
            errmax = 0.0
            for j in range(1, nv + 1):
                errmax = max(errmax, abs(yerr[j] / yscal[j]))
                err[j] = abs(yerr[j] / yscal[j])
            errmax = errmax / eps
            if errmax < ONE:
                x0 = x0 + h
                hdid = h
                if i == intpar_mod.nuse:
                    hnext = h * SHRINK
                elif i == intpar_mod.nuse - 1:
                    hnext = h * GROW
                else:
                    hnext = h * float(NSEQ[intpar_mod.nuse - 1]) / float(NSEQ[i])
                return x0, hdid, hnext, rl, kount, ksaha

        # All IMAX sequences failed: shrink the step and retry (Fortran
        # lines 62–72).
        h = 0.25 * h / 2.0**int((intpar_mod.imax - intpar_mod.nuse) // 2)
        # Suspected bug (see the page): this tests x (the hydrogen
        # abundance), not the independent variable x0.
        if x + h == x:
            text = (f" ERROR IN BSSTEP: STEP SIZE UNDERFLOW AT X{x:14.6E}"
                    f" H{h:12.3E}\n RUN STOPPED")  # format 900
            print(text)
            iounits_mod.write(iounits_mod.ishort, text)
            raise SystemExit(1)    # STOP 1, the run-time error code
