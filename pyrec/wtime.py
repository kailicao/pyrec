"""Python version of wtime.f (YREC 09c3b9b): the time step allowed by the
change of the angular velocity."""

from pyrec import nmlpar_mod
from pyrec import oldrot_mod


def wtime(delts, m, omega, deltsw, dwmax):
    """Limit the step of a rotating model by the largest relative change
    of Omega since the previous model: the step is divided by that change
    over dtwind, the factor kept within [1/atime(13), atime(13)]; with no
    change at all, the step is 1e20 (no limit).

    Translates: wtime.f#WTIME (09c3b9b)

    Arguments: delts, the previous step; m, the number of shells; omega
    (index 0 unused); deltsw and dwmax are passed in with any value.

    Returns (deltsw, dwmax).
    """
    wold = oldrot_mod.wold
    ii = 1
    # Question (see the page): 0/0 where Omega and its old value are 0.
    dwmax = 2.0 * abs(omega[ii] - wold[ii]) / (omega[ii] + wold[ii])
    iwmax = ii       # Question (see the page): found but not used.
    for i in range(ii + 1, m + 1):      # loop 50
        testw = 2.0 * abs(omega[i] - wold[i]) / (omega[i] + wold[i])
        if testw > dwmax:
            dwmax = testw
            iwmax = i
    fac = dwmax / nmlpar_mod.dtwind
    if fac == 0.0:
        deltsw = 1.0e20
        return deltsw, dwmax           # GOTO 999
    temp = nmlpar_mod.atime[13]
    if fac > temp:
        fac = temp
    if fac < 1.0 / temp:
        fac = 1.0 / temp
    deltsw = delts / fac
    return deltsw, dwmax
