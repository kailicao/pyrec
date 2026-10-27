"""Python version of lirate88.f (YREC development, 09c3b9b)."""

import math

from pyrec import constants_mod
from pyrec import dimensions_mod  # noqa: F401  (JSON, the declared size)
from pyrec import newrat_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldrat_mod

TLIM = 6.0  # DATA; Question (see the page): hard-coded cut in log T.


def lirate88(hcomp, hd, ht, m, j):
    """Proton-burning rates of Li6, Li7 and Be9 (Caughlan & Fowler 1988).

    Each rate is the Caughlan-Fowler fit times its user scale factor and
    rho X. With ``j`` = 1 the rates come from the arguments and fill
    ``newrat_mod``; with ``j`` = 2 from the old model in ``oldmod_mod`` and
    fill ``oldrat_mod``. From the first point outward where both the new
    and old log T are below TLIM, all six rate arrays are zeroed.

    Translates: lirate88.f#LIRATE88 (09c3b9b)

    Parameters
    ----------
    hcomp : ndarray, shape (16, JSON + 1)
        Composition, ``hcomp[1, i]`` = X (used with ``j`` = 1).
    hd, ht : ndarray
        log rho and log T at each point, 1..m (with ``j`` = 2, ``ht``
        serves only for the cut).
    m : int
        Number of points.
    j : int
        1 for the current model, 2 for the old one.

    Sets ``newrat_mod.rli6, rli7, rbe9`` or ``oldrat_mod.rli60, rli70,
    rbe90``. Nothing is returned.
    """
    cln = constants_mod.cln
    cc13 = constants_mod.cc13
    cc23 = constants_mod.cc23
    rli6 = newrat_mod.rli6
    rli7 = newrat_mod.rli7
    rbe9 = newrat_mod.rbe9
    rli60 = oldrat_mod.rli60
    rli70 = oldrat_mod.rli70
    rbe90 = oldrat_mod.rbe90
    hto = oldmod_mod.hto

    i = m + 1  # the Fortran's I after loop 100 runs to its end
    for ii in range(1, m + 1):  # loop 100
        if ht[ii] < TLIM and hto[ii] < TLIM:
            i = ii
            break  # to label 110
        # Question (see the page): with j = 2 the rates come from the old
        # model, so both calls in main.f may give the starting rates.
        if j == 1:
            rhox = math.exp(cln * hd[ii]) * hcomp[1, ii]
            t9 = math.exp(cln * (ht[ii] - 9.0))
        else:
            rhox = math.exp(cln * oldmod_mod.hdo[ii]) * oldmod_mod.hcompp[1, ii]
            t9 = math.exp(cln * (hto[ii] - 9.0))
        t913 = t9 ** cc13
        t923 = t913 * t913
        t943 = t923 * t923
        t953 = t943 * t913
        t932 = t9 ** 1.5
        t934 = t9 ** 7.5e-1
        # Question (see the page): the old-typo notes in the Fortran's
        # comments on these coefficients.
        fli6 = (3.73e10 / t923 * math.exp(-8.413 / t913 - (t9 / 5.50) ** 2)
                * (1.0 + 0.050 * t913 - 0.061 * t923 - 0.021 * t9
                   + 0.006 * t943 + 0.005 * t953)
                + 1.33e10 / t932 * math.exp(-17.763 / t9)
                + 1.29e9 / t9 * math.exp(-21.820 / t9))
        t9a = t9 / (1.0 + 0.759 * t9)
        c56 = 1.25 * cc23
        fli7 = (1.096e9 / t923 * math.exp(-8.472 / t913)
                - 4.830e8 * t9a ** c56 / t932 * math.exp(-8.472 / t9a ** cc13)
                + 1.06e10 / t932 * math.exp(-30.442 / t9))
        fbe91 = (1.33e+07 / t923 * math.exp(-10.359 / t913 - (t9 / 0.846) ** 2)
                 * (1.0 + 4.0e-2 * t913 + 1.52 * t923 + 4.28e-1 * t9
                    + 2.15 * t943 + 1.54 * t953)
                 + 9.64e+4 / t932 * math.exp(-3.445 / t9)
                 + 2.72e+6 / t932 * math.exp(-10.62 / t9))
        fx = (2.11e+11 / t923 * math.exp(-10.359 / t913 - (t9 / 5.2e-1) ** 2)
              * (1.0 + 4.0e-2 * t913 + 1.09 * t923 + 3.07e-1 * t9
                 + 3.21 * t943 + 2.3 * t953))
        ex = math.exp(-3.046 / t9) / t9
        fbe92 = fx + 5.79e+8 * ex + 8.5e+8 / t934 * math.exp(-5.8 / t9)
        fbe93 = fx + 4.51e+8 * ex + 6.7e+8 / t934 * math.exp(-5.16 / t9)
        fli6 = nmlpar_mod.sli6 * fli6
        fli7 = nmlpar_mod.sli7 * fli7
        fbe91 = nmlpar_mod.sbe91 * fbe91
        fbe92 = nmlpar_mod.sbe92 * fbe92
        fbe93 = nmlpar_mod.sbe93 * fbe93
        fsbe9 = fbe91 + fbe92 + fbe93
        if j == 1:
            rli6[ii] = fli6 * rhox
            rli7[ii] = fli7 * rhox
            rbe9[ii] = fsbe9 * rhox
        else:
            rli60[ii] = fli6 * rhox
            rli70[ii] = fli7 * rhox
            rbe90[ii] = fsbe9 * rhox

    # label 110: zero both sets from the cut outward (question on the page).
    for k in range(i, m + 1):  # loop 120
        rli6[k] = 0.0
        rli7[k] = 0.0
        rbe9[k] = 0.0
        rli60[k] = 0.0
        rli70[k] = 0.0
        rbe90[k] = 0.0
