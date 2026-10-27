"""Python version of xtime.f (YREC 09c3b9b): the time step allowed by
hydrogen burning."""

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec.eqburn import eqburn
from pyrec.rates import rates


def xtime(hd, hcomp, hl, hs1, hs2, ht, hydlum, jcore, jxmid, m,
          deltsh, hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11,
          hr12, hr13, hf1, hf2):
    """The hydrogen-burning time step: while the central X is at least
    atime(1), the time to burn a given amount of the core's hydrogen;
    afterwards the time to burn a fraction atime(6) of the hydrogen, and,
    if X falls at the shell's midpoint jxmid, at most the time to deplete
    atime(7) of X there.

    Translates: xtime.f#XTIME (09c3b9b)

    Arguments: the model hd, hcomp[species, shell], hl, hs1, hs2, ht
    (index 0 unused); the hydrogen-burning luminosity hydlum; jcore, the
    convective core's edge; jxmid, the shell's midpoint; m, the number of
    shells; deltsh is passed in with any value; hr1..hr13, hf1, hf2, the
    rate arrays that rates fills (at element jxmid).

    Returns deltsh (s).
    """
    atime = nmlpar_mod.atime
    c = constants_mod
    if hcomp[1, 1] >= atime[1]:
        # Question (see the page): X at jcore, not at the centre.
        dx = min(atime[2], atime[3] * hcomp[1, jcore])
        deltsh = (6.00e18 / c.clsun) * (hs1[jcore] / hl[jcore]) * dx
        return deltsh
    # Question (see the page): one solar mass, not the star's.
    dx = atime[6] * hcomp[1, m] * (c.cmsun / c.clsun)
    deltsh = 6.00e18 * dx / hydlum
    ibegin = jxmid
    iend = jxmid
    dl = hd[iend]
    tl = ht[iend]
    x = hcomp[1, iend]
    y = hcomp[2, iend]
    z = hcomp[3, iend]
    xhe3 = hcomp[4, iend]
    xc12 = hcomp[5, iend]
    xc13 = hcomp[6, iend]
    xn14 = hcomp[7, iend]
    xo16 = hcomp[9, iend]
    xo18 = hcomp[11, iend]
    # (HCOMP 8, 10, 12-15 are read into XN15, XO17, XH2, XLI6, XLI7, XBE9
    # in the Fortran, but not used.)
    rates(dl, tl, x, y, xhe3, xc12, xc13, xn14, xo16, xo18, iend,
          hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12,
          hr13, hf1, hf2)
    dcdt = dodt = dxdt = dydt = 0.0
    dcdt, dodt, dxdt, dydt, xc12, xo16, x, z = eqburn(
        hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr10, hr12, hs2, ht, ibegin,
        iend, dcdt, dodt, dxdt, dydt, xc12, xo16, x, z)
    if dxdt < 0.0 and atime[7] > 0.0:
        # Question (see the page): dX/dt taken as per Gyr.
        deltsh2 = abs(c.csecyr * 1.0e9 * atime[7] / dxdt)
        deltsh = min(deltsh, deltsh2)
    return deltsh
