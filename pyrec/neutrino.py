"""Python version of neutrino.f (YREC 09c3b9b): the neutrino energy-loss
rate for a four-species composition (Grant Newsham 9/06)."""

import numpy as np

from pyrec.nulosses import nulosses

IONMAX = 4


def neutrino(temp, den, x, y, zc, zo, snu, dsnudt, dsnudd):
    """Set up hydrogen, helium, carbon 12 and oxygen 16 and return the
    neutrino energy-loss rate with its derivatives in temperature and
    density, from NULOSSES.

    Translates: neutrino.f#NEUTRINO (09c3b9b)

    Arguments: temp (K) and den (g/cm^3), unlogged; the mass fractions x,
    y, zc, zo; snu, dsnudt and dsnudd are passed in with any value.

    Returns (snu, dsnudt, dsnudd).
    """
    xmass = np.zeros(IONMAX + 1)
    ymass = np.zeros(IONMAX + 1)
    # Question (see the page): the other metals are left out, and the four
    # fractions are not renormalized.
    aion = np.array([0.0, 1.0, 4.0, 12.0, 16.0])
    zion = np.array([0.0, 1.0, 2.0, 6.0, 8.0])
    xmass[1] = x
    xmass[2] = y
    xmass[3] = zc
    xmass[4] = zo
    snu, dsnudt, dsnudd = nulosses(temp, den, snu, xmass, ymass, aion, zion,
                                   dsnudt, dsnudd)
    return snu, dsnudt, dsnudd
