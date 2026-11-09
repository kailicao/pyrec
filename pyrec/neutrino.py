"""Python version of neutrino.f (YREC 7be4133): the neutrino energy-loss
rate for a four-species composition (Grant Newsham 9/06)."""

import numpy as np

from pyrec.nulosses import nulosses

# PARAMETER (Fortran lines 18–19).
IONMAX = 4


def neutrino(temp: float, den: float, x: float, y: float, zc: float, zo: float, snu: float, dsnudt: float,
             dsnudd: float) -> tuple[float, float, float]:
    """Set up hydrogen, helium, carbon 12 and oxygen 16 and return the
    neutrino energy-loss rate with its derivatives in temperature and
    density, from NULOSSES.

    Translates: neutrino.f#NEUTRINO (7be4133)

    Arguments: temp (K) and den (g/cm^3), unlogged; the mass fractions x,
    y, zc, zo; snu, dsnudt and dsnudd are passed in with any value.

    Returns (snu, dsnudt, dsnudd).

    """

    # The four species: H, He, C12, O16 (Fortran lines 20–55).
    # Question (see the page): the other metals are left out, and the four
    # fractions are not renormalized.
    xmass = np.zeros(IONMAX + 1)
    ymass = np.zeros(IONMAX + 1)
    aion = np.array([0.0, 1.0, 4.0, 12.0, 16.0])
    zion = np.array([0.0, 1.0, 2.0, 6.0, 8.0])
    xmass[1] = x
    xmass[2] = y
    xmass[3] = zc
    xmass[4] = zo

    # The losses (Fortran line 58).
    snu, dsnudt, dsnudd = nulosses(temp, den, snu, xmass, ymass, aion, zion,
                                   dsnudt, dsnudd)
    return snu, dsnudt, dsnudd
