"""Python version of nulosses.f (YREC 7be4133): the neutrino energy-loss
rate of a mixture, from Timmes' code (Grant Newsham 9/06)."""

import numpy as np

from pyrec.sneut import sneut

# PARAMETER (Fortran lines 37–38).
IONMAX = 4


def nulosses(temp: float, den: float, snu: float, xmass: np.ndarray, ymass: np.ndarray, aion: np.ndarray,
             zion: np.ndarray, dsnudt: float, dsnudd: float) -> tuple[float, float, float]:
    """The mean nucleon number and charge of the mixture (azbar), then the
    neutrino loss rate and its derivatives in T and rho (sneut).

    Translates: nulosses.f#NULOSSES (7be4133)

    Arguments: temp (K) and den (g/cm^3), unlogged; the mass fractions
    xmass, atomic masses aion and charges zion of the IONMAX species
    (index 0 unused); ymass, set in place to the molar abundances; snu,
    dsnudt and dsnudd are passed in with any value.

    Returns (snu, dsnudt, dsnudd).

    """

    # abar and zbar (Fortran lines 45–47).
    abar = zbar = 0.0
    abar, zbar = azbar(xmass, aion, zion, IONMAX, ymass, abar, zbar)

    # The neutrino losses (Fortran lines 51–53); the derivatives in abar
    # and zbar are dropped.
    # Question (see the page): the Fermi-temperature cutoff inside sneut.
    dsnuda = dsnudz = 0.0
    snu, dsnudt, dsnudd, dsnuda, dsnudz = sneut(
        temp, den, abar, zbar, snu, dsnudt, dsnudd, dsnuda, dsnudz)
    return snu, dsnudt, dsnudd


def azbar(xmass: np.ndarray, aion: np.ndarray, zion: np.ndarray, ionmax: int, ymass: np.ndarray, abar: float,
          zbar: float) -> tuple[float, float]:
    """The mean nucleon number abar = 1/sum(X_i/A_i) and mean charge
    zbar = abar*sum(Z_i X_i/A_i) of the mixture.

    Translates: nulosses.f#AZBAR (7be4133)

    Arguments: xmass, aion, zion (index 0 unused); ionmax, the number of
    species; ymass, set in place to X_i/A_i; abar and zbar are passed in
    with any value.

    Returns (abar, zbar).

    """

    # Fortran lines 87–95.
    zbarxx = 0.0
    ytot1 = 0.0
    for i in range(1, ionmax + 1):
        ymass[i] = xmass[i] / aion[i]
        ytot1 = ytot1 + ymass[i]
        zbarxx = zbarxx + zion[i] * ymass[i]
    abar = 1.0 / ytot1
    zbar = zbarxx * abar
    return abar, zbar
