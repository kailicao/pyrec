"""Python version of nulosses.f (YREC 09c3b9b): the neutrino energy-loss
rate of a mixture, from Timmes' code (Grant Newsham 9/06)."""

from pyrec.sneut import sneut

IONMAX = 4


def nulosses(temp, den, snu, xmass, ymass, aion, zion, dsnudt, dsnudd):
    """The mean nucleon number and charge of the mixture (azbar), then the
    neutrino loss rate and its derivatives in T and rho (sneut).

    Translates: nulosses.f#NULOSSES (09c3b9b)

    Arguments: temp (K) and den (g/cm^3), unlogged; the mass fractions
    xmass, atomic masses aion and charges zion of the IONMAX species
    (index 0 unused); ymass, set in place to the molar abundances; snu,
    dsnudt and dsnudd are passed in with any value.

    Returns (snu, dsnudt, dsnudd).
    """
    abar = zbar = 0.0
    abar, zbar = azbar(xmass, aion, zion, IONMAX, ymass, abar, zbar)
    dsnuda = dsnudz = 0.0
    # Question (see the page): the Fermi-temperature cutoff inside sneut.
    snu, dsnudt, dsnudd, dsnuda, dsnudz = sneut(
        temp, den, abar, zbar, snu, dsnudt, dsnudd, dsnuda, dsnudz)
    return snu, dsnudt, dsnudd


def azbar(xmass, aion, zion, ionmax, ymass, abar, zbar):
    """The mean nucleon number abar = 1/sum(X_i/A_i) and mean charge
    zbar = abar*sum(Z_i X_i/A_i) of the mixture.

    Translates: nulosses.f#AZBAR (09c3b9b)

    Arguments: xmass, aion, zion (index 0 unused); ionmax, the number of
    species; ymass, set in place to X_i/A_i; abar and zbar are passed in
    with any value.

    Returns (abar, zbar).
    """
    zbarxx = 0.0
    ytot1 = 0.0
    for i in range(1, ionmax + 1):
        ymass[i] = xmass[i] / aion[i]
        ytot1 = ytot1 + ymass[i]
        zbarxx = zbarxx + zion[i] * ymass[i]
    abar = 1.0 / ytot1
    zbar = zbarxx * abar
    return abar, zbar
