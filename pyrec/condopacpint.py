"""Python version of condopacpint.f (YREC 7be4133): the conductive opacity
of an H/He/metal mixture from Potekhin's conductivities."""

import math
import numpy as np

from pyrec.condopacp import condopacp

# DATA: the atomic masses of H1, He4 and the metals (taken as oxygen).
AH1 = 1.008
AHE4 = 4.004
AOX = 16.0


def condopacpint(dl: float, tl: float, x: float, z: float, oc: float, ocl: float, qodc: float, qotc: float,
                 fxion: np.ndarray, lcondo: bool) -> tuple[float, float, float, float, bool]:
    """The conductive opacity of a mixture: Potekhin's conductivities for
    Z_ion = 1, 2 and 8, summed harmonically with weights by ionized number
    fraction, then kappa = 16 sigma T^3 / (3 rho K); extrapolated linearly
    below log rho = -6.

    Translates: condopacpint.f#CONDOPACPINT (7be4133)

    Arguments: dl, log rho; tl, log T; x, z, the mass fractions of H and
    metals; fxion[1..3], the ionization fractions (H+, He+, He2+); oc,
    ocl, qodc, qotc, lcondo, passed in (with any value).

    Returns (oc, ocl, qodc, qotc, lcondo): the opacity and its log, its
    derivatives in log rho and log T, and True if computed (otherwise oc
    is 10**9.9998, with qodc = qotc = oc).

    """
    # The number fractions of H1, He4 and the metals (Fortran lines 24–28).
    y = 1.0 - x - z
    wth1base = (x / AH1) / (x / AH1 + y / AHE4 + z / AOX)
    wthe4base = (y / AHE4) / (x / AH1 + y / AHE4 + z / AOX)
    wtoxbase = (z / AOX) / (x / AH1 + y / AHE4 + z / AOX)

    # Weighted by ionization; the metals taken as fully ionized (Fortran lines 30–41).
    # Question (see the page): (f1 + 2 f2)/2 for helium.
    # Question (see the page): all metals as oxygen, fully ionized.
    wth1 = wth1base * fxion[1]
    wthe4 = wthe4base * (fxion[2] + 2.0 * fxion[3]) / 2.0
    wtox = wtoxbase

    # No significant conductive opacity, to begin with (Fortran lines 43–50).
    lcondo = False
    ocl = 9.9998
    oc = 10.0 ** ocl
    qodc = oc
    qotc = oc

    # Within the table's range (Fortran lines 52–54).
    rl = dl - 3.0 * (tl - 6.0)
    if tl >= 3.0 and tl <= 9.0 and dl <= 9.75:
        if dl >= -6.0:
            # The three pure compositions from the table (Fortran lines 55–58).
            condxl, qodcx, qotcx = condopacp(1.0, tl, dl, 0.0, 0.0, 0.0)
            condyl, qodcy, qotcy = condopacp(2.0, tl, dl, 0.0, 0.0, 0.0)
            condzl, qodcz, qotcz = condopacp(8.0, tl, dl, 0.0, 0.0, 0.0)

            # Their harmonic sum, log K (Fortran lines 60–66).
            condx = 10.0 ** (-condxl)
            condy = 10.0 ** (-condyl)
            condz = 10.0 ** (-condzl)
            condl = -math.log10(wth1 * condx + wthe4 * condy + wtox * condz)

            # The mixture's derivatives (Fortran lines 68–74).
            # Suspected bug (see the page): the factor condl (log K).
            qodc = condl * (wth1 * qodcx + wthe4 * qodcy + wtox * qodcz)
            # Suspected bug (see the page): the density derivatives again.
            qotc = condl * (wth1 * qodcx + wthe4 * qodcy + wtox * qodcz)

            # The conductive opacity (Fortran lines 76–80).
            ocl = -3.5194 + 3.0 * tl - dl - condl
            qodc = -1.0 - qodc
            qotc = 3.0 - qotc
            lcondo = True
            oc = 10.0 ** ocl

        if dl < -6.0 and rl >= 0.0:
            # Extrapolate from log rho = -6: the table there (Fortran lines 82–86).
            condxl, qodc, qotc = condopacp(1.0, tl, -6.0, 0.0, 0.0, 0.0)
            condyl, crap, crap2 = condopacp(2.0, tl, -6.0, 0.0, 0.0, 0.0)
            condzl, crap, crap2 = condopacp(8.0, tl, -6.0, 0.0, 0.0, 0.0)

            # Their harmonic sum (Fortran lines 88–92).
            condx = 10.0 ** (-condxl)
            condy = 10.0 ** (-condyl)
            condz = 10.0 ** (-condzl)
            # Suspected bug (see the page): no minus sign here, unlike above.
            condl = math.log10(wth1 * condx + wthe4 * condy + wtox * condz)

            # The opacity at log rho = -6, extrapolated to dl (Fortran lines 96–101).
            ocl = -3.5194 + 3.0 * tl + 6.0 - condl
            qodc = -1.0 - qodc
            qotc = 3.0 - qotc
            # Suspected bug (see the page): the step (-6 - dl) is reversed.
            ocl = ocl + qodc * (-6.0 - dl)
            oc = 10.0 ** ocl
            lcondo = True

    return oc, ocl, qodc, qotc, lcondo
