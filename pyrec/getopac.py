"""Python version of getopac.f (YREC development, 7be4133)."""

import math
import sys
import numpy as np

from pyrec import comp_mod
from pyrec import newopac_mod
from pyrec import nmlpar_mod
from pyrec.condopacpint import condopacpint
from pyrec.getalex06 import getalex06
from pyrec.getopal95 import getopal95
from pyrec.gtpurz import gtpurz


def _stop(code, *items):
    """WRITE(ISHORT,*) the message, then STOP with the exit code: 3 for an
    input error (no table chosen), 1 for a run-time error.

    The unit ISHORT (from iounits_mod) is the run's short output file;
    printing to standard error stands in for it here.

    """

    print(*items, file=sys.stderr)
    sys.exit(code)


def getopac(dl: float, tl: float, x: float, z: float, o: float, ol: float, qod: float, qot: float,
            fxion: np.ndarray) -> tuple[float, float, float, float]:
    """Opacity and its logarithmic derivatives for a given rho, T, X, Z.

    Combines the molecular (Alexander 2006) and interior (OPAL 1995)
    tables with a ramp between TMOLMIN and TMOLMAX, interpolates towards
    the pure-metal table in metal-rich helium-burning cores, and adds
    conductive opacity (Potekhin, or the Hubbard-Lampe fit).

    Translates: getopac.f#GETOPAC (7be4133)

    Parameters
    ----------
    dl, tl : float
        log10 density and log10 temperature.
    x, z : float
        Hydrogen and metal mass fractions.
    o, ol, qod, qot : float
        Output only.
    fxion : ndarray, shape (4,)
        Ionization fractions ``fxion[1..3]``, passed to the conductive
        opacity.

    Returns
    -------
    o, ol, qod, qot : float
        kappa, log10 kappa, dlog kappa/dlog rho, dlog kappa/dlog T.

    """

    nml = nmlpar_mod

    # Get atmosphere (molecular) opacity (Fortran lines 78–95).
    lgotatm = False
    if tl <= nml.tmolmax:
        if nml.lalex06:
            so, sol, sqod, sqot = getalex06(dl, tl, x, z, 0.0, 0.0, 0.0, 0.0)
            lgotatm = True
        else:
            _stop(3, 'NO MOLECULAR OPACITY TABLE CHOSEN',
                  ' RUN STOPPED. X Z TL=', x, z, tl)

    # Get interior opacity if needed (else go to label 1000) (Fortran line 101).
    if not (tl < nml.tmolmin and lgotatm):

        if z > 0.1 and tl > 7.0:
            # Helium-burning region: interpolate towards the pure-Z table
            # (Fortran lines 109–138).
            if not nml.lpurez:
                # Suspected bug (see the page): the message still says
                # 5 x 10^7 K although the threshold is now log T > 7.0.
                _stop(1, ' ERROR: Z>0.10 T > 5 X 10^7 K',
                      ' NEED PURE Z TABLE TO CONTINUE. Z,LOG T=', z, tl)
            oz, olz, qodz, qotz = gtpurz(dl, tl, 0.0, 0.0, 0.0, 0.0)
            if nml.lopal95:
                # Interpolate to the maximum Z in the table.
                zit = 0.1
                o, ol, qod, qot = getopal95(dl, tl, x, zit, o, ol, qod, qot)
            else:
                _stop(3, 'NO OPACITY TABLE CHOSEN',
                      ' RUN STOPPED. X Z TL=', x, z, tl)
            slope = (ol - olz) / (zit - 1.0)
            ol = olz + (z - 1.0) * slope
            o = 10.0 ** ol
            slope = (qod - qodz) / (zit - 1.0)
            qod = qodz + (z - 1.0) * slope
            slope = (qot - qotz) / (zit - 1.0)
            qot = qotz + (z - 1.0) * slope

        # Z too large, or Z outside a single table (Fortran lines 140–147).
        # Suspected bug (see the page): with OPAL 1995 off, the second
        # test turns the "no table chosen" input error into this stop.
        elif (z > 0.12 or (abs(z - comp_mod.zenv) > nml.optol
                           and not newopac_mod.l2z and not nml.lopal95)):
            _stop(1, ' Z>0.12 T < 10^7 K',
                  ' OUTSIDE OPAL OPACITY TABLE RANGE OR Z',
                  ' OUTSIDE SINGLE TABLE USED.Z,ZENV,LOG T=',
                  z, comp_mod.zenv, tl)

        # OPAL 1995 (Fortran lines 154–155).
        elif nml.lopal95:
            o, ol, qod, qot = getopal95(dl, tl, x, z, o, ol, qod, qot)

        else:
            # Final trap: should not be able to get here (Fortran lines 162–166).
            _stop(3, 'NO OPACITY TABLE CHOSEN',
                  ' RUN STOPPED. X Z TL=', x, z, tl)

    # label 1000: ramp between surface and interior opacity (Fortran lines 173–188).
    if lgotatm and tl <= nml.tmolmax:
        if tl >= nml.tmolmin:
            rmpwt = (tl - nml.tmolmin) / (nml.tmolmax - nml.tmolmin)
            o = rmpwt * o + (1.0 - rmpwt) * so
            ol = math.log10(o)
            qod = rmpwt * qod + (1.0 - rmpwt) * sqod
            qot = rmpwt * qot + (1.0 - rmpwt) * sqot
        else:
            o = so
            ol = sol
            qod = sqod
            qot = sqot

    # Conductive opacity correction: Potekhin (Fortran lines 190–195).
    if nml.lcondopacp:
        # Potekhin conductive opacity.
        oc, ocl, qodc, qotc, lcondo = condopacpint(
            dl, tl, x, z, 0.0, 0.0, 0.0, 0.0, fxion, False)
    else:
        lcondo = False

    if not lcondo:
        # No Potekhin opacity: try Hubbard-Lampe (Fortran lines 197–213).
        # Question (see the page): this fallback is silent.
        if tl < 4.2:
            return o, ol, qod, qot
        elif dl < 2.0 * tl - 13.0:
            return o, ol, qod, qot
        else:
            cond = (math.log10(1.0 - 0.6 * x) - 14.6196
                    - (3.5853 + 0.1386 * dl) * dl
                    + (5.1324 - 0.3219 * tl) * tl + 0.3901 * dl * tl)
            oc = 10.0 ** cond
            qodc = 0.3901 * tl - 0.2772 * dl - 3.5853
            qotc = 0.3901 * dl - 0.6438 * tl + 5.1324

    # Save the radiative opacity, then add the opacities: 1/O = 1/OX + 1/OC
    # (Fortran lines 215–222).
    ox = o
    qodx = qod
    qotx = qot
    o = ox * oc / (ox + oc)
    ol = math.log10(o)
    qod = qodx + qodc - (ox * qodx + oc * qodc) / (ox + oc)
    qot = qotx + qotc - (ox * qotx + oc * qotc) / (ox + oc)

    return o, ol, qod, qot
