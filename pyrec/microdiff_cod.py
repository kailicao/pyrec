"""Python version of microdiff_cod.f (YREC 09c3b9b): settling coefficients
of one species on the equally spaced diffusion grid."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import nmlpar_mod
from pyrec.thdiff import thdiff

JSON = dimensions_mod.json
NSPEC = 4  # DATA NSPEC/4/: H, He, the species, electrons

# The Fortran's local arrays. MICRODIFF_COD has SAVE, so they keep their
# values between calls: HYDCON, filled on the hydrogen call, and the
# FACS ... ADS of points skipped on a later call (question 3 of the page).
atomwt = np.zeros(NSPEC + 1)
zspec = np.zeros(NSPEC + 1)
xfrac = np.zeros(NSPEC + 1)
concen = np.zeros(NSPEC + 1)
tcl = np.zeros((NSPEC + 1, NSPEC + 1))
tap = np.zeros(NSPEC + 1)
tat = np.zeros(NSPEC + 1)
tac = np.zeros((NSPEC + 1, NSPEC + 1))
facs = np.zeros(JSON + 1)
aps = np.zeros(JSON + 1)
ats = np.zeros(JSON + 1)
ahs = np.zeros(JSON + 1)
ads = np.zeros(JSON + 1)
hydcon = np.zeros(JSON + 1)


def microdiff_cod(npt, espec, er, eden, etem, ehqpr, edelr,
                  ecod1, ecod2, edxdr, thea, thez, j):
    """Compute the settling coefficients of species column j at each of the
    npt points of the equally spaced grid (Thoul et al. 1994 coefficients
    from THDIFF), for MICRODIFF_RUN.

    Translates: microdiff_cod.f#MICRODIFF_COD (09c3b9b)

    Arguments: npt points with radius er (Rsun), density eden, temperature
    etem (K), dlnP/dr ehqpr and the gradient edelr; espec[k, i], the mass
    fractions of H (k = 1), He (2) and the species (3); its weight thea and
    charge thez; j, the column diffused (1 for hydrogen, 3 otherwise).
    Arrays have index 0 unused.

    Sets in place: ecod1, ecod2 (the coefficients) and, on the hydrogen call
    (j = 1), edxdr, the hydrogen concentration gradient read by later calls.
    Returns nothing.
    """
    nml = nmlpar_mod
    # Set up the atomic weight and charge matrices.
    atomwt[1] = 1.008
    atomwt[2] = 4.004
    atomwt[3] = thea
    atomwt[4] = 5.486e-4
    zspec[1] = 1.0
    zspec[2] = 2.0
    zspec[3] = thez
    zspec[4] = -1.0
    # Set LN_LAMBDA, CON_RAD and CON_TEMP.
    ln_lambda = 2.2
    con_rad = 1.0 / 6.9598e10  # set but not used
    con_temp = 1.0e-7
    # Calculate diffusion coefficients for each layer.
    for i in range(1, npt + 1):  # DO 5
        xfrac[1] = espec[1, i]
        xfrac[2] = espec[2, i]
        xfrac[3] = espec[3, i]
        # Make sure XFRAC = 0.0 isn't used for diffusion coefficients.
        for k in (1, 2, 3):
            if xfrac[k] < 1.0e-24:
                xfrac[k] = 1.0e-24
        # Concentrations from mass fractions.
        zxa = 0.0
        for ii in range(1, NSPEC):
            zxa = zxa + zspec[ii] * xfrac[ii] / atomwt[ii]
        for ii in range(1, NSPEC):
            concen[ii] = xfrac[ii] / (atomwt[ii] * zxa)
        concen[NSPEC] = 1.0
        # Save the hydrogen concentration when X is diffused.
        if j == 1:
            hydcon[i] = concen[1]
        # Skip the Thoul routine where the species is absent here and at the
        # next point (the first depleted shell is kept, so that AD is right).
        # Question (see the page): only ADS is reset; FACS ... AHS keep old
        # values.
        if espec[j, i] == 0.0 and i != npt:
            if espec[j, i + 1] == 0.0:
                ads[i] = 0.0
                continue  # GOTO 5
        rho = eden[i]
        t = etem[i]
        # Electron density NE from the mass density.
        ac = 0.0
        for ii in range(1, NSPEC + 1):
            ac = ac + atomwt[ii] * concen[ii]
        ne = rho / (1.6726e-24 * ac)
        # Inter-ionic distance AO.
        ni = 0.0
        for ii in range(1, NSPEC):
            ni = ni + concen[ii] * ne
        ao = (0.23873 / ni)**constants_mod.cc13
        # Debye length LAMBDAD.
        cz = 0.0
        for ii in range(1, NSPEC + 1):
            cz = cz + concen[ii] * zspec[ii]**2
        lambdad = 6.9010 * math.sqrt(t / (ne * cz))
        lam = max(lambdad, ao)  # LAMBDA, used in the Coulomb logarithm
        # Coulomb logarithms.
        for ii in range(1, NSPEC + 1):
            for jj in range(1, NSPEC + 1):
                xij = 2.3939e3 * t * lam / abs(zspec[ii] * zspec[jj])
                tcl[ii, jj] = 0.81245 * math.log(1.0 + 0.18769 * xij**1.2)
        # The diffusion coefficients.
        thdiff(NSPEC, atomwt, zspec, xfrac, tcl, tap, tat, tac)
        hru_i = er[i]
        htu_i = t * con_temp
        # JvS 01/26: the FGRY and FGRZ factors.
        # Question (see the page): FGRZ also scales the light elements (j = 3).
        fac = hru_i**2 * htu_i**2.5 / ln_lambda
        if j == 1:
            fac = nml.fgry * hru_i**2 * htu_i**2.5 / ln_lambda
        if j == 3:
            fac = nml.fgrz * hru_i**2 * htu_i**2.5 / ln_lambda
        # The first diffusion terms for hydrogen, the third for the others;
        # stored so that the hydrogen gradient can be finished first.
        facs[i] = fac
        aps[i] = -tap[j]
        ats[i] = -tat[j] * edelr[i]
        ahs[i] = -tac[j, 1]
        ads[i] = -tac[j, j]
    # The velocity is Ap dlnP/dR + At dlnT/dR + Ax dlnX/dR + Ad dlnD/dR
    # (Thoul et al. 1994); for hydrogen the third term is the fourth.
    for i in range(1, npt + 1):  # DO 10
        fac = facs[i]
        ap = aps[i]
        at = ats[i]
        ah = ahs[i]
        ad = ads[i]
        if j == 1:
            # Hydrogen call: the gradient of ln c_H, zero at the ends, then
            # zeroed locally to remove the second self-diffusion term.
            if i == 1:
                dlncdr = 0.0
            elif i == npt:
                dlncdr = 0.0
            else:
                coni = math.log(hydcon[i])
                conip1 = math.log(hydcon[i + 1])
                conim1 = math.log(hydcon[i - 1])
                dradi = er[i + 1] - er[i]
                t1 = (coni - conim1) / dradi
                t2 = (conip1 - coni) / dradi
                dlncdr = 0.5 * (t1 + t2)
            edxdr[i] = dlncdr
            dlncdr = 0.0
        else:
            # Question (see the page): EDXDR is stale if hydrogen is not
            # diffused.
            dlncdr = edxdr[i]
        ecod1[i] = fac * (ehqpr[i] * (ap + at) + dlncdr * ah) * espec[j, i]
        ecod2[i] = fac * ad
