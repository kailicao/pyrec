"""Python version of mu.f (YREC 09c3b9b): the ions' and electrons'
mean-molecular-weight terms for the equation of state."""

from pyrec import comp_mod
from pyrec import constants_mod

# DATA ATOMWT: 1/A for H, He4 and C12 (excess metals), and 2/4.004
# (index 0 unused).
# Question (see the page): atomwt[4] is never used.
ATOMWT = (0.0, 0.9921, 0.24975, 0.08322, 0.4995)


def mu(t, p, d, x, z, rmu, amu, emu, beta):
    """Compute the ions' term amu (AMUENV corrected for the change of X and
    Z, excess metals as carbon), the electrons' term emu from the gas
    pressure beta*p, and rmu = CGAS (amu + emu).

    Translates: mu.f#MU (09c3b9b)

    Arguments: t, p, d (rho), x, z and beta; rmu, amu and emu are passed in
    with any value.

    Returns (rmu, amu, emu).
    """
    cgas = constants_mod.cgas
    dfx1 = x - comp_mod.xenv
    dfx12 = z - comp_mod.zenv
    if abs(dfx1) + abs(dfx12) < 1.0e-5:
        # The envelope's abundances.
        amu = comp_mod.amuenv
    else:
        dfx1 = dfx1 * ATOMWT[1]
        dfx12 = dfx12 * ATOMWT[3]
        dfx4 = (comp_mod.xenv + comp_mod.zenv - x - z) * ATOMWT[2]
        # Question (see the page): excess Z (metals) taken as carbon 12.
        amu = comp_mod.amuenv + dfx1 + dfx4 + dfx12
    ee = ((beta * p) / (d * t * cgas * amu)) - 1.0
    emu = ee * amu
    rmu = cgas * (amu + emu)
    return rmu, amu, emu
