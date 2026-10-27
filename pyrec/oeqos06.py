"""Python version of oeqos06.f (YREC 09c3b9b): the interface to the 2006
OPAL equation of state (LLP, October 17, 2006)."""

import math

from pyrec import constants_mod
from pyrec import eos06misc_mod
from pyrec import iounits_mod
from pyrec.esac06 import esac06
from pyrec.mu import mu
from pyrec.rhoofp06 import rhoofp06

# The table's dimensions (unused here).
MX, MV, NR, NT = 5, 10, 169, 197


def oeqos06(t10, t, p10, p, d10, d, xc, z, beta, betai, beta14, rmu, amu,
            emu, qdt, qdp, qcp, dela):
    """From log T and log P, find the density in the 2006 OPAL tables
    (RHOOFP06, radiation included), the thermodynamic derivatives (ESAC06),
    the radiation pressure fraction and the mean-molecular-weight terms
    (MU). Fails (alt1) outside 0.00187 <= T6 <= 200 or where the table
    fails, so that the caller falls back on another equation of state.

    Translates: oeqos06.f#OEQOS06 (09c3b9b)

    Arguments: t10 = log T, p10 = log P, the composition xc (X) and z; the
    others are passed in with any value.

    Returns (t, p, d10, d, beta, betai, beta14, rmu, amu, emu, qdt, qdp,
    qcp, dela, alt1): T and P unlogged, log rho and rho, beta = 1 - aT^4/3P,
    1/beta and 1 - beta, the mean-molecular-weight terms, dln rho/dln T at
    constant P, dln rho/dln P at constant T, c_P, nabla_ad, and alt1, True
    when the Fortran takes its alternate return (RETURN 1).
    """
    def out(alt1):
        return (t, p, d10, d, beta, betai, beta14, rmu, amu, emu, qdt, qdp,
                qcp, dela, alt1)

    iorder = 9   # all first and second order data
    irad = 1     # adds radiation corrections
    t = 10.0**t10
    p = 10.0**p10
    t6 = t / 1.0e6
    p12 = p / 1.0e12
    if t6 < 0.001870 or t6 > 200.0:
        return out(True)                                 # RETURN 1
    x = xc
    ztab = z
    r = rhoofp06(x, t6, p12, irad)
    if r <= -998.0:
        return out(True)                                 # RETURN 1
    d = r
    d10 = math.log10(d)
    if esac06(x, t6, r, iorder, irad):
        iounits_mod.write(iounits_mod.ishort, " WARNNING... OPAL TBL FAIL")
        return out(True)                                 # label 999
    # Question (see the page): the pressure check is commented out.
    eos = eos06misc_mod.eos
    qdp = 1.0 / eos[6]                 # EOS(6): dlogP/dlogRho at const T6
    qdt = -eos[7] / eos[6]             # EOS(7): dlogP/dlogT6 at const Rho
    qcp = 1.0e6 * eos[5] * eos[8] / eos[6]   # EOS(5): dE/dT6, EOS(8): Gamma1
    dela = 1.0 / eos[9]                # EOS(9): Gamma2/(Gamma2-1)
    beta14 = (2.521971383e-3 * t6 * t6) * (t6 * t6 / p12)
    beta = 1.0 - beta14
    betai = 1.0 / beta
    rmu = p * beta / (d * t)
    rmu2, amu, emu = mu(t, p, d, xc, z, 0.0, amu, emu, beta)
    if emu <= 0.0:
        emu = 0.0
        amu = rmu / constants_mod.cgas - emu
    return out(False)
