"""Python version of rhoofp06.f (YREC 09c3b9b): the density for a given X,
T6 and pressure, from the OPAL 2006 equation of state."""

from pyrec import aeos06_mod
from pyrec import eos06misc_mod
from pyrec.esac06 import esac06

MX = 5
MV = 10
NR = 169
NT = 197
SIGMACC = 1.8914785e-3    # DATA: a/3 in the table's units (Mb per T6**4)


def rhoofp06(x, t6, p, irad):
    """Invert the OPAL 2006 EOS: find the density whose gas pressure
    matches p (less the radiation pressure when irad is 1), by bisection
    in the table's X and T6 grids, then up to 11 steps of regula falsi on
    esac06 called without radiation.

    Translates: rhoofp06.f#RHOOFP06 (09c3b9b)

    Arguments: x (X), t6 (T in 1e6 K), p (pressure, Mb), irad (1: p
    includes radiation). None is set.

    Returns the density, or -999.0 on failure (out of the table, no
    convergence, or esac06's alternate return).
    """
    a = aeos06_mod
    eos = eos06misc_mod.eos
    nra = eos06misc_mod.nra
    rat = SIGMACC
    pr = 0.0
    if irad == 1:
        pr = 4.0 / 3.0 * rat * t6**4    # Mb
    pnr = p - pr

    # Bisection for the composition index (label 8).
    ilo = 2
    ihi = MX
    while ihi - ilo > 1:
        imd = (ihi + ilo) // 2
        if x <= a.xa[imd] + 1.0e-7:
            ihi = imd
        else:
            ilo = imd
    mlo = ilo

    # Bisection for the temperature index (label 11), from the top.
    ilo = NT
    ihi = 2
    while ilo - ihi > 1:
        imd = (ihi + ilo) // 2
        if t6 == a.t6list[1, imd]:
            ilo = imd
            break                                   # go to 14
        if t6 <= a.t6list[1, imd]:
            ihi = imd
        else:
            ilo = imd
    klo = ilo                                       # label 14

    pmax = a.xz[mlo, 1, klo, nra[klo]] * t6 * a.rho[nra[klo]]
    pmin = a.xz[mlo, 1, klo, 1] * t6 * a.rho[1]
    if pnr > 1.25 * pmax or pnr < pmin:
        return -999.0                               # go to 999

    rhog1 = a.rho[nra[klo]] * pnr / pmax
    if esac06(x, t6, rhog1, 1, 0):
        return -999.0
    p1 = eos[1]
    if p1 > pnr:
        p2 = p1
        rhog2 = rhog1
        rhog1 = 0.2 * rhog1
        if rhog1 < 1.0e-14:
            rhog1 = 1.0e-14
        if esac06(x, t6, rhog1, 1, 0):
            return -999.0
        p1 = eos[1]
    else:
        rhog2 = 5.0 * rhog1
        if rhog2 > a.rho[nra[klo]]:
            rhog2 = a.rho[nra[klo]]
        if esac06(x, t6, rhog2, 1, 0):
            return -999.0
        p2 = eos[1]

    # Regula falsi (label 1), at most 11 steps.
    icount = 0
    while True:
        icount = icount + 1
        rhog3 = rhog1 + (rhog2 - rhog1) * (pnr - p1) / (p2 - p1)
        if esac06(x, t6, rhog3, 1, 0):
            return -999.0
        p3 = eos[1]
        # Question (see the page): the tolerance 1e-5 (MR 2025-10-10).
        if abs((p3 - pnr) / pnr) < 1.0e-5:
            return rhog3
        if p3 > pnr:
            rhog2 = rhog3
            p2 = p3
        else:
            rhog1 = rhog3
            p1 = p3
        if icount >= 11:
            return -999.0                           # no convergence
