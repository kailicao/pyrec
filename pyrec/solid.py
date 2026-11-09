"""Python version of solid.f (YREC 7be4133): solid-body rotation of a range
of shells at fixed total angular momentum."""

from pyrec import disk_mod
import numpy as np
from pyrec import errmom_mod
from pyrec import nmlpar_mod
from pyrec.momi import momi
from pyrec.shape import shape


def solid(hd: np.ndarray, hjm: np.ndarray, hr: np.ndarray, hs: np.ndarray, hs2: np.ndarray, jstart: int, jend: int,
          eta2: np.ndarray, hi: np.ndarray, omega: np.ndarray, qiw: np.ndarray, r0: np.ndarray, m: int) -> None:
    """Make shells jstart..jend rotate as a solid body with their total
    angular momentum J kept: Newton iteration on Omega I(Omega) = J, or,
    under disk locking (the surface zone, LDISK and SAGE <= TDISK), the
    fixed rate PDISK. Then reset each shell's specific angular momentum.

    Translates: solid.f#SOLID (7be4133)

    Arguments: hd, hr, hs (log rho, log r, log m); hjm (specific angular
    momenta, reset in place); hs2 (shell masses); jstart, jend; eta2, hi,
    omega, qiw, r0, set in place (through SHAPE and MOMI); m (the number of
    shells). None is returned.

    """

    nml = nmlpar_mod

    # The zone's total angular momentum and the mean of the previous rates
    # (Fortran lines 49–57).
    czj = 0.0
    wsum = 0.0
    icount = 0
    for j in range(jstart, jend + 1):
        czj = czj + hs2[j] * hjm[j]
        wsum = wsum + omega[j]

    # Disk locking in the surface convection zone, and the first guess
    # (Fortran lines 58–67).
    lock = False
    if (jend == m and jstart != jend and nml.ldisk
            and disk_mod.sage <= nml.tdisk):
        lock = True
    if not lock:
        wguess = wsum / float(jend - jstart + 1)
    else:
        # Question (see the page): PDISK is a rate.
        wguess = nml.pdisk

    # Newton's iteration on Omega I(Omega) = J (Fortran lines 68–96).
    while True:
        for j in range(jstart, jend + 1):  # label 20
            omega[j] = wguess
        shape(hd, hr, hs, jstart, jend, omega, eta2, r0)
        momi(eta2, hr, hs, hs2, jstart, jend, omega, r0, hi, qiw)
        czi = 0.0
        czqiw = 0.0
        for j in range(jstart, jend + 1):
            czi = czi + hi[j]
            czqiw = czqiw + qiw[j]
        czjnew = wguess * czi
        delj = czj - czjnew
        if lock:
            break
        # Suspected bug (see the page): divides by CZJ, zero for a zone
        # without rotation.
        # Question (see the page): no warning after 20 steps.
        if abs(delj / czj) > errmom_mod.toleri and icount < 20:
            icount = icount + 1
            # Suspected bug (see the page): QIW is A dI/dA, so Newton's
            # denominator would be czi + 2 czqiw.
            delw = delj / (czi + wguess * czqiw)
            wguess = wguess + delw
            continue
        break

    # Each shell's specific angular momentum at the final rate (Fortran
    # lines 96–99).
    for j in range(jstart, jend + 1):  # label 45
        hjm[j] = hi[j] * omega[j] / hs2[j]
