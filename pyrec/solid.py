"""Python version of solid.f (YREC 09c3b9b): solid-body rotation of a range
of shells at fixed total angular momentum."""

from pyrec import disk_mod
from pyrec import errmom_mod
from pyrec import nmlpar_mod
from pyrec.momi import momi
from pyrec.shape import shape


def solid(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega, qiw, r0, m):
    """Make shells jstart..jend rotate as a solid body with their total
    angular momentum J kept: Newton iteration on Omega I(Omega) = J, or,
    under disk locking (the surface zone, LDISK and SAGE <= TDISK), the
    fixed rate PDISK. Then reset each shell's specific angular momentum.

    Translates: solid.f#SOLID (09c3b9b)

    Arguments: hd, hr, hs (log rho, log r, log m); hjm (specific angular
    momenta, reset in place); hs2 (shell masses); jstart, jend; eta2, hi,
    omega, qiw, r0, set in place (through SHAPE and MOMI); m (the number of
    shells). None is returned.
    """
    nml = nmlpar_mod
    czj = 0.0
    wsum = 0.0
    icount = 0
    for j in range(jstart, jend + 1):
        czj = czj + hs2[j] * hjm[j]
        wsum = wsum + omega[j]
    lock = False
    if (jend == m and jstart != jend and nml.ldisk
            and disk_mod.sage <= nml.tdisk):
        lock = True
    if not lock:
        wguess = wsum / float(jend - jstart + 1)
    else:
        # Question (see the page): PDISK is a rate.
        wguess = nml.pdisk
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
        # Question (see the page): TOLERI is never set; CZJ may be zero;
        # no warning after 20 steps.
        if abs(delj / czj) > errmom_mod.toleri and icount < 20:
            icount = icount + 1
            delw = delj / (czi + wguess * czqiw)
            wguess = wguess + delw
            continue
        break
    for j in range(jstart, jend + 1):  # label 45
        hjm[j] = hi[j] * omega[j] / hs2[j]
