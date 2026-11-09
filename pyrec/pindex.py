"""Python version of pindex.f (YREC 7be4133): which shells of a model are
printed."""

import numpy as np

from pyrec import nmlpar_mod

# IEND is kept between calls (SAVE). Suspected bug (see the page): the
# LSHELL branch can read it before setting it.
iend = 0


def pindex(jxbeg: int, jxend: int, lshell: bool, m: int, id_: np.ndarray, idm: int) -> int:
    """List the shells to print in id_[1..idm]: 1, every NPRTPT-th shell
    and m; with lshell, also every shell of the hydrogen-burning shell
    jxbeg..jxend.

    Translates: pindex.f#PINDEX (7be4133)

    Arguments: jxbeg, jxend, the burning shell's edges; lshell; m, the
    number of shells; id_, set in place (ID in the Fortran); idm (set).
    Returns idm, the length of the list.

    """

    global iend
    nprtpt = nmlpar_mod.nprtpt

    # The first shell (Fortran lines 44–45).
    id_[1] = 1
    idm = 2

    # Every zone of the H-burning shell (Fortran lines 46–72).
    if lshell:

        # Every nprtpt-th shell below it (Fortran lines 48–55).
        if nprtpt < jxbeg:
            ibeg = max(2, nprtpt)
            iend = (jxbeg // nprtpt) * nprtpt
            for j in range(ibeg, iend + 1, nprtpt):
                id_[idm] = j
                idm = idm + 1

        # The shell itself (Fortran lines 56–64). Suspected bug (see the
        # page): IEND may be left from an earlier call.
        if iend == jxbeg:
            ixbeg = jxbeg + 1
        else:
            ixbeg = jxbeg
        for j in range(ixbeg, jxend + 1):
            id_[idm] = j
            idm = idm + 1

        # Every nprtpt-th shell above it (Fortran lines 65–72).
        if nprtpt < m:
            ibeg = (jxend // nprtpt + 1) * nprtpt
            iend = (m // nprtpt) * nprtpt
            # Suspected bug (see the page): m may be listed twice here.
            for j in range(ibeg, iend + 1, nprtpt):
                id_[idm] = j
                idm = idm + 1

    # The general case: every nprtpt-th shell (Fortran lines 73–82).
    elif nprtpt < m:
        ibeg = max(2, nprtpt)
        iend = (m // nprtpt) * nprtpt
        if iend == m:
            iend = iend - nprtpt
        for j in range(ibeg, iend + 1, nprtpt):
            id_[idm] = j
            idm = idm + 1

    # The last shell (Fortran line 83).
    id_[idm] = m
    return idm
