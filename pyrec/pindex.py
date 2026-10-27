"""Python version of pindex.f (YREC 09c3b9b): which shells of a model are
printed."""

from pyrec import nmlpar_mod

# IEND is kept between calls (SAVE): question 1 reads it before setting it.
iend = 0


def pindex(jxbeg, jxend, lshell, m, id_, idm):
    """List the shells to print in id_[1..idm]: 1, every NPRTPT-th shell
    and m; with lshell, also every shell of the hydrogen-burning shell
    jxbeg..jxend.

    Translates: pindex.f#PINDEX (09c3b9b)

    Arguments: jxbeg, jxend, the burning shell's edges; lshell; m, the
    number of shells; id_, set in place (ID in the Fortran); idm (set).
    Returns idm, the length of the list.
    """
    global iend
    nprtpt = nmlpar_mod.nprtpt
    id_[1] = 1
    idm = 2
    if lshell:
        if nprtpt < jxbeg:
            ibeg = max(2, nprtpt)
            iend = (jxbeg // nprtpt) * nprtpt
            for j in range(ibeg, iend + 1, nprtpt):
                id_[idm] = j
                idm = idm + 1
        # Question (see the page): IEND may be left from an earlier call.
        if iend == jxbeg:
            ixbeg = jxbeg + 1
        else:
            ixbeg = jxbeg
        for j in range(ixbeg, jxend + 1):
            id_[idm] = j
            idm = idm + 1
        if nprtpt < m:
            ibeg = (jxend // nprtpt + 1) * nprtpt
            iend = (m // nprtpt) * nprtpt
            # Question (see the page): m may be listed twice here.
            for j in range(ibeg, iend + 1, nprtpt):
                id_[idm] = j
                idm = idm + 1
    elif nprtpt < m:
        ibeg = max(2, nprtpt)
        iend = (m // nprtpt) * nprtpt
        if iend == m:
            iend = iend - nprtpt
        for j in range(ibeg, iend + 1, nprtpt):
            id_[idm] = j
            idm = idm + 1
    id_[idm] = m
    return idm
