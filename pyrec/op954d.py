"""Python version of op954d.f (YREC development, 09c3b9b): the OPAL 1995
opacity from the full set of tables, by interpolation in Z, X, T and R."""

import numpy as np

from pyrec import llot95a_mod
from pyrec import op95interp_mod

# OP954D has SAVE: its work arrays are kept between calls.
c = np.zeros(5)
qc = np.zeros(5)
cc = np.zeros(5)
cct = np.zeros(5)
ccd = np.zeros(5)
cz = np.zeros(5)
czt = np.zeros(5)
czd = np.zeros(5)


def op954d(o, ol, qod, qot):
    """The OPAL95 opacity at the current point from the full set of tables
    (MHP 7/98): the log R and log T interpolation of OP952D for each of the
    4 x 4 nearest (Z, X) tables, then four-point interpolation in X and in
    Z, with the indices and weights GETOPAL95 left in OP95INTERP_MOD.

    Translates: op954d.f#OP954D (09c3b9b)

    Arguments: o, ol, qod, qot (set only; any value).

    Returns (o, ol, qod, qot): the opacity kappa, log10 kappa, and
    dlog kappa/dlog rho and dlog kappa/dlog T at constant rho.
    """
    tab = llot95a_mod
    ip = op95interp_mod
    cappa2 = tab.cappa2
    bd, bdd, bt, btt = ip.bd, ip.bdd, ip.bt, ip.btt
    # Find the opacity at each of the 4 nearby values of Z.
    for k in range(1, 5):
        # Question (see the page): the slots NZ(IZ)+10 at Z = 0.06, 0.08.
        iz = tab.nz[ip.jz + k - 1]
        # Find the opacity at each of the 4 nearby values of X.
        for j in range(1, 5):
            ixz = ip.jx[k, j] + iz
            # For each X, kappa at 4 values of T.
            for i in range(1, 5):
                it = ip.jt + i - 1
                # For each T, kappa at 4 values of R.
                id_ = ip.jd[i]
                # log kappa at the desired rho for each of the 4 T.
                c[i] = (bd[i, 1] * cappa2[ixz, it, id_]
                        + bd[i, 2] * cappa2[ixz, it, id_ + 1]
                        + bd[i, 3] * cappa2[ixz, it, id_ + 2]
                        + bd[i, 4] * cappa2[ixz, it, id_ + 3])
                # dlog kappa/dlog R for each of the 4 T.
                qc[i] = (bdd[i, 1] * cappa2[ixz, it, id_]
                         + bdd[i, 2] * cappa2[ixz, it, id_ + 1]
                         + bdd[i, 3] * cappa2[ixz, it, id_ + 2]
                         + bdd[i, 4] * cappa2[ixz, it, id_ + 3])
            # Rho outside the table: extrapolate linearly with the
            # derivative at the boundary.
            if ip.lextl:
                deltar = ip.rho_t63l - ip.rho_t63ll
                for i in range(1, 5):
                    c[i] = c[i] + deltar * qc[i]
            elif ip.lexth:
                for i in range(1, 5):
                    if ip.lexh[i]:
                        deltar = ip.rho_t63l - ip.rlh[i]
                        c[i] = c[i] + deltar * qc[i]
            # Interpolate in T: log kappa, dlog kappa/dlog T, and
            # dlog kappa/dlog R.
            cc[j] = bt[1] * c[1] + bt[2] * c[2] + bt[3] * c[3] + bt[4] * c[4]
            cct[j] = (btt[1] * c[1] + btt[2] * c[2] + btt[3] * c[3]
                      + btt[4] * c[4])
            ccd[j] = (bt[1] * qc[1] + bt[2] * qc[2] + bt[3] * qc[3]
                      + bt[4] * qc[4])
        # Interpolate in X: log kappa, QOT and QOD.
        bx = ip.bx
        cz[k] = (bx[k, 1] * cc[1] + bx[k, 2] * cc[2] + bx[k, 3] * cc[3]
                 + bx[k, 4] * cc[4])
        czt[k] = (bx[k, 1] * cct[1] + bx[k, 2] * cct[2] + bx[k, 3] * cct[3]
                  + bx[k, 4] * cct[4])
        czd[k] = (bx[k, 1] * ccd[1] + bx[k, 2] * ccd[2] + bx[k, 3] * ccd[3]
                  + bx[k, 4] * ccd[4])
    # Interpolate in Z.
    bz = ip.bz
    ol = bz[1] * cz[1] + bz[2] * cz[2] + bz[3] * cz[3] + bz[4] * cz[4]
    qot = bz[1] * czt[1] + bz[2] * czt[2] + bz[3] * czt[3] + bz[4] * czt[4]
    qod = bz[1] * czd[1] + bz[2] * czd[2] + bz[3] * czd[3] + bz[4] * czd[4]
    # From the derivative at fixed R to the derivative at fixed rho.
    qot = qot - 3.0 * qod
    o = 10.0**ol
    return o, ol, qod, qot
