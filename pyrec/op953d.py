"""Python version of op953d.f (YREC 7be4133): the OPAL95 opacity from the
table at the model's Z."""

import numpy as np

from pyrec import llot95_mod
from pyrec import op95interp_mod

# The Fortran's local arrays (SAVE).
c = np.zeros(5)
qc = np.zeros(5)
cc = np.zeros(5)
cct = np.zeros(5)
ccd = np.zeros(5)


def op953d(o: float, ol: float, qod: float, qot: float) -> tuple[float, float, float, float]:
    """Get the OPAL95 opacity from the table pre-interpolated to the
    model's Z (MHP 7/98): as OP952D (log R, then log T) for each of the four
    nearest X tables, then four-point interpolation in X, with the indices
    and weights GETOPAL95 left in op95interp_mod.

    Translates: op953d.f#OP953D (7be4133)

    Arguments: o, ol, qod, qot, passed in with any value.

    Returns (o, ol, qod, qot): kappa, log kappa, dlog kappa/dlog rho and
    dlog kappa/dlog T (at constant rho).

    """

    im = op95interp_mod
    cappa = llot95_mod.cappa

    # Find the opacity at each of the 4 nearby values of X (Fortran lines 39–75).
    for j in range(1, 5):
        ix = im.jx[1, j]
        # For each X, log kappa and its slope in log R at four values of T
        # (Fortran lines 42–52).
        for i in range(1, 5):
            it = im.jt + i - 1
            id_ = im.jd[i]
            c[i] = (im.bd[i, 1] * cappa[ix, it, id_]
                    + im.bd[i, 2] * cappa[ix, it, id_ + 1]
                    + im.bd[i, 3] * cappa[ix, it, id_ + 2]
                    + im.bd[i, 4] * cappa[ix, it, id_ + 3])
            qc[i] = (im.bdd[i, 1] * cappa[ix, it, id_]
                     + im.bdd[i, 2] * cappa[ix, it, id_ + 1]
                     + im.bdd[i, 3] * cappa[ix, it, id_ + 2]
                     + im.bdd[i, 4] * cappa[ix, it, id_ + 3])

        # Outside the table in rho: extrapolate with the slope at the edge
        # (Fortran lines 54–68).
        if im.lextl:
            deltar = im.rho_t63l - im.rho_t63ll
            for i in range(1, 5):
                c[i] = c[i] + deltar * qc[i]
        elif im.lexth:
            for i in range(1, 5):
                if im.lexh[i]:
                    deltar = im.rho_t63l - im.rlh[i]
                    c[i] = c[i] + deltar * qc[i]

        # Interpolate in T (Fortran lines 70–74).
        cc[j] = (im.bt[1] * c[1] + im.bt[2] * c[2] + im.bt[3] * c[3]
                 + im.bt[4] * c[4])
        cct[j] = (im.btt[1] * c[1] + im.btt[2] * c[2] + im.btt[3] * c[3]
                  + im.btt[4] * c[4])
        ccd[j] = (im.bt[1] * qc[1] + im.bt[2] * qc[2] + im.bt[3] * qc[3]
                  + im.bt[4] * qc[4])

    # Interpolate in X (Fortran lines 77–84).
    bx = im.bx
    ol = bx[1, 1] * cc[1] + bx[1, 2] * cc[2] + bx[1, 3] * cc[3] + bx[1, 4] * cc[4]
    qot = (bx[1, 1] * cct[1] + bx[1, 2] * cct[2] + bx[1, 3] * cct[3]
           + bx[1, 4] * cct[4])
    qod = (bx[1, 1] * ccd[1] + bx[1, 2] * ccd[2] + bx[1, 3] * ccd[3]
           + bx[1, 4] * ccd[4])

    # From the derivative at fixed R to the derivative at fixed rho
    # (Fortran lines 86–87).
    qot = qot - 3.0 * qod
    o = 10.0**ol
    return o, ol, qod, qot
