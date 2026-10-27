"""Python version of op952d.f (YREC 09c3b9b): the OPAL95 opacity from the
table at the surface X and Z."""

import numpy as np

from pyrec import llot95e_mod
from pyrec import op95interp_mod

# The Fortran's local arrays (SAVE).
c = np.zeros(5)
qc = np.zeros(5)


def op952d(o, ol, qod, qot):
    """Get the OPAL95 opacity from the table pre-interpolated to the
    surface composition (MHP 7/98): four-point Lagrange interpolation in
    log R (R = rho/T6**3), with linear extrapolation beyond the table, then
    in log T, with the indices and weights GETOPAL95 left in
    op95interp_mod.

    Translates: op952d.f#OP952D (09c3b9b)

    Arguments: o, ol, qod, qot, passed in with any value.

    Returns (o, ol, qod, qot): kappa, log kappa, dlog kappa/dlog rho and
    dlog kappa/dlog T (at constant rho).
    """
    im = op95interp_mod
    cappae = llot95e_mod.cappae
    # Log kappa at the desired rho, and its slope in log R, at each of the
    # four nearby values of T.
    for i in range(1, 5):
        it = im.jt + i - 1
        id_ = im.jd[i]
        c[i] = (im.bd[i, 1] * cappae[it, id_] + im.bd[i, 2] * cappae[it, id_ + 1]
                + im.bd[i, 3] * cappae[it, id_ + 2]
                + im.bd[i, 4] * cappae[it, id_ + 3])
        qc[i] = (im.bdd[i, 1] * cappae[it, id_]
                 + im.bdd[i, 2] * cappae[it, id_ + 1]
                 + im.bdd[i, 3] * cappae[it, id_ + 2]
                 + im.bdd[i, 4] * cappae[it, id_ + 3])
    # Outside the table in rho: extrapolate with the slope at the edge.
    if im.lextl:
        deltar = im.rho_t63l - im.rho_t63ll
        for i in range(1, 5):
            c[i] = c[i] + deltar * qc[i]
    elif im.lexth:
        for i in range(1, 5):
            if im.lexh[i]:
                deltar = im.rho_t63l - im.rlh[i]
                c[i] = c[i] + deltar * qc[i]
    # Interpolate in T: log kappa, dlog kappa/dlog T and dlog kappa/dlog R.
    ol = im.bt[1] * c[1] + im.bt[2] * c[2] + im.bt[3] * c[3] + im.bt[4] * c[4]
    qot = (im.btt[1] * c[1] + im.btt[2] * c[2] + im.btt[3] * c[3]
           + im.btt[4] * c[4])
    qod = (im.bt[1] * qc[1] + im.bt[2] * qc[2] + im.bt[3] * qc[3]
           + im.bt[4] * qc[4])
    # From the derivative at fixed R to the derivative at fixed rho.
    qot = qot - 3.0 * qod
    o = 10.0**ol
    return o, ol, qod, qot
