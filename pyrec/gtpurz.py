"""Python version of gtpurz.f (YREC development, 7be4133): the opacity of
pure metals (Z = 1) from the Livermore table."""

import numpy as np

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import zlaol_mod
from pyrec import zslaol_mod
from pyrec.cspline import cspline
from pyrec.locate import locate
from pyrec.splint import splint
from pyrec.xrng4 import xrng4

# GTPURZ has SAVE: its work arrays are module-level.
top = np.zeros(105)
tr = np.zeros(105)
tdor = np.zeros(105)
dodr = np.zeros(53)
top1 = np.zeros(53)
tt = np.zeros(53)
tdot = np.zeros(53)


def gtpurz(dl: float, tl: float, o: float, ol: float, qod: float, qot: float) -> tuple[float, float, float, float]:
    """Opacity of a pure-metal mixture, with its derivatives.

    A spline in log rho at each of four temperatures (with linear
    extrapolation within the margin ln TOLLAOL), then a natural spline in
    log T; chord slopes for the derivatives; log kappa capped at 35.

    Translates: gtpurz.f#GTPURZ (7be4133)

    Parameters
    ----------
    dl, tl : float
        log10 density and log10 temperature.
    o, ol, qod, qot : float
        Set here (passed in with any value).

    Returns
    -------
    o, ol, qod, qot
        The opacity, its log10, dlog(kappa)/dlog(rho), dlog(kappa)/dlog(T).

    """

    zl = zlaol_mod
    zs = zslaol_mod

    # The margin for extrapolation in density (Fortran lines 66–68).
    # Question (see the page): a natural log, as in GTLAOL.
    tollog = np.log(nmlpar_mod.tollaol)
    rtl = tl
    rdl = dl

    # The four temperatures around T (Fortran lines 69–72).
    jt = locate(zl.zot, zl.numtz, rtl, 0)
    numts = 0
    jt1, jt2 = xrng4(jt, zl.numtz, 0, 0)

    # A spline in density at each of the four temperatures (Fortran lines 73–113).
    for it in range(jt1, jt2 + 1):
        nums = zs.numrsz[it]
        if nums >= 4:
            # This temperature's splines (Fortran lines 76–80).
            for ir in range(1, nums + 1):
                top[ir] = zs.zslaoll[ir, it]
                tr[ir] = zs.zsrhol[ir, it]
                tdor[ir] = zs.zsdorl[ir, it]

            # Inside the table: the spline (Fortran lines 81–86).
            if tr[1] < rdl < tr[nums]:
                rol, klo, khi = splint(tr, top, nums, tdor, rdl, 0.0, 0, 0)
                numts = numts + 1
                top1[numts] = rol
                tt[numts] = zl.zot[it]
                dodr[numts] = (top[khi] - top[klo]) / (tr[khi] - tr[klo])
            elif tr[1] - tollog < rdl <= tr[1]:
                # Linear extrapolation below the first density (Fortran lines 87–94).
                slope = (top[2] - top[1]) / (tr[2] - tr[1])
                rol = top[1] + slope * (rdl - tr[1])
                numts = numts + 1
                top1[numts] = rol
                tt[numts] = zl.zot[it]
                dodr[numts] = slope
            elif tr[nums] <= rdl < tr[nums] + tollog:
                # Linear extrapolation above the last density (Fortran lines 95–104).
                slope = (top[nums - 1] - top[nums]) / (tr[nums - 1] - tr[nums])
                rol = top[nums] + slope * (rdl - tr[nums])
                numts = numts + 1
                top1[numts] = rol
                tt[numts] = zl.zot[it]
                dodr[numts] = slope
        else:
            # Fewer than four entries at this temperature: stop (Fortran lines 106–111).
            iounits_mod.write(iounits_mod.ishort,
                              f" OUTSIDE Z OPACITY TABLE, IN DENSITY.  "
                              f"LOG(RHO)={dl:12.3E} LOG(T)={tl:12.3E}")  # format 120
            raise SystemExit(1)    # STOP 1, the run-time error code

    # A natural spline in temperature, or a stop (Fortran lines 114–126).
    if numts >= 4:
        cspline(tt, top1, numts, 1.0e30, 1.0e30, tdot)
        rol, klo, khi = splint(tt, top1, numts, tdot, rtl, 0.0, 0, 0)
        slope = (dodr[khi] - dodr[klo]) / (tt[khi] - tt[klo])
        qod = dodr[klo] + slope * (rtl - tt[klo])
        qot = (top1[khi] - top1[klo]) / (tt[khi] - tt[klo])
    else:
        iounits_mod.write(iounits_mod.ishort,
                          f" OUTSIDE Z OPACITY TABLE, IN TEMPERATURE.  "
                          f"LOG(RHO)={dl:12.3E} LOG(T)={tl:12.3E}")      # format 121
        raise SystemExit(1)    # STOP 1, the run-time error code

    # The opacity, its log capped at 35 (Fortran lines 127–133).
    if rol > 35:
        o = 1.0e35
        ol = 35.0
    else:
        o = 10.0 ** rol
        ol = rol
    return o, ol, qod, qot
