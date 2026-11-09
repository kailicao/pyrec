"""Python version of tpgrad.f (YREC 7be4133): the Schwarzschild test, the
temperature gradient and the convective velocity at a point."""

import math

from pyrec import constants_mod
from pyrec import dpmix_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import pualpha_mod

VTOL = 1.0e-10  # DATA VTOL: the cubic's tolerance


def tpgrad(tl: float, t: float, pl: float, p: float, d: float, rl: float, sl: float, b: float, o: float, qdt: float,
           qdp: float, qot: float, qod: float, qcp: float, dele: float, delr: float, dela: float, qdtt: float,
           qdtp: float, qat: float, qap: float, qact: float, qacp: float, qacr: float, qcpt: float, qcpp: float,
           vel: float, lderiv: bool, lconv: bool, fpl: float, ftl: float, teffl: float, cmixl: float) -> tuple[float,
           float, float, float, float, float, bool]:
    """Decide whether a point is convective and give its temperature
    gradient: the radiative gradient, the Schwarzschild test against the
    adiabatic one, and in convective regions the mixing-length cubic for
    the actual gradient and the convective velocity; with lderiv, the
    gradient's logarithmic derivatives for the Henyey solution.

    Translates: tpgrad.f#TPGRAD (7be4133)

    Arguments: log T, T, log P, P, rho (d), log r, log m (sl), L/Lsun (b),
    the opacity o and the equation-of-state quantities (qdt = -delta, qdp,
    qot, qod, qcp = c_P, dela = nabla_ad, and their derivatives qdtt, qdtp,
    qat, qap, qcpt, qcpp), lderiv, the rotation factors fpl, ftl, log Teff
    and cmixl (the mixing length in pressure scale heights). dele (the
    Fortran's DEL), delr, qact, qacp, qacr, vel and lconv are outputs; the
    derivatives keep their values where the Fortran does not set them.

    Returns (dele, delr, qact, qacp, qacr, vel, lconv).

    """

    c = constants_mod
    nml = nmlpar_mod

    # The pulsation output's MLT values, always zero (Fortran lines 77–79).
    pualpha_mod.alfmlt = 0.0
    pualpha_mod.phmlt = 0.0
    pualpha_mod.cmxmlt = 0.0

    # The radiative gradient and the Schwarzschild test (Fortran lines
    # 80–83).
    # Question (see the page): the margin 1e-6.
    delr = o * b * math.exp(c.cln * (pl - sl - 4.0 * tl + c.clsunl - c.cgl
                                     + c.cdelrl)) * ftl / fpl
    deldel = delr - dela
    if deldel <= 1.0e-6:
        # Radiative (adiabatic inside an adiabatic overshoot region)
        # (Fortran lines 84–92).
        lconv = False
        dele = delr
        vel = 0.0
        if (nml.ladov and dpmix_mod.iovim >= dpmix_mod.iov1
                and dpmix_mod.iovim <= dpmix_mod.iov2
                and dpmix_mod.iovim != -1):
            dele = dela
        return dele, delr, qact, qacp, qacr, vel, lconv  # label 200

    # Convective (Fortran lines 94–96).
    lconv = True
    dele = dela

    # Cores: no mixing-length theory (Fortran lines 97–107).
    # Question (see the page): skipping MLT in cores.
    if tl > nml.tgcut:
        vel = 1.0e-11
        if lderiv:
            qact = qat
            qacp = qap
            qacr = 0.0
        return dele, delr, qact, qacp, qacr, vel, lconv

    # Spots reduce the radiative flux in convective regions (G Somers 9/14)
    # (Fortran lines 108–136).
    # Question (see the page): delr itself and its derivatives are left
    # uncorrected.
    # Suspected bug (see the page): the radiative/convective test is made
    # without the correction, and a deldel the correction makes negative
    # stops the run below.
    if nml.spotf != 0.00:
        if nml.lsdepth:
            ateffl = teffl - 0.25 * math.log10(
                nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
            deepx = 1.0 - (1.0 - nml.spotx) * (10.0**ateffl) / (10.0**tl)
        else:
            deepx = nml.spotx
        deldel = delr / (nml.spotf * deepx**4.0 + 1.0 - nml.spotf) - dela

    # Gravity, the scale height and the element's optical thickness (Fortran
    # lines 137–141).
    # Question (see the page): which variant of the theory.
    g = math.exp(c.cln * (c.cgl + sl - rl - rl))
    presht = p / (d * g)
    phi = cmixl * d * o * presht
    phi2 = phi * phi
    phiphi = 1.0 / (1.0 + c.cmixl2 * phi2)

    # Solve the cubic a3 v^3 + v^2 + a1 v - 1 = 0: a1, or a stop on a
    # negative square root (Fortran lines 142–159).
    test = deldel * g * (-qdt) / presht
    if test > 0.0:
        a1 = c.cmixl3 * phiphi * t**3 * o / (qcp * math.sqrt(test))
    else:
        print(f"***WARNING***\nTPGRAD TRIED TO TAKE SQUARE ROOT"
              f" OF NEGATIVE NUMBER - RUN ENDED\n"
              f" DELDEL {deldel:12.3E} QDT {qdt:12.3E} DISC. {test:12.3E}")
        print("".join(f"{v:12.3E}" for v in
                      (o, b, pl, sl, tl, ftl, fpl, delr, dela, rl)))
        raise SystemExit(912)

    # The Newton iteration (Fortran lines 160–171).
    v = 1.0 / a1
    a3 = 0.75 * phi2 * phiphi / a1
    a3p = 3.0 * a3
    if a3 > 1.0e+3:
        v = a3**(-c.cc13)
    converged = False
    for it in range(1, 26):  # loop 10
        v = min(v, 1.0)
        vp = a1 + v * (2.0 + v * a3p)
        vd = (-1.0 + v * (a1 + v * (1.0 + v * a3))) / vp
        vd = vd * (1.0 + vd * (1.0 + v * a3p) / vp)
        v = v - vd
        if abs(vd) < VTOL:
            converged = True
            break  # to label 30

    # No convergence: assume adiabatic convection (Fortran lines 172–184).
    if not converged:
        iounits_mod.write(iounits_mod.ishort,
                          f" -----CUBIC NON-CONVERGENCE(PL,TL,CAPPA,CP,QDT)="
                          f"{pl:10.6f}{tl:10.6f}\n"
                          f"{o:12.6E}{qcp:12.6E}{qdt:12.6E}")  # format 20
        # Assume adiabatic convection.
        vel = 1.0e-11
        if lderiv:
            qact = qat
            qacp = qap
            qacr = 0.0
        return dele, delr, qact, qacp, qacr, vel, lconv

    # label 30: the gradient (Fortran lines 185–186).
    ddel = deldel * v * (v + a1)
    dele = dela + ddel

    # The convective velocity; its radicand is positive, since a
    # non-positive test above stopped (MHP 10/30/26) (Fortran lines
    # 187–194).
    vel = v * math.sqrt(g * (-qdt) * (cmixl**2) * presht * deldel / 8.0)
    delpm = dele - v * v * deldel  # noqa: F841  (unused in the Fortran too)
    rrr = 10.0**rl  # noqa: F841  (unused in the Fortran too)

    # Derivatives of the convective gradient (Fortran lines 195–223).
    if lderiv:
        qdelat = qat * dela
        qdelap = qap * dela
        tempot = qod * qdt + qot
        tempop = qod * qdp
        qddelt = ((-4.0 + tempot) * delr - qdelat) / deldel
        qddelp = ((+1.0 + tempop) * delr - qdelap) / deldel
        temp1 = 2.0 * c.cmixl2 * phi2 * phiphi
        qa1t = (tempot * (1.0 - temp1) - qcpt
                - 0.5 * (qdt + qddelt + qdtt / qdt) + 3.0)
        qa1p = (tempop * (1.0 - temp1) - temp1 - qcpp
                - 0.5 * (qdp + qddelp + qdtp / qdt) + 0.5)
        qa1r = -2.0 * temp1 + 2.0
        qa3t = 2.0 * phiphi * tempot - qa1t
        qa3p = 2.0 * phiphi * (1.0 + tempop) - qa1p
        qa3r = 4.0 * phiphi - qa1r
        temp1 = a1 * v
        temp2 = a3 * v * v * v
        temp3 = 1.0 / vp
        qvt = -(temp1 * qa1t + temp2 * qa3t) * temp3
        qvp = -(temp1 * qa1p + temp2 * qa3p) * temp3
        qvr = -(temp1 * qa1r + temp2 * qa3r) * temp3
        temp1 = deldel * (a1 + v + v)
        temp2 = deldel * a1 * v
        deli = 1.0 / dele
        qact = (qdelat + qddelt * ddel + temp1 * qvt + temp2 * qa1t) * deli
        qacp = (qdelap + qddelp * ddel + temp1 * qvp + temp2 * qa1p) * deli
        qacr = (temp1 * qvr + temp2 * qa1r) * deli

    return dele, delr, qact, qacp, qacr, vel, lconv  # label 200
