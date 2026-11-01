"""Python version of tpgrad.f (YREC 09c3b9b): the Schwarzschild test, the
temperature gradient and the convective velocity at a point."""

import math

from pyrec import constants_mod
from pyrec import dpmix_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import pualpha_mod

VTOL = 1.0e-10  # DATA VTOL: the cubic's tolerance


def tpgrad(tl, t, pl, p, d, rl, sl, b, o, qdt, qdp, qot, qod, qcp, dele,
           delr, dela, qdtt, qdtp, qat, qap, qact, qacp, qacr, qcpt, qcpp,
           vel, lderiv, lconv, fpl, ftl, teffl, cmixl):
    """Decide whether a point is convective and give its temperature
    gradient: the radiative gradient, the Schwarzschild test against the
    adiabatic one, and in convective regions the mixing-length cubic for
    the actual gradient and the convective velocity; with lderiv, the
    gradient's logarithmic derivatives for the Henyey solution.

    Translates: tpgrad.f#TPGRAD (09c3b9b)

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
    pualpha_mod.alfmlt = 0.0
    pualpha_mod.phmlt = 0.0
    pualpha_mod.cmxmlt = 0.0
    # Question (see the page): the pulsation output's MLT values stay zero.
    delr = o * b * math.exp(c.cln * (pl - sl - 4.0 * tl + c.clsunl - c.cgl
                                     + c.cdelrl)) * ftl / fpl
    deldel = delr - dela
    if deldel <= 1.0e-6:
        # Radiative (adiabatic inside an adiabatic overshoot region).
        lconv = False
        dele = delr
        vel = 0.0
        if (nml.ladov and dpmix_mod.iovim >= dpmix_mod.iov1
                and dpmix_mod.iovim <= dpmix_mod.iov2
                and dpmix_mod.iovim != -1):
            dele = dela
        return dele, delr, qact, qacp, qacr, vel, lconv  # label 200
    # Convective.
    lconv = True
    dele = dela
    # Cores: no mixing-length theory.
    if tl > nml.tgcut:
        vel = 1.0e-11
        if lderiv:
            qact = qat
            qacp = qap
            qacr = 0.0
        return dele, delr, qact, qacp, qacr, vel, lconv
    # Spots reduce the radiative flux in convective regions (G Somers 9/14).
    # Question (see the page): delr itself, its derivatives and the
    # radiative/convective test are left uncorrected.
    if nml.spotf != 0.00:
        if nml.lsdepth:
            ateffl = teffl - 0.25 * math.log10(
                nml.spotf * nml.spotx**4.0 + 1.0 - nml.spotf)
            deepx = 1.0 - (1.0 - nml.spotx) * (10.0**ateffl) / (10.0**tl)
        else:
            deepx = nml.spotx
        deldel = delr / (nml.spotf * deepx**4.0 + 1.0 - nml.spotf) - dela
    g = math.exp(c.cln * (c.cgl + sl - rl - rl))
    presht = p / (d * g)
    phi = cmixl * d * o * presht
    phi2 = phi * phi
    phiphi = 1.0 / (1.0 + c.cmixl2 * phi2)
    # Solve the cubic a3 v^3 + v^2 + a1 v - 1 = 0.
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
    # label 30
    ddel = deldel * v * (v + a1)
    dele = dela + ddel
    # The convective velocity.
    test = g * (-qdt) * presht * deldel
    if test > 0.0:
        vel = v * math.sqrt(g * (-qdt) * (cmixl**2) * presht * deldel / 8.0)
    else:
        # Question (see the page): the Fortran's 1.0D0-11 is 1 - 11 = -10,
        # not 1e-11.
        vel = 1.0 - 11
    # Question (see the page): delpm and rrr are not used.
    delpm = dele - v * v * deldel
    rrr = 10.0**rl
    if lderiv:
        # Derivatives of the convective gradient.
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
