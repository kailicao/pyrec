"""Python version of getfc.f (YREC development, 09c3b9b): the factor f_c
relating the chemical diffusion coefficient to that of angular momentum
transport by meridional circulation (Zahn 1992)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import intvar_mod
from pyrec import nmlpar_mod

# Locals kept between calls (SAVE).
alpha = np.zeros(dimensions_mod.json + 1)
v = np.zeros(dimensions_mod.json + 1)


def getfc(hd, hru, hv, imin, imax, omega, vfc):
    """Compute, for the shells imin..imax, Zahn's (1992) ratio f_c of the
    chemical to the angular-momentum diffusion coefficient, capped at 1
    (equations on the page); zero where the circulation velocity is <= 0,
    and everywhere if FC <= 0.

    Translates: getfc.f#GETFC (09c3b9b)

    Arguments: hd (log rho at zone centres), hru (radii at zone centres),
    hv (circulation velocities U at zone edges), imin, imax, omega. Array
    changed in place: vfc (owned by seculr, read by codiff). Nothing is
    returned.
    """
    c = constants_mod
    iv = intvar_mod
    fc = nmlpar_mod.fc
    if fc <= 0.0:
        for i in range(imin, imax + 1):
            vfc[i] = 0.0
        return                                     # go to 9999
    # alpha = d ln(r^2 Omega) / d ln r between neighbouring centres.
    # Question (see the page): no factor 1/2.
    for i in range(imin, imax + 1):
        wm = 0.5 * (omega[i] + omega[i - 1])
        alpha[i] = ((omega[i] * hru[i]**2 - omega[i - 1] * hru[i - 1]**2)
                    / (hru[i] - hru[i - 1]) / (wm * iv.rm[i]))
    # v = (1/(3 rho r)) d(rho r^2 U)/dr, i.e. 2V; U interpolated to the
    # zone centres.
    fx = (hru[imin] - iv.rm[imin]) / (iv.rm[imin + 1] - iv.rm[imin])
    vm = hv[imin] + fx * (hv[imin + 1] - hv[imin])
    dr2vm = math.exp(c.cln * hd[imin]) * hru[imin]**2 * vm
    for i in range(imin + 1, imax):
        fx = (hru[i] - iv.rm[i]) / (iv.rm[i + 1] - iv.rm[i])
        vp = hv[i] + fx * (hv[i + 1] - hv[i])
        dr2vp = math.exp(c.cln * hd[i]) * hru[i]**2 * vp
        rhor = iv.dm[i] * iv.rm[i]
        v[i] = c.cc13 * (dr2vp - dr2vm) / rhor / (hru[i] - hru[i - 1])
        dr2vm = dr2vp
    # The first and last shells, by the product rule.
    rhor = iv.dm[imin] * iv.rm[imin]
    drhor2 = ((math.exp(c.cln * hd[imin]) * hru[imin]**2
               - math.exp(c.cln * hd[imin - 1]) * hru[imin - 1]**2)
              / (hru[imin] - hru[imin - 1]))
    dudr = (hv[imin + 1] - hv[imin]) / (iv.rm[imin + 1] - iv.rm[imin])
    v[imin] = c.cc13 * (hv[imin] * drhor2 / rhor + iv.rm[imin] * dudr)
    rhor = iv.dm[imax] * iv.rm[imax]
    drhor2 = ((math.exp(c.cln * hd[imax]) * hru[imax]**2
               - math.exp(c.cln * hd[imax - 1]) * hru[imax - 1]**2)
              / (hru[imax] - hru[imax - 1]))
    dudr = (hv[imax] - hv[imax - 1]) / (iv.rm[imax] - iv.rm[imax - 1])
    v[imax] = c.cc13 * (hv[imax] * drhor2 / rhor + iv.rm[imax] * dudr)
    # f_c = FC / max(|v - alpha U| / U, FC).
    for i in range(imin, imax + 1):
        if hv[i] <= 0.0:
            # Question (see the page): no mixing where U <= 0.
            vfc[i] = 0.0
        else:
            ftest = max(abs(v[i] - alpha[i] * hv[i]) / hv[i], fc)
            vfc[i] = fc / ftest
