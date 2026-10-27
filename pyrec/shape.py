"""Python version of shape.f (YREC 09c3b9b): the distortion of each shell
of a rotating star (Radau's eta_2 and the characteristic radius r_0)."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod


def _deta(dr, rho, rhobar, eta, r):
    """One Runge-Kutta increment of Radau's equation over dr, at radius r."""
    return dr * (6.0 - 6.0 * rho * (eta + 1.0) / rhobar
                 - eta * (eta - 1.0)) / r


def shape(hd, hr, hs, jstart, jend, omega, eta2, r0):
    """For shells jstart..jend: integrate Radau's equation outward for
    eta_2, solve for the characteristic radius r_0 (the shell's radius
    being r_0 (1 - A P_2)) by Newton's method, and iterate the two.

    Translates: shape.f#SHAPE (09c3b9b)

    Arguments: hd, hr, hs (log rho, log r, log m, cgs); jstart, jend;
    omega (angular velocity per shell); eta2 and r0, changed in place
    (when jstart > 1, eta2[jstart-1] and r0[jstart-1] are read as the
    starting values). None is returned.
    """
    c = constants_mod
    nml = nmlpar_mod
    cg = math.exp(c.cln * c.cgl)
    c1 = 0.6
    c2 = 2.0 / 3.5e1
    c3 = 3.0 * c1
    c4 = 4.0 * c2
    jbegin = jstart
    if jstart == 1:
        # The first shell: eta_2 from its leading term near the centre.
        # Question (see the page): the centre.
        rhobar = math.exp(c.cln * (hs[1] - c.c4pi3l - 3.0 * hr[1]))
        rphi = math.exp(c.cln * hr[1])
        rho = math.exp(c.cln * hd[1])
        eta2[1] = 6.0 * (1.0 - rho / rhobar)
        gm = cg * math.exp(c.cln * hs[1])
        rphi3 = rphi ** 3
        r03 = rphi3
        fact = 5.0 * c.cc13 * omega[1] ** 2 / (gm * (2.0 + eta2[1]))
        # Question (see the page): silent non-convergence.
        for j in range(1, nml.itfp2 + 1):
            a = fact * r03
            delr03 = ((rphi3 - r03 * (1.0 + c1 * a ** 2 - c2 * a ** 3))
                      / (1.0 + c3 * a ** 2 - c4 * a ** 3))
            r03 = r03 + delr03
            if abs(delr03 / r03) <= nml.acfpft:
                break
        r0[1] = r03 ** c.cc13  # label 20
        if jend == 1:
            return
        rhobp = rhobar
        rphip = rphi
        rhop = rho
        jbegin = 2
    else:
        rhobp = math.exp(c.cln * (hs[jstart - 1] - c.c4pi3l
                                  - 3.0 * hr[jstart - 1]))
        rphip = math.exp(c.cln * hr[jstart - 1])
        rhop = math.exp(c.cln * hd[jstart - 1])

    for i in range(jbegin, jend + 1):
        rhobar = math.exp(c.cln * (hs[i] - c.c4pi3l - 3.0 * hr[i]))
        rho = math.exp(c.cln * hd[i])
        rphi = math.exp(c.cln * hr[i])
        gm = cg * math.exp(c.cln * hs[i])
        # Runge-Kutta step in the volume radius.
        dr = rphi - rphip
        deta1 = _deta(dr, rhop, rhobp, eta2[i - 1], rphip)
        rhoa = 0.5 * (rho + rhop)
        rhoba = 0.5 * (rhobar + rhobp)
        etaa = eta2[i - 1] + 0.5 * deta1
        rphia = rphip + 0.5 * dr
        deta2 = _deta(dr, rhoa, rhoba, etaa, rphia)
        etaa = eta2[i - 1] + 0.5 * deta2
        deta3 = _deta(dr, rhoa, rhoba, etaa, rphia)
        etaa = eta2[i - 1] + deta3
        deta4 = _deta(dr, rho, rhobar, etaa, rphi)
        eta2[i] = eta2[i - 1] + c.cc13 * (0.5 * deta1 + deta2 + deta3
                                          + 0.5 * deta4)
        rphi3 = rphi ** 3
        r03 = rphi3
        r0est = rphi
        acc = nml.acfpft ** c.cc13
        for k in range(1, nml.itfp1 + 1):
            fact = 5.0 * c.cc13 * omega[i] ** 2 / (gm * (2.0 + eta2[i]))
            for j in range(1, nml.itfp2 + 1):
                a = fact * r03
                delr03 = ((rphi3 - r03 * (1.0 + c1 * a ** 2 - c2 * a ** 3))
                          / (1.0 + c3 * a ** 2 - c4 * a ** 3))
                r03 = r03 + delr03
                if abs(delr03 / r03) < nml.acfpft:
                    break
            r0[i] = r03 ** c.cc13  # label 50
            err = r0[i] - r0est
            # Question (see the page): the two tolerance tests use one ACC.
            if abs(err) / r0[i] <= acc:
                break
            # Runge-Kutta step again, in r_0.
            # Question (see the page): RHOA, RHOBA are the first pass's.
            eta2t = eta2[i]
            dr = r0[i] - r0[i - 1]
            deta1 = _deta(dr, rhop, rhobp, eta2[i - 1], r0[i - 1])
            etaa = eta2[i - 1] + 0.5 * deta1
            r0a = r0[i - 1] + 0.5 * dr
            deta2 = _deta(dr, rhoa, rhoba, etaa, r0a)
            etaa = eta2[i - 1] + 0.5 * deta2
            deta3 = _deta(dr, rhoa, rhoba, etaa, r0a)
            etaa = eta2[i - 1] + deta3
            deta4 = _deta(dr, rho, rhobar, etaa, r0[i])
            eta2[i] = eta2[i - 1] + c.cc13 * (0.5 * deta1 + deta2 + deta3
                                              + 0.5 * deta4)
            err = eta2t - eta2[i]
            if abs(err) <= acc:
                break
            r0est = r0[i]
        # label 70
        rhop = rho
        rhobp = rhobar
        rphip = rphi
