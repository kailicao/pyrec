"""Python version of vcirc.f (YREC 09c3b9b): the velocities of the
rotational mixing processes (meridional circulation, the GSF instability,
secular and dynamical shear) at each interface of a radiative region."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import difad2_mod
from pyrec import difad3_mod
from pyrec import dimensions_mod
from pyrec import egridder_mod
from pyrec import intvar_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import oldrot2_mod
from pyrec import quadd_mod
from pyrec import quadru_mod
from pyrec import taukh_mod
from pyrec import temp2_mod
from pyrec import vfact_mod

JSON = dimensions_mod.json

# Local arrays kept between calls (SAVE): the midpoint values (formerly
# COMMON/INTVR2/), VMU of the previous iteration (formerly COMMON/PREVMU/)
# and VMU2.
viscmi = np.zeros(JSON + 1)
thdifmi = np.zeros(JSON + 1)
amumi = np.zeros(JSON + 1)
wm = np.zeros(JSON + 1)
vmup = np.zeros(JSON + 1)
vmu2 = np.zeros(JSON + 1)

# The speed of light. Question (see the page): its DATA statement is
# commented out, so the Fortran's CC is undefined; 0.0 stands for it here
# (TKHS is then zero, and the division by it raises ZeroDivisionError
# where the Fortran would give an infinity).
CC = 0.0


def _interp4(f, faci, i, j0):
    # 4-point Lagrangian interpolation in the logarithm of f, from shells
    # j0..j0+3 with the weights FACI(1..4, I).
    return math.exp(math.log(f[j0]) * faci[1, i]
                    + math.log(f[j0 + 1]) * faci[2, i]
                    + math.log(f[j0 + 2]) * faci[3, i]
                    + math.log(f[j0 + 3]) * faci[4, i])


def vcirc(hr, hru, imin, imax, it, lcz, m, omega, ldo, qwlnr, qwrmax, hv,
          cltot, dt, hp, faci):
    """Estimate the velocities of meridional circulation, the GSF
    instability and secular and dynamical shear at each radiative
    interface imin..imax (between shells i-1 and i), inhibited by mu
    gradients and scaled by the user's factors; their sum hv sets the
    diffusion coefficient for angular momentum and mixing.

    Translates: vcirc.f#VCIRC (09c3b9b)

    Arguments: log r hr and r hru, the range imin..imax, the iteration it,
    the convective flags lcz, the number of shells m, omega, ldo (set),
    the shear qwlnr and its dynamical limit qwrmax, hv (filled in place),
    the total luminosity cltot, the time step dt, log P hp, and the
    interpolation weights faci[1..4, i]. The velocities also go to
    temp2_mod, difad2_mod, difad3_mod and quadd_mod.rat.

    Returns ldo: True if any hv is positive.
    """
    n = nmlpar_mod
    c = constants_mod
    cln = c.cln
    t2 = temp2_mod
    d2 = difad2_mod
    d3 = difad3_mod
    o2 = oldrot2_mod
    tk = taukh_mod
    vf = vfact_mod
    iv = intvar_mod
    rat = quadd_mod.rat
    quad = quadru_mod.quad
    gg = quadru_mod.gg
    amum = mdphy_mod.amum
    viscm = mdphy_mod.viscm
    thdifm = mdphy_mod.thdifm
    ves, vgsf, vss, vmu = t2.ves, t2.vgsf, t2.vss, t2.vmu
    rm = iv.rm

    ldo = False
    # Store the old velocity estimates for later use.
    if it > 1:
        for i in range(imin, imax + 1):  # label 5
            vmup[i] = vmu[i]
            t2.vgsf0[i] = vgsf[i]
            t2.ves0[i] = ves[i]
            t2.vss0[i] = vss[i]
        if n.ldifad:
            for i in range(imin, imax + 1):
                d2.vesa0[i] = d2.vesa[i]
                d2.vesd0[i] = d2.vesd[i]
    for i in range(1, m + 1):  # label 10
        hv[i] = 0.0
        ves[i] = 0.0
        vgsf[i] = 0.0
        vss[i] = 0.0
        vmu[i] = 0.0

    # Viscosity and thermal diffusivity at the interfaces, by 4-point
    # Lagrangian interpolation in their logarithms.
    # First interface.
    if imin < 3:
        viscmi[2] = _interp4(viscm, faci, 2, 1)
        thdifmi[2] = _interp4(thdifm, faci, 2, 1)
        i0 = 3
    else:
        i0 = imin
    # Last interface.
    if imax == m:
        viscmi[m] = _interp4(viscm, faci, m, m - 3)
        thdifmi[m] = _interp4(thdifm, faci, m, m - 3)
        i1 = m - 1
    else:
        i1 = imax
    # General case.
    for i in range(i0, i1 + 1):  # label 20
        viscmi[i] = _interp4(viscm, faci, i, i - 2)
        thdifmi[i] = _interp4(thdifm, faci, i, i - 2)
    # Linear interpolation for omega and mu.
    for i in range(2, m + 1):
        amumi[i] = 0.5 * (amum[i] + amum[i - 1])
        wm[i] = 0.5 * (omega[i] + omega[i - 1])

    # Meridional circulation.
    if n.ies == 2:
        # Kippenhahn & Moellenhoff (1974), restored MHP 4/26.
        for i in range(imin, imax + 1):  # label 30
            if lcz[i] and lcz[i - 1]:
                continue
            dr = hru[i] - hru[i - 1]
            dv = ((omega[i] * hru[i])**2
                  - (omega[i - 1] * hru[i - 1])**2) / dr
            gmid = 0.5 * (gg[i] + gg[i - 1])
            phistd = 2.0 * c.cc23 * wm[i]**2 * rm[i] / gmid
            phis2 = c.cc13 * wm[i]**2 * (hru[i]**2 / gg[i]
                                         - hru[i - 1]**2 / gg[i - 1]) / dr
            phi2 = (quad[i - 1] / gg[i - 1] - quad[i] / gg[i]) / dr
            rat[i] = (phis2 + phi2) / phistd
            # Question (see the page): overwritten at once, dropping the
            # quadrupole term; VESA, VESD and THN are not set here.
            rat[i] = phis2 / phistd
            ves[i] = abs(vf.fact1[i] * (vf.fact2[i] * rat[i] * wm[i]**2
                                        + tk.fact6[i] * dv))
    else:
        # Zahn (1991), with the quadrupole term of Zahn (1992).
        for i in range(imin, imax + 1):  # label 31
            if lcz[i] and lcz[i - 1]:
                continue
            ves[i] = wm[i]**2 * (tk.fes1[i] + wm[i]**2 * tk.fes2[i])
            dr = hru[i] - hru[i - 1]
            gmid = 0.5 * (gg[i] + gg[i - 1])
            phistd = 2.0 * c.cc23 * wm[i]**2 * rm[i] / gmid
            phis2 = c.cc13 * wm[i]**2 * (hru[i]**2 / gg[i]
                                         - hru[i - 1]**2 / gg[i - 1]) / dr
            phi2 = (quad[i - 1] / gg[i - 1] - quad[i] / gg[i]) / dr
            qw = wm[i] * (omega[i] - omega[i - 1]) / dr
            rat[i] = (phis2 + phi2) / phistd
            if rat[i] < 0.0:
                print(f"{i:5d} RAT {rat[i]:10.3E} POT SPH,CYL,QUA"
                      f"{phistd:10.3E}{phis2:10.3E}{phi2:10.3E}\n"
                      f" QUA,G,R {quad[i]:12.3E}{quad[i - 1]:12.3E}"
                      f"{gg[i]:12.3E}{gg[i - 1]:12.3E}{hru[i]:12.3E}"
                      f"{hru[i - 1]:12.3E}")  # FORMAT 303
            if not n.ldifad:
                ves[i] = abs(rat[i] * ves[i] - tk.fes3[i] * qw)
            else:
                # MHP 06/02: the d theta/dt term; MHP 05/02: factor 1/5.
                o2.thn[i] = o2.festime[i] * (o2.thm[i] * qw - o2.thp[i]) / dt
                d2.vesa[i] = 0.2 * (rat[i] * ves[i] + o2.thn[i])
                q1 = vf.fact7[i] * wm[i]**2
                q2 = vf.fact8[i] * qw
                q0 = tk.fes3[i] * wm[i]**2
                d2.vesd[i] = 0.2 * (q0 + q1 + q2)
                # Second- and third-order terms.
                d3.vesd2[i] = 0.2 * d3.facd2[i] * wm[i]
                d3.vesd3[i] = 0.2 * d3.facd3[i] * wm[i]
                ves[i] = (rat[i] * ves[i] + o2.thn[i]
                          + (tk.fes3[i] + vf.fact7[i]) * qw
                          + vf.fact8[i] * ((omega[i] - omega[i - 1]) / dr)**2)
                egridder_mod.eqqcod[i] = egridder_mod.qqcod[i] * wm[i]

    # Inhibition by mu gradients.
    if n.imu == 3:
        # Kippenhahn & Moellenhoff (1974): fictitious mu currents that
        # oppose the circulation, restored MHP 4/26.
        for i in range(imin, imax + 1):
            if lcz[i] and lcz[i - 1]:
                continue  # GOTO 32
            # Square root of tKH* (Kippenhahn, IAU 66, p. 23); the factor
            # of g is in FACT4.
            tkhs = math.sqrt(viscmi[i] * CC / thdifmi[i])
            qmu = (amum[i] - amum[i - 1]) / (hru[i] - hru[i - 1])
            vmu[i] = abs(n.fmu * vf.fact1[i] * vf.fact4[i] * amumi[i] * qmu
                         / tkhs)
            vmu2[i] = n.fmu * iv.hgm[i] * abs(qmu) / amumi[i] / wm[i]**2
    else:
        # Zahn & Maeder (1998): (del - del_ad) replaced by
        # (del + del_mu - del_ad).
        # Question (see the page): convective interfaces are not skipped,
        # and FESTIME and FV0 are rescaled at every call.
        for i in range(imin, imax + 1):
            qmu = math.log(amum[i]) - math.log(amum[i - 1])
            qp = cln * (hp[i] - hp[i - 1])
            ddel = o2.ddelm[i] + qmu / qp
            ddtest = max(o2.ddelm[i], 1.0e-3)
            ddtest2 = max(ddel, 1.0e-3)
            qqq = ddtest / ddtest2
            ves[i] = ves[i] * qqq
            if n.ldifad:
                d2.vesa[i] = d2.vesa[i] * qqq
                d2.vesd[i] = d2.vesd[i] * qqq
                d3.vesd2[i] = d3.vesd2[i] * qqq
                d3.vesd3[i] = d3.vesd3[i] * qqq
            q1 = vf.fact7[i] * wm[i]**2 * qqq
            q0 = tk.fes3[i] * wm[i]**2 * qqq
            dr = hru[i] - hru[i - 1]
            v2 = 0.2 * (q0 + q1) * (omega[i] - omega[i - 1]) / dr / wm[i]  # unused
            o2.festime[i] = o2.festime[i] * qqq
            d3.fv0[i] = d3.fv0[i] * qqq

    # GSF instability: the Kippenhahn (1980) estimate. IGSF = 0: inhibited
    # by viscosity and mu gradients; 1: not inhibited; 2: by viscosity only.
    for i in range(imin, imax + 1):  # label 40
        if lcz[i] and lcz[i - 1]:
            continue
        if n.igsf == 2 or n.igsf == 0:
            qwrmx = 2.0 * math.sqrt(viscmi[i] / thdifmi[i]) * qwrmax[i]
            if abs(qwlnr[i]) < qwrmx:
                vgsf[i] = 0.0
                continue
            fxx = math.sqrt((abs(qwlnr[i]) - qwrmx) / qwrmx)
        else:
            fxx = 1.0
        dr = hru[i] - hru[i - 1]
        if n.igsf == 0:
            qwrmx = 2.0 * math.sqrt(iv.hgm[i] * abs(amum[i] - amum[i - 1])
                                    / dr / amumi[i])
            if abs(qwlnr[i]) < qwrmx:
                vgsf[i] = 0.0
                continue
            # MHP 05/02: tests to avoid dividing by zero.
            qq = 1.0e-2 * qwlnr[i]
            if abs(qwrmx) < abs(qq):
                fxx2 = 1.0
            elif abs(qwrmx) < 1.0e-20:
                fxx2 = 1.0
            else:
                fxx2 = math.sqrt((abs(qwlnr[i]) - qwrmx) / qwrmx)
            fxx = min(fxx, fxx2)
        rmid = 0.5 * (hru[i] + hru[i - 1])
        dlnwdr = abs(math.log(omega[i]) - math.log(omega[i - 1])) / dr
        # If d ln(j)/dr < 0, use the larger of the two criteria; otherwise
        # only d omega/dz nonzero. The factor 1/r: all the velocities are
        # later multiplied by r.
        if omega[i] * hru[i]**2 < omega[i - 1] * hru[i - 1]**2:
            dlnjmdr = abs(2.0 / rmid + (math.log(omega[i])
                                        - math.log(omega[i - 1])) / dr)
            fx = max(2.0 * dlnjmdr, 0.25 * dlnwdr)
            vgsf[i] = tk.fgsfk[i] * fx * thdifmi[i] * wm[i]**2 / rmid
        else:
            vgsf[i] = 0.25 * tk.fgsfk[i] * thdifmi[i] * dlnwdr * wm[i]**2 \
                / rmid
        vgsf[i] = abs(fxx * vgsf[i])

    # Secular and dynamical shear (Endal & Sofia, paper II).
    for i in range(imin, imax + 1):  # label 60
        if abs(qwlnr[i]) > qwrmax[i]:
            # Dynamical shear: a (large) velocity, and the GSF velocity
            # rescaled to the gradient at the limit.
            qwr = abs(qwlnr[i])
            vss[i] = 8.0 / 4.5e1 * thdifmi[i] * (qwr / qwrmax[i])**2 / rm[i]
            dr = hru[i] - hru[i - 1]
            rmid = 0.5 * (hru[i] + hru[i - 1])
            if omega[i] < omega[i - 1]:
                wmin = omega[i] + qwrmax[i] * cln * (hr[i] - hr[i - 1])
                wmid = omega[i] + qwrmax[i] * cln * (hr[i] - math.log10(rmid))  # unused
            else:
                wmin = omega[i] - qwrmax[i] * cln * (hr[i] - hr[i - 1])
                wmid = omega[i] - qwrmax[i] * cln * (hr[i] - math.log10(rmid))  # unused
            # Question (see the page): DLNWDR0 is zero if the two omegas
            # are equal (ZeroDivisionError here, an infinity or NaN in the
            # Fortran).
            dlnwdr0 = abs(math.log(omega[i]) - math.log(omega[i - 1])) / dr
            dlnjmdr0 = abs(2.0 / rmid + (math.log(omega[i])
                                         - math.log(omega[i - 1])) / dr)
            dlnjmdr = abs(2.0 / rmid + (math.log(omega[i])
                                        - math.log(wmin)) / dr)
            dlnwdr = abs(math.log(omega[i]) - math.log(wmin)) / dr
            if omega[i] * hru[i]**2 < wmin * hru[i - 1]**2:
                fx = max(2.0 * dlnjmdr, 0.25 * dlnwdr)
                fx2 = max(2.0 * dlnjmdr0, 0.25 * dlnwdr0)
                vgsf[i] = vgsf[i] * fx / fx2
            else:
                vgsf[i] = vgsf[i] * dlnwdr / dlnwdr0
            # Question (see the page): to unit 6 at every unstable interface.
            print(f" DYNAMICAL SHEAR-SHELL{i:5d} WTOP{omega[i]:11.3E}"
                  f" WBOT{omega[i - 1]:11.3E} LIMIT{wmin:11.3E}")  # 9911
            continue
        # Secular shear: the maximum gradient it allows, with the square
        # root of the Prandtl number (MHP 3/92); velocity from Zahn 1991.
        qwrmx = math.sqrt(viscmi[i] / thdifmi[i] * 1.25e-1 * n.rcrit) \
            * qwrmax[i]
        if abs(qwlnr[i]) > qwrmx:
            if abs((amum[i] - amum[i - 1]) / amumi[i]) < 1.0e-10:
                qwrmx2 = 0.0
                qwr = abs(qwlnr[i]) - qwrmx
                vss[i] = 8.0 / 4.5e1 * thdifmi[i] * (qwr / qwrmax[i])**2 \
                    / rm[i]
            else:
                # The limit set by the mu gradient; the larger one is used.
                qwrmx2 = 2.0 * math.sqrt(max(1.0e-20, vf.fact5[i] * abs(
                    (amum[i] - amum[i - 1]) / amumi[i])))
                if abs(qwlnr[i]) > qwrmx2:
                    qwrmx = max(qwrmx2, qwrmx)
                    qwr = abs(qwlnr[i]) - qwrmx
                    vss[i] = 8.0 / 4.5e1 * thdifmi[i] \
                        * (qwr / qwrmax[i])**2 / rm[i]

    # The user's scale factors (and, with IMU = 3, the mu currents).
    if n.imu == 3:
        if it > 1:
            for i in range(1, m + 1):
                vmu[i] = 0.5 * (vmu[i] + vmup[i])
        for i in range(imin, imax + 1):
            # Question (see the page): VMU2 is stale at convective
            # interfaces.
            fm = 1.0 + vmu2[i]
            rmid = 0.5 * (hru[i] + hru[i - 1])
            fcc = math.sqrt(n.fc * n.fesc)
            vest = max(0.0, n.fes * ves[i]
                       - math.sqrt(n.fes * ves[i] * rmid) * fcc * vmu[i])
            if vest <= 0.0:
                ves[i] = abs(n.fes * ves[i] / fm)
            else:
                ves[i] = vest
            fcc = math.sqrt(n.fc * n.fgsfc)
            if vgsf[i] > 0.0:
                vgsft = max(0.0, n.fgsf * vgsf[i]
                            - math.sqrt(n.fgsf * vgsf[i] * rmid) * fcc
                            * vmu[i])
            else:
                vgsft = 0.0
            if vgsft <= 0.0:
                # Question (see the page): FES, not FGSF.
                vgsf[i] = abs(n.fes * vgsf[i] / fm)
            else:
                vgsf[i] = vgsft
            vss[i] = max(0.0, n.fss * vss[i])
    else:
        for i in range(imin, imax + 1):
            ves[i] = abs(n.fes * ves[i])
            vgsf[i] = n.fgsf * vgsf[i]
            vss[i] = n.fss * vss[i]
    if n.ldifad:
        for i in range(imin, imax + 1):
            d2.vesd[i] = abs(n.fes * d2.vesd[i])
            d2.vesa[i] = n.fes * d2.vesa[i]
    if n.ldifad:
        # MHP 05/02: the GSF and shear velocities, times r, join the
        # diffusive part.
        for i in range(imin, imax + 1):
            d2.vesd[i] = d2.vesd[i] + rm[i] * (abs(vgsf[i]) + abs(vss[i]))
            d3.dss[i] = rm[i] * abs(vss[i])
            d3.dgsf[i] = rm[i] * abs(vgsf[i])

    # Average the previous and new estimates after the first iteration.
    if it > 1:
        for i in range(imin, imax + 1):  # label 70
            vgsf[i] = 0.5 * (vgsf[i] + t2.vgsf0[i])
            ves[i] = 0.5 * (ves[i] + t2.ves0[i])
            vss[i] = 0.5 * (vss[i] + t2.vss0[i])
        if n.ldifad:
            for i in range(imin, imax + 1):
                d2.vesa[i] = 0.5 * (d2.vesa[i] + d2.vesa0[i])
                d2.vesd[i] = 0.5 * (d2.vesd[i] + d2.vesd0[i])
    for i in range(imin, imax + 1):  # label 80
        hv[i] = vgsf[i] + ves[i] + vss[i]
        if hv[i] < 1.0e-20:
            hv[i] = 0.0
        if hv[i] > 0.0:
            ldo = True

    # 9/93: no transport in the core (L(r)/L < 0.99), to avoid numerical
    # problems in the H-burning shell.
    # Question (see the page): the threshold is hard-coded.
    if n.lnoj:
        for i in range(imin, imax + 1):
            dell = iv.alm[i] / cltot
            if dell < 9.9e-1:
                hv[i] = 0.0
                vgsf[i] = 0.0
                ves[i] = 0.0
                vss[i] = 0.0
                d2.vesa[i] = 0.0
                d2.vesd[i] = 0.0
            else:
                break  # GOTO 81

    return ldo
