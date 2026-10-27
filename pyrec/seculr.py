"""Python version of seculr.f (YREC development, 09c3b9b)."""

import math
from types import SimpleNamespace

import numpy as np

from pyrec import constants_mod
from pyrec import dimensions_mod
from pyrec import disk_mod
from pyrec import egrid_mod
from pyrec import intvar_mod
from pyrec import nmlpar_mod
from pyrec import oldrot2_mod
from pyrec.checkc import checkc
from pyrec.checkj import checkj
from pyrec.codiff import codiff
from pyrec.dadcoeft import dadcoeft
from pyrec.dcoeft import dcoeft
from pyrec.eq2mod import eq2mod
from pyrec.getfc import getfc
from pyrec.getqua import getqua
from pyrec.mcowind import mcowind
from pyrec.mixcz import mixcz
from pyrec.mwind import mwind
from pyrec.ndifcom import ndifcom
from pyrec.rotgrid import rotgrid
from pyrec.setupv import setupv
from pyrec.solid import solid
from pyrec.tridia import tridia
from pyrec.vcirc import vcirc

# Fortran SAVEs every local of SECULR: the work arrays below, and the
# scalars IBEG, IEND, DR, SUMDJ, IDM, ... keep their values between
# calls. Some of them are read before they are set in a call (question 1
# on the page), so they live here rather than inside the function.
_JSON = dimensions_mod.json
_s = SimpleNamespace(
    qwrmax=np.zeros(_JSON + 1), em=np.zeros(_JSON + 1),
    cod=np.zeros(_JSON + 1), ecod2=np.zeros(_JSON + 1),
    ei=np.zeros(_JSON + 1), ej=np.zeros(_JSON + 1),
    ecod=np.zeros(_JSON + 1), hru=np.zeros(_JSON + 1),
    rmid=np.zeros(_JSON + 1), qwlnr=np.zeros(_JSON + 1),
    dj=np.zeros(_JSON + 1), hjmp=np.zeros(_JSON + 1),
    wsav=np.zeros(_JSON + 1), ew=np.zeros(_JSON + 1),
    id=np.zeros(_JSON + 1, dtype=int),
    faci=np.zeros((4 + 1, _JSON + 1)),   # formerly COMMON/INTFAC/
    vfc=np.zeros(_JSON + 1),             # formerly COMMON/VARFC/
    ibeg=0, iend=0, istart=0, dr=0.0, sumdj=0.0, idm=0,
    ldum=False, ldum2=False, lokad=False,
    wind1=0.0, wind2=0.0, lfix=False,
)


def seculr(dt, hd, hg, hi, hl, hp, hr, hs, hs1, hs2, ht,
           hjmsav, lcz, m, omega,
           hjm, eta2, qiw, r0, hcomp, imin, imax,
           bl, hstot, smass, teffl, lredo, iredo,
           hicz, hsbot, hstop, omegas, lczsur,
           cod2, hv, lok, hcomps):
    """Secular transport of angular momentum and composition, one step.

    Transfers angular momentum (and the associated transport of
    composition) by secular rotational instabilities: secular and
    dynamical shear, Eddington circulation and the GSF instability. Sets
    up, iterates and checks the diffusion of Omega and of the abundances,
    with the wind applied to the surface convection zone.

    Translates: seculr.f#SECULR (09c3b9b)

    Parameters
    ----------
    dt : float
        Time step.
    hd, hg, hi, hl, hp, hr, hs, hs1, hs2, ht : ndarray, shape (JSON + 1,)
        The model, by shell.
    hjmsav : ndarray, shape (JSON + 1,)
        Specific angular momentum at the start of the step.
    lcz : ndarray of bool, shape (JSON + 1,)
        Convection flags; ``lcz[1]`` and ``lcz[m]`` are changed during the
        call and restored.
    m : int
        Number of shells.
    omega, hjm : ndarray, shape (JSON + 1,)
        Angular velocity and specific angular momentum, updated in place.
    eta2, qiw, r0 : ndarray, shape (JSON + 1,)
        Rotation quantities passed to SOLID and CHECKJ.
    hcomp : ndarray, shape (16, JSON + 1)
        Abundances ``hcomp[k, j]``, updated in place.
    imin, imax : int
        Range of radiative interfaces.
    bl, hstot, smass, teffl : float
        log L/Lsun, log total mass, mass in solar units, log Teff (wind).
    lredo, iredo, hicz, omegas, lok : bool, int, float, float, bool
        Input values are kept where the routine does not set them.
    hsbot, hstop : float
        Masses at the bottom and top of the surface convection zone.
    lczsur : bool
        The star has a surface convection zone.
    cod2, hv : ndarray, shape (JSON + 1,)
        Composition diffusion coefficients and circulation velocities,
        set in place.
    hcomps : ndarray, shape (16, JSON + 1)
        Abundances at the start of the step (formerly COMMON/OLDAB/).

    Returns
    -------
    dt : float
        The step, halved by CHECKC (or CHECKJ) when it asks for a cut.
    lredo : bool
        The step must be cut.
    iredo : int
        Why, as set by CHECKJ or CHECKC.
    hicz : float
        Moment of inertia of the surface convection zone.
    omegas : float
        Surface angular velocity.
    lok : bool
        Converged.
    """
    c = constants_mod
    nml = nmlpar_mod
    s = _s
    qwlnr, hru, faci, vfc = s.qwlnr, s.hru, s.faci, s.vfc
    wsav, hjmp = s.wsav, s.hjmp

    # Printing of the composition changes is switched off (G Somers 11/14).
    lprt = False
    # Disk locking (MHP 9/94).
    lock = False
    if nml.ldisk and disk_mod.sage <= nml.tdisk:
        lock = True
    # The first and last shells are treated as convective for diffusion.
    lcz1 = lcz[1]
    lczm = lcz[m]
    # Set up factors used to compute diffusion velocities.
    setupv(hd, hg, hl, hp, hr, hs1, ht, imin, imax, m, hru, s.qwrmax, faci)
    cltot = c.clsun * hl[m]
    # Angular velocity gradients, centered logarithmic derivative (label 10).
    for i in range(imin, imax + 1):
        rlnc = math.log10(intvar_mod.rm[i])
        dr1 = 1.0 / c.cln / (hr[i] - rlnc) + 1.0 / c.cln / (rlnc - hr[i - 1])
        qwlnr[i] = 0.25 * (omega[i] - omega[i - 1]) * dr1
    getqua(hd, hg, hru, omega, m)
    vfc[1:m + 1] = 0.0
    # Stability of the model and circulation velocities.
    it = 1
    ldo = vcirc(hr, hru, imin, imax, it, lcz, m, omega, False, qwlnr,
                s.qwrmax, hv, cltot, dt, hp, faci)
    if nml.lvfc:
        getfc(hd, hru, hv, imin, imax, omega, vfc)
    # Store the initial angular momentum distribution (label 20).
    hjmp[1:m + 1] = hjmsav[1:m + 1]
    wsav[1:m + 1] = omega[1:m + 1]
    # Initial surface angular velocity, for the angular momentum loss.
    wbeg = oldrot2_mod.wmst[m]
    lok = False

    converged = False
    for it in range(1, nml.itdif2 + 1):  # label 100
        omegas = omega[m]
        if it > 1:
            # New Omega: average of the initial and the latest values.
            for i in range(1, m + 1):
                omega[i] = 0.5 * (omega[i] + wsav[i])
            # New angular velocity gradients (label 30).
            for i in range(imin, imax + 1):
                rlnc = math.log10(intvar_mod.rm[i])
                dr1 = (1.0 / c.cln / (hr[i] - rlnc)
                       + 1.0 / c.cln / (rlnc - hr[i - 1]))
                qwlnr[i] = 0.25 * (omega[i] - omega[i - 1]) * dr1
            ldo = vcirc(hr, hru, imin, imax, it, lcz, m, omega, ldo, qwlnr,
                        s.qwrmax, hv, cltot, dt, hp, faci)
            if nml.lvfc:
                getfc(hd, hru, hv, imin, imax, omega, vfc)
            # Reset angular momentum and species 1-4 to their start values
            # (labels 40, 50).
            for i in range(1, m + 1):
                hjm[i] = hjmsav[i]
                for j in range(1, 5):
                    hcomp[j, i] = hcomps[j, i]
            # Remix convection zones to their proper depth (MHP 10/91).
            # Question (see the page): LIMIX is set but no longer used.
            limix = False  # noqa: F841
            mixcz(hcomp, hs2, lcz, m)

        # Wind loss for an isolated surface convection zone.
        if lczsur and hv[imax] == 0.0 and not lock:
            hicz = 0.0
            for i in range(imax, m + 1):  # label 55
                hicz = hicz + hi[i]
            ljdot = nml.ljdot0
            mwind(bl, dt, hsbot, hstop, imax, m, ljdot, omegas,
                  smass, teffl, hicz, hjm)
        # Remove the torque from the entire star (IMPJMOD = 1 as LSOLID).
        elif lczsur and (nml.lsolid or nml.impjmod == 1):
            ljdot = nml.ljdot0
            hsbotj = 0.0
            hstopj = math.exp(c.cln * hstot)
            imaxj = 1
            mwind(bl, dt, hsbotj, hstopj, imaxj, m, ljdot, omegas,
                  smass, teffl, hicz, hjm)
            jstart = 1
            jend = m
            solid(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega,
                  qiw, r0, m)

        # No instability: done (label 9999).
        if not ldo:
            return dt, lredo, iredo, hicz, omegas, lok

        # Treat the central and surface zones as convective.
        lcz[1] = True
        lcz[m] = True
        # Diffusion coefficients for all unstable regions.
        codiff(hru, m, s.rmid, s.cod, cod2, vfc)

        # Each unstable region is solved separately. IBEG is the zone
        # below the first non-zero v, IEND the zone above the last.
        s.ldum = False
        s.istart = imin
        while True:  # label 60
            ltest = False
            if nml.lcodm and nml.codm > 0.0:
                # Treat the entire domain as unstable (MHP 8/13).
                # Question (see the page): not gated on solid-body models.
                ltest = True
                s.ibeg = imin - 1
                s.iend = imax
                s.istart = imax + 1
            else:
                found_end = False
                for j in range(s.istart, imax + 1):  # label 70
                    if hv[j] > 0.0:
                        ltest = True
                        if not s.ldum:
                            # Start of an unstable region.
                            s.ldum = True
                            s.ibeg = j - 1
                    elif s.ldum:
                        # End of an unstable region.
                        s.iend = j - 1
                        s.ldum = False
                        s.istart = j + 1
                        found_end = True
                        break
                if not found_end:
                    # The last interface is unstable: set IEND properly.
                    if s.ldum:
                        s.iend = imax
                    s.istart = imax + 1
            # label 80: no non-zero v encountered, exit.
            if not ltest:
                break  # to label 90

            # Equal spacing in CHI for the region.
            s.dr, s.ldum2 = rotgrid(
                s.cod, cod2, hd, hi, hjmsav, hl, hp, hr, hru, hs, hs1, hs2,
                hstot, s.ibeg, s.iend, lcz, m, wsav, s.dr, s.ecod,
                s.ecod2, s.ei, s.ej, s.em, s.ew, s.ldum2)
            # Two convection zones separated by one radiative zone: skip.
            if s.ldum2:
                if s.istart <= imax:
                    continue  # to label 60
                break  # to label 90

            # Skip the angular momentum evolution for solid-body models
            # (MHP 3/09; IMPJMOD = 1 as LSOLID).
            if not nml.lsolid and nml.impjmod != 1:
                # Surface convection zone in the region: wind terms.
                if lczsur and s.iend == imax:
                    if lock:
                        s.wind1 = 0.0
                        s.wind2 = 0.0
                        s.lfix = True
                    else:
                        s.lfix = False
                        if nml.ljdot0:
                            hicz = s.ei[egrid_mod.ntot]
                            omegas = omega[m]
                            s.wind1, s.wind2 = mcowind(
                                bl, dt, hicz, it, omegas, smass,
                                teffl, wbeg, s.wind1, s.wind2)
                        else:
                            s.wind1 = 0.0
                            s.wind2 = 0.0
                else:
                    s.lfix = False
                    s.wind1 = 0.0
                    s.wind2 = 0.0
                s.lokad = False
                if not nml.ldifad:
                    # Diffusion equation for Omega at the end of the step,
                    # solved on the equally spaced grid.
                    dcoeft(s.ecod, s.dr, dt, s.ei, s.ej, s.ew,
                           egrid_mod.ntot, s.wind1, s.wind2, s.lfix)
                    s.sumdj = tridia(egrid_mod.ntot, s.ei, s.dj, s.sumdj)
                else:
                    # Omega and its derivatives in a band matrix.
                    dt, s.wind2, s.sumdj, s.lokad = dadcoeft(
                        s.dr, dt, s.ei, s.ew, egrid_mod.ntot, s.wind1,
                        s.wind2, s.dj, s.ecod2, s.sumdj, s.lfix, s.lokad)
                # Back to the original grid.
                eq2mod(s.dj, s.ej, hs2, s.ibeg, s.iend, lcz, m, s.sumdj, hjm)

            # Composition diffusion of species 1-4 until the final
            # iteration. Question (see the page): the Fortran comment
            # names only H, He3 and He4.
            jbeg = 1
            jend = 4
            s.ibeg, s.iend = ndifcom(
                dt, cod2, s.ecod2, s.em, hd, hl, hp, hr, hs, hs1, hs2,
                hstot, hv, s.ibeg, s.iend, imax, imin, lcz, lok, m,
                hcomp, jbeg, jend)
            # Next region, if any.
            if s.istart <= imax:
                continue  # to label 60
            break

        # label 90: check the solution, update Omega, decide whether
        # another iteration is needed.
        lcz[1] = lcz1
        lcz[m] = lczm
        dt, iredo, lok, lredo, s.idm = checkj(
            hd, hjmp, hjmsav, hr, hs, hs2, hv, imin, imax, it, lcz,
            m, dt, eta2, hi, hjm, iredo, lok, lredo, omega, qiw, r0,
            wsav, s.id, s.idm, s.lokad)
        # Check composition diffusion; recompute mean molecular weight.
        if not lredo:
            dt, iredo, lok, lredo = checkc(hcomp, it, lprt, m, dt,
                                           iredo, lok, lredo)
        if nml.lnoj:
            lok = True
        if lok:
            converged = True
            break  # to label 200
        if lredo:
            return dt, lredo, iredo, hicz, omegas, lok  # label 9999
    if not converged:
        # Question (see the page): non-convergence counts as success.
        lok = True

    # label 200: composition diffusion of the remaining species.
    lcz[1] = True
    lcz[m] = True
    jbeg = 5
    if nml.lexcom:
        jend = 15
    else:
        jend = 11
    s.ibeg, s.iend = ndifcom(
        dt, cod2, s.ecod2, s.em, hd, hl, hp, hr, hs, hs1, hs2,
        hstot, hv, s.ibeg, s.iend, imax, imin, lcz, lok, m,
        hcomp, jbeg, jend)
    lcz[1] = lcz1
    lcz[m] = lczm

    # Output of the composition changes to unit IMODPT (iounits_mod.imodpt;
    # standard output stands in for it). Never runs: lprt is False.
    if lprt:
        print('COMPOSITION CHANGE FROM ROTATIONAL MIXING')
        print('    H        HE3       C12       C13       N14       LI7'
              '       BE9')
        for i in range(1, s.idm + 1):
            k = s.id[i]
            print(f'{k:5d}' + ''.join(
                f'{hcomp[j, k] - hcomps[j, k]:10.3E}'
                for j in (1, 4, 5, 6, 7, 14, 15)))

    # label 9999
    return dt, lredo, iredo, hicz, omegas, lok
