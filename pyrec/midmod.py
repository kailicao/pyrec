"""Python version of midmod.f (YREC 09c3b9b): the structure at an
intermediate time of a model step, for the angular-momentum transport
loop."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import deuter_mod
from pyrec import dimensions_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import oldmod_mod
from pyrec import oldphy_mod
from pyrec import oldrot_mod
from pyrec import scrtch_mod
from pyrec import temp_mod
from pyrec.dburnm import dburnm
from pyrec.getrot import getrot
from pyrec.ovrot import ovrot
from pyrec.rotmix import rotmix
from pyrec.solid import solid

JSON = dimensions_mod.json

# The Fortran's local arrays (SAVE): the previous intermediate step's values
# (DDELP, DDELN, HRP, LCP, formerly COMMON/PREVMID/), the changed-shell flags
# and the deuterium rates.
ddelp = np.zeros(JSON + 1)
ddeln = np.zeros(JSON + 1)
hrp = np.zeros(JSON + 1)
lcp = np.zeros(JSON + 1, dtype=bool)
lchg = np.zeros(JSON + 1, dtype=bool)
dratem = np.zeros(JSON + 1)
dratem0 = np.zeros(JSON + 1)


def midmod(delts, dt, fx, hcomp, hd, hg, hjm, hl, hp, hr, hs, hs1,
           hs2, hstot, ht, lc, lfirstit, m, eta2m, hdm, hgm, hicz, him,
           hlm, hpm, hrm, hsbot, hstop, htm, imin, imax, lallcz, lcm,
           lczm, lczsur, omegam, r0m, qiwm, mrzone, mxzone, nrzone,
           nzone, hcomps):
    """Build the structure at the fraction fx of the model step delts by
    linear interpolation in time between the old model (the *_mod arrays)
    and the new one (the arguments); burn and mix the composition over the
    diffusion step dt; find the convection-zone edges; and solve for the
    rotation consistent with the specific angular momentum hjm.

    Translates: midmod.f#MIDMOD (09c3b9b)

    Arguments: the new model (hcomp[k, j], hd, hg, hjm, hl, hp, hr, hs, hs1,
    hs2, ht, lc; log total mass hstot; m shells); lfirstit, True on the
    first diffusion step of the model step. Arrays (index 0 unused) set in
    place: hcomp (reset, burnt and mixed), hjm (only in the dead block),
    eta2m, hdm, hgm, him, hlm, hpm, hrm, htm, lcm, lczm, omegam, r0m, qiwm,
    mrzone, mxzone, hcomps (the composition after mixing) and the
    interpolated physics in mdphy_mod. Scalars set are passed in with any
    value.

    Returns (fx, hicz, hsbot, hstop, imin, imax, lallcz, lczsur, nrzone,
    nzone): fx (changed only in the dead block), the moment of inertia of
    the surface convection zone (unchanged without one), the masses at its
    bottom and top, the edges imin and imax, whether the star is all
    convective, whether a surface zone exists, and the numbers of radiative
    and mixed zones.
    """
    c = constants_mod
    nml = nmlpar_mod
    om_ = oldmod_mod
    op = oldphy_mod
    md = mdphy_mod
    sc = scrtch_mod
    tp = temp_mod
    iend = 11
    if nml.lexcom:
        iend = 15
    # On the first diffusion step, start from the old composition; save the
    # radius, convective state and DELA - DELR of the prior step (MHP 06/02).
    if lfirstit:
        for j in range(1, m + 1):  # DO 20
            lcp[j] = om_.lco[j]
            hrp[j] = om_.hro[j]
            ddelp[j] = op.delao[j] - op.delro[j]
            md.amum[j] = op.amuo[j]
            for i in range(1, iend + 1):
                hcomp[i, j] = om_.hcompp[i, j]
    else:
        for j in range(1, m + 1):
            hrp[j] = hrm[j]
            ddelp[j] = ddeln[j]
            lcp[j] = lcm[j]
    lnewcz = False
    # Interpolate in the model variables and the auxiliary physics.
    # Question (see the page): AMUM is incremented, not interpolated.
    fac2 = dt / delts
    for j in range(1, m + 1):  # DO 40
        hdm[j] = om_.hdo[j] + fx * (hd[j] - om_.hdo[j])
        hlm[j] = om_.hlo[j] + fx * (hl[j] - om_.hlo[j])
        hpm[j] = om_.hpo[j] + fx * (hp[j] - om_.hpo[j])
        hrm[j] = om_.hro[j] + fx * (hr[j] - om_.hro[j])
        htm[j] = om_.hto[j] + fx * (ht[j] - om_.hto[j])
        hgm[j] = oldrot_mod.hgo[j] + fx * (hg[j] - oldrot_mod.hgo[j])
        md.delam[j] = op.delao[j] + fx * (sc.sdel[3, j] - op.delao[j])
        md.delm[j] = op.delo[j] + fx * (sc.sdel[2, j] - op.delo[j])
        md.delrm[j] = op.delro[j] + fx * (sc.sdel[1, j] - op.delro[j])
        md.esumm[j] = op.esumo[j] + fx * (sc.sesum[j] - op.esumo[j])
        md.viscm[j] = op.visco[j] + fx * (tp.visc[j] - op.visco[j])
        md.thdifm[j] = op.thdifo[j] + fx * (tp.thdif[j] - op.thdifo[j])
        md.cpm[j] = op.cpo[j] + fx * (tp.cp[j] - op.cpo[j])
        md.qdtm[j] = op.qdto[j] + fx * (tp.sqdt[j] - op.qdto[j])
        md.om[j] = op.oo[j] + fx * (sc.so[j] - op.oo[j])
        md.amum[j] = md.amum[j] + fac2 * (tp.hamu[j] - op.amuo[j])
        # The total energy generation: nuclear plus the gravitational terms (MHP 6/00).
        temp = sc.sesum[j] + sc.seg[6, j] + sc.seg[7, j]
        md.epsm[j] = op.epso[j] + fx * (temp - op.epso[j])
    # Shells that change state over the step are convective if
    # DELA < DELR at the intermediate time.
    for i in range(1, m + 1):  # DO 50
        ddeln[i] = md.delam[i] - md.delrm[i]
        if lc[i] == lcp[i]:
            lcm[i] = lc[i]
            md.velm[i] = op.velo[i] + fx * (sc.svel[i] - op.velo[i])
            lchg[i] = False
        else:
            lchg[i] = True
            lnewcz = True
            if md.delam[i] < md.delrm[i]:
                lcm[i] = True
                md.velm[i] = max(op.velo[i], sc.svel[i])
            else:
                lcm[i] = False
                md.velm[i] = 0.0
    # Question (see the page): this reset makes the block below dead code.
    lnewcz = False
    if lnewcz:
        fx = _receding_cz(m, fx, hjm, hr, hs2, lcm)
    # Convective overshoot applied to normal convection zones.
    nrzone, nzone = ovrot(hcomp, hdm, hpm, hrm, hs, htm, lcm, m, lczm,
                          mrzone, mxzone, nrzone, nzone)
    # Reburn the original mixture for the short step, at the current T (MHP 6/00).
    rotmix(dt, hcomp, hs2, htm, m, mrzone, mxzone, nrzone, nzone,
           hstot, hdm, hs, hrm, hpm, lczm, hs1)
    for i in range(1, m + 1):
        for j in range(1, iend + 1):
            hcomps[j, i] = hcomp[j, i]
    # Deuterium burning (MHP 05/02).
    if om_.hcompp[12, m] > 1.0e-14:
        dr = deuter_mod
        if lfirstit:
            for i in range(1, m + 1):
                dratem[i] = dr.drate0[i] + fac2 * (dr.drate[i] - dr.drate0[i])
                dratem0[i] = dr.drate0[i]
        else:
            for i in range(1, m + 1):
                dratem0[i] = dratem[i]
                dratem[i] = dratem[i] + fac2 * (dr.drate[i] - dr.drate0[i])
        # Radiative zones, shell by shell, until T drops below the cutoff.
        # Question (see the page): the first cool shell ends both loops.
        done = False
        for k in range(1, nrzone + 1):
            for j in range(mrzone[k, 1], mrzone[k, 2] + 1):
                if htm[j] <= nml.tcut[1]:
                    done = True  # GOTO 190
                    break
                ibegin = j
                iend = j
                dburnm(ibegin, iend, m, hs2, hcomp, dt, dratem, dratem0,
                       fac2)
            if done:
                break
        # label 190. Convection zones, zone by zone.
        for k in range(1, nzone + 1):
            ibegin = mxzone[k, 1]
            iend = mxzone[k, 2]
            dburnm(ibegin, iend, m, hs2, hcomp, dt, dratem, dratem0, fac2)
    # Extent of a central convection zone.
    # Question (see the page): IMIN is the last convective shell.
    if lczm[1]:
        i = m + 1
        for k in range(2, m + 1):  # DO 60
            if not lczm[k]:
                i = k
                break
        imin = max(2, i - 1)  # label 65
    else:
        imin = 2
    # Extent of a surface convection zone.
    lallcz = False
    if lczm[m]:
        i = 0
        lallcz = True
        for k in range(m - 1, 0, -1):  # DO 70
            if not lczm[k]:
                i = k
                lallcz = False
                break
        imax = i + 1  # label 80
        # HSTOP, HSBOT: the masses at the top and bottom of the zone.
        hstop = math.exp(c.cln * hstot)
        if imax > 1:
            hsbot = 0.5 * (hs1[imax] + hs1[imax - 1])
        else:
            hsbot = 0.0
        # Question (see the page): true for any surface zone.
        if (hstop - hsbot) / c.cmsun > 0.0:
            lczsur = True
        else:
            lczsur = False
    else:
        imax = m
        lczsur = False
    # A local convective overturn time is not computed here (MHP 3/09).
    # IMPJMOD = 1 is the same as LSOLID (JNT 09/25).
    if not nml.lsolid and nml.impjmod != 1:
        # First guess: thin spherical shells, I/M = 2/3 r**2, Omega = J/I.
        for i in range(1, m + 1):
            him[i] = c.cc23 * math.exp(c.cln * 2.0 * hrm[i])
            omegam[i] = hjm[i] / him[i]
            him[i] = hs2[i] * him[i]
        getrot(hdm, hjm, hrm, hs, hs2, lczm, m, eta2m, him, omegam, qiwm,
               r0m)
        # Moment of inertia of the surface convection zone.
        if lczsur:
            hicz = him[m]
            for i in range(m - 1, 0, -1):
                if not lczm[i]:
                    break  # GOTO 110
                hicz = hicz + him[i]
    else:
        jstart = 1
        jend = m
        solid(hdm, hjm, hrm, hs, hs2, jstart, jend, eta2m, him, omegam,
              qiwm, r0m, m)
        # Question (see the page): the whole star for solid rotation.
        if lczsur:
            hicz = him[1]
            for i in range(2, m + 1):
                hicz = hicz + him[i]
    return (fx, hicz, hsbot, hstop, imin, imax, lallcz, lczsur, nrzone,
            nzone)


def _receding_cz(m, fx, hjm, hr, hs2, lcm):
    # Correct the specific angular momentum of shells that left a convection
    # zone during the step (MHP 06/02). Never run (LNEWCZ is reset first);
    # translated as written, its apparent bugs included (question 1 of the
    # page). Returns FX, which the Fortran overwrites.
    c = constants_mod
    hro = oldmod_mod.hro
    lcc = False
    ldo = False
    for i in range(1, m):
        if lchg[i]:
            if not lcc:
                lcc = True  # start of a changed region
                i0 = i
        if not lchg[i]:
            if lcc:
                i1 = i - 1  # end of a changed region
                lcc = False
                ldo = True
                print(i0, i1, ddeln[i0], ddeln[i0 - 1], ddeln[i1 + 1])
                # Only a receding convection zone changes J.
                if ddeln[i0] > 0.0:
                    if ddeln[i0 - 1] < 0.0:
                        ind = -1  # subtracted from a core zone
                    elif ddeln[i1 + 1] < 0.0:
                        ind = 1  # subtracted from an envelope zone
                    else:
                        ldo = False
                else:
                    ldo = False
        if ldo:
            ldo = False
            if ind < 0:
                # Core zone: locate its edges.
                icz1 = i0 - 1
                icz0 = 1
                for ii in range(i0 - 1, 0, -1):
                    if not lcm[ii]:
                        icz0 = ii + 1
                        break  # GOTO 11
                js = range(i1, i0 - 1, -1)
            else:
                # Envelope zone: locate its edges.
                icz0 = i1 + 1
                icz1 = m
                for ii in range(i1 + 2, m + 1):
                    if not lcm[ii]:
                        icz1 = ii - 1
                        break  # GOTO 12
                js = range(i0, i1 + 1)
            for j in js:
                if ddeln[j] < 0.0:
                    break  # GOTO 13
                czi = 0.0
                czj = 0.0
                # Fractions of the step the shell was radiative, convective.
                fx = ddeln[j] / (ddeln[j] - ddelp[j])
                f2 = 1.0 - fx
                zone = range(icz0, j) if ind < 0 else range(j + 1, icz1 + 1)
                for ii in zone:
                    hrr = hrp[ii] + f2 * (hr[ii] - hro[ii])
                    czi = czi + c.cc23 * hs2[ii] * 10.0**(2.0 * hrr)
                    czj = czj + hjm[ii] * hs2[ii]
                # The zone's J/M when the shell was released; the
                # difference goes back into the zone.
                hrsh = 2.0 * (hrp[j] + f2 * (hr[j] - hrp[j]))
                hjm2 = c.cc23 * (czj / czi) * 10.0**hrsh
                delj = (hjm[j] - hjm2) * hs2[j]
                ratio = (czj - delj) / czj
                hjm[j] = hjm2
                rr = hjm2 / hjm[j]
                print(j, fx, rr, delj, ratio)
                rescale = range(icz0, j) if ind < 0 else range(j + 1, m + 1)
                for ii in rescale:
                    hjm[ii] = hjm[ii] * ratio
            # label 13
    return fx
