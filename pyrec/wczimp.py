"""Python version of wczimp.f (YREC 09c3b9b): the rotation of a convection
zone at fixed angular momentum, with options to impose a rotation profile
on the interior below the surface convection zone."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec.solid import solid


def wczimp(hd, hjm, hr, hs, hs2, istart, iend, eta2, hi, omega, qiw, r0,
           m):
    """Set the rotation of the convection zone istart..iend at fixed total
    angular momentum: solid body, constant specific angular momentum, or
    Omega = W0 * r**WALPCZ; IMPJMOD = 2 or 3 also makes the interior below
    the surface zone a solid body, decoupled from or coupled to the zone.

    Translates: wczimp.f#WCZIMP (09c3b9b)

    Arguments: the model's arrays hd, hjm (specific angular momentum), hr
    (log radius), hs, hs2 (shell masses), eta2, hi (moments of inertia),
    omega, qiw, r0; the zone's first and last shells istart, iend; the
    number of shells m. hjm, omega, hi, qiw, eta2 and r0 are updated in
    place (through solid). Nothing is returned.
    """
    n = nmlpar_mod
    cln = constants_mod.cln
    walpcz = n.walpcz

    # JVS 2026/02/13: GT reverted to GE.
    # Question (see the page): IMPJMOD = 2 and 3 are reached only when
    # WALPCZ < 0, since this test comes first.
    if walpcz >= 0.0 or n.lsolid or n.impjmod == 1 or iend < m:
        # Solid-body rotation in convective regions.
        solid(hd, hjm, hr, hs, hs2, istart, iend, eta2, hi, omega, qiw,
              r0, m)

    elif n.impjmod == 2:
        # The interior below the surface zone rotates as a solid body,
        # detached from the convection zone.
        cjtot = hjm[istart] * hs2[istart]
        cmtot = math.exp(cln * walpcz * hr[istart]) * hi[istart]
        for j in range(istart + 1, iend + 1):  # label 30
            cjtot = cjtot + hjm[j] * hs2[j]
            cmtot = cmtot + math.exp(cln * walpcz * hr[j]) * hi[j]
        w0 = cjtot / cmtot
        for j in range(istart, iend + 1):  # label 31
            omega[j] = w0 * math.exp(cln * walpcz * hr[j])
            hjm[j] = omega[j] * hi[j] / hs2[j]
            solid(hd, hjm, hr, hs, hs2, j, j, eta2, hi, omega, qiw, r0, m)
        if istart > 1:
            solid(hd, hjm, hr, hs, hs2, 1, istart - 1, eta2, hi, omega,
                  qiw, r0, m)

    elif n.impjmod == 3:
        # The interior rotates as a solid body at the rate of the zone's
        # base, sharing one angular-momentum budget with the zone.
        rbase = math.exp(cln * walpcz * hr[istart])
        cjtot = hjm[1] * hs2[1]
        cmtot = rbase * hi[1]
        if istart > 1:
            for j in range(2, istart + 1):  # label 32
                cjtot = cjtot + hjm[j] * hs2[j]
                cmtot = cmtot + rbase * hi[j]
        for j in range(istart + 1, iend + 1):  # label 33
            cjtot = cjtot + hjm[j] * hs2[j]
            cmtot = cmtot + math.exp(cln * walpcz * hr[j]) * hi[j]
        w0 = cjtot / cmtot
        for j in range(istart, iend + 1):  # label 34
            omega[j] = w0 * math.exp(cln * walpcz * hr[j])
            hjm[j] = omega[j] * hi[j] / hs2[j]
            solid(hd, hjm, hr, hs, hs2, j, j, eta2, hi, omega, qiw, r0, m)
        if istart > 1:
            for j in range(1, istart):  # label 35
                omega[j] = omega[istart]
                hjm[j] = omega[j] * hi[j] / hs2[j]
            # Question (see the page): SOLID may change the interior's
            # Omega again.
            solid(hd, hjm, hr, hs, hs2, 1, istart - 1, eta2, hi, omega,
                  qiw, r0, m)

    elif walpcz <= -2.0 and n.impjmod == 0:
        # Constant specific angular momentum in the convection zone.
        cjtot = hjm[istart] * hs2[istart]
        cmtot = hs2[istart]
        for j in range(istart + 1, iend + 1):  # label 36
            cjtot = cjtot + hjm[j] * hs2[j]
            cmtot = cmtot + hs2[j]
        cjnew = cjtot / cmtot
        for j in range(istart, iend + 1):  # label 40
            hjm[j] = cjnew
            omega[j] = cjnew * hs2[j] / hi[j]
            solid(hd, hjm, hr, hs, hs2, j, j, eta2, hi, omega, qiw, r0, m)

    else:
        # General law Omega = W0 * r**WALPCZ in the convection zone, W0 set
        # by conserving the zone's total angular momentum.
        cjtot = hjm[istart] * hs2[istart]
        cmtot = math.exp(cln * walpcz * hr[istart]) * hi[istart]
        for j in range(istart + 1, iend + 1):  # label 50
            cjtot = cjtot + hjm[j] * hs2[j]
            cmtot = cmtot + math.exp(cln * walpcz * hr[j]) * hi[j]
        w0 = cjtot / cmtot
        for j in range(istart, iend + 1):  # label 60
            omega[j] = w0 * math.exp(cln * walpcz * hr[j])
            hjm[j] = omega[j] * hi[j] / hs2[j]
            solid(hd, hjm, hr, hs, hs2, j, j, eta2, hi, omega, qiw, r0, m)
