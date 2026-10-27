"""Python version of taucal.f (YREC 09c3b9b): the convective overturn time
near the base of a surface convection zone lying within the envelope."""

import math

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod

# Locals kept between calls (SAVE). Question (see the page): IMAX and
# ICZTOP may be used without being set in this call; 0 stands for the
# compiler's usual zero initial value.
imax = 0
icztop = 0


def taucal(hs2, hs1, lcz, hr, hp, hd, hg, m, svel, del1, del2):
    """Set ovrtrn_mod.taucz, the convective overturn time near the base of
    the surface convection zone of the envelope: H_P / v one pressure scale
    height above the base or, for a thin zone or two zones, the zone's
    width over its mass-weighted mean velocity.

    Translates: taucal.f#TAUCAL (09c3b9b)

    Arguments: hs2 (shell masses), hs1 (masses), lcz (convective flags),
    hr, hp, hd (log r, log P, log rho), hg (gravity; changed in place at
    the base when there is no rotation), m, svel (convective velocities),
    del1, del2 (two temperature gradients). None is returned.
    """
    global imax, icztop
    c = constants_mod
    lallcz = False
    lflag = False
    ltwocz = False
    lczex = False  # a convection zone exists
    # Down from the surface: the first convective shell, the zone's base
    # IMAX, and a second, deeper zone (top ICZTOP).
    for i in range(m - 1, 0, -1):
        if not lcz[i] and not lflag and lczex:  # leaving the zone
            imax = i + 1
            lflag = True
        elif lflag and lcz[i]:  # two zones with a radiative one between
            lflag = False
            ltwocz = True
            icztop = i
        elif lcz[i] and not lczex:  # first convective shell from the top
            lczex = True
    if imax < 1:
        lallcz = True
    if lczex:
        hstop = hs1[m]
        if imax > 1:
            hsbot = 0.5 * (hs1[imax] + hs1[imax - 1])
        else:
            hsbot = hs1[1]
        # Question (see the page): LCZSUR is true for any zone and unused.
        if (hstop - hsbot) / c.cmsun > 0.0:
            lczsur = True
        else:
            lczsur = False
    else:
        # Question (see the page): no zone, but a turnover time below.
        imax = m
        lczsur = False

    if not lallcz:
        if not nmlpar_mod.lrot:
            hg[imax] = hs1[imax] * math.exp(c.cln * (c.cgl - 2.0 * hr[imax]))
            hg[imax - 1] = hs1[imax - 1] * math.exp(
                c.cln * (c.cgl - 2.0 * hr[imax - 1]))
        # The base r_b, where DEL1 - DEL2 changes sign, and H_P there.
        dd2 = del1[imax - 1] - del2[imax - 1]
        dd1 = del1[imax] - del2[imax]
        fx = dd2 / (dd2 - dd1)
        envrl = hr[imax - 1] + fx * (hr[imax] - hr[imax - 1])
        envr = math.exp(c.cln * envrl)
        ps2 = math.exp(c.cln * (hp[imax] - hd[imax])) / hg[imax]
        ps1 = math.exp(c.cln * (hp[imax - 1] - hd[imax - 1])) / hg[imax - 1]
        psca = ps1 + fx * (ps2 - ps1)
        rtestl = math.log10(envr + psca)
        tol = 0.5
        # Question (see the page): .OR., not .AND.; ICZTOP may be unset.
        if (abs(1.0 - ((math.exp(c.cln * (rtestl - hr[imax])) - 1.0)
                       / (math.exp(c.cln * (hr[icztop] - hr[imax])) - 1.0)))
                < tol or ltwocz):
            # The mass-weighted mean convective velocity.
            # Question (see the page): the width from the shells.
            denom = 0.0
            umer = 0.0
            czwide = abs(math.exp(c.cln * hr[icztop])
                         - math.exp(c.cln * hr[imax]))
            for j in range(imax, icztop + 1):
                denom = denom + hs2[j]
                umer = umer + 0.5 * (svel[j] + svel[j + 1]) * hs2[j]
            cveavg = umer / denom
            ovrtrn_mod.taucz = czwide / cveavg
        else:
            # The velocity one scale height above the base.
            for k in range(imax + 1, m + 1):
                if hr[k] > rtestl:
                    fx = (rtestl - hr[k - 1]) / (hr[k] - hr[k - 1])
                    cvel = svel[k - 1] + fx * (svel[k] - svel[k - 1])
                    break  # to label 85
            else:
                cvel = svel[m]
            ovrtrn_mod.taucz = psca / cvel  # label 85
    else:
        ovrtrn_mod.taucz = 0.0
    # Question (see the page): a debugging print at every call.
    print(" TauCal")
