"""Python version of tauint.f (YREC 09c3b9b): the convective overturn time
at the base of the surface convection zone, on the interior profile (the
older method; TAUINTNEW is the newer)."""

import math

from pyrec import constants_mod
from pyrec import envstruct_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod


def tauint(hs1, lcz, hr, hp, hd, hg, m, svel, del1, del2):
    """Set ovrtrn_mod.taucz: H_P / v one pressure scale height above the
    base of the surface convection zone (searching the envelope of the
    last ENVINT call when the profile ends first); for a fully convective
    star, H_P / v where H_P = r; 0 with no surface zone.

    Translates: tauint.f#TAUINT (09c3b9b)

    Arguments: hs1 (masses), lcz (convective flags), hr, hp, hd (log r,
    log P, log rho), hg (gravity; changed in place where computed, when
    there is no rotation), m, svel (convective velocities), del1, del2
    (the radiative and adiabatic gradients). None is returned.
    """
    c = constants_mod
    env = envstruct_mod
    lrot = nmlpar_mod.lrot
    lallcz = False
    if lcz[m] or lcz[m - 1]:
        # IMAX: the lowest shell of the convective run from the surface.
        for i in range(m - 1, 0, -1):
            if not lcz[i]:
                break  # to label 81
        else:
            lallcz = True
            i = 0
        imax = i + 1  # label 81
        hstop = hs1[m]
        if imax > 1:
            hsbot = 0.5 * (hs1[imax] + hs1[imax - 1])
        else:
            hsbot = 0.0
        # Question (see the page): true for any zone.
        if (hstop - hsbot) / c.cmsun > 0.0:
            lczsur = True
        else:
            lczsur = False
    else:
        imax = m
        lczsur = False

    if not lczsur:
        ovrtrn_mod.taucz = 0.0
        return

    if not lallcz:
        # A partly convective star: the base r_b and H_P there.
        if not lrot:
            hg[imax] = hs1[imax] * math.exp(c.cln * (c.cgl - 2.0 * hr[imax]))
            hg[imax - 1] = hs1[imax - 1] * math.exp(
                c.cln * (c.cgl - 2.0 * hr[imax - 1]))
        dd2 = del1[imax - 1] - del2[imax - 1]
        dd1 = del1[imax] - del2[imax]
        fx = dd2 / (dd2 - dd1)
        envrcl = hr[imax - 1] + fx * (hr[imax] - hr[imax - 1])
        envrc = math.exp(c.cln * envrcl)
        ps2 = math.exp(c.cln * (hp[imax] - hd[imax])) / hg[imax]
        ps1 = math.exp(c.cln * (hp[imax - 1] - hd[imax - 1])) / hg[imax - 1]
        psca = ps1 + fx * (ps2 - ps1)
        rtestl = math.log10(envrc + psca)
        found = False
        for k in range(imax + 1, m + 1):
            if hr[k] > rtestl:
                fx = (rtestl - hr[k - 1]) / (hr[k] - hr[k - 1])
                cvel = svel[k - 1] + fx * (svel[k] - svel[k - 1])
                found = True
                break  # to label 85
        if not found:
            # One pressure scale height overshoots the interior: stitch on
            # the envelope for more room.
            # Question (see the page): between shell M and envelope point
            # K, skipping the points before K.
            for k in range(2, env.numenv + 1):
                if env.envr[k] > rtestl and env.evels[k] > 0.0:
                    fx = (rtestl - hr[m]) / (env.envr[k] - hr[m])
                    cvel = svel[m] + fx * (env.evels[k] - svel[m])
                    break  # to label 85
            else:
                cvel = svel[m]
        ovrtrn_mod.taucz = psca / cvel  # label 85
        return

    # A fully convective star: where H_P = r.
    if not lrot:
        hg[1] = hs1[1] * math.exp(c.cln * (c.cgl - 2.0 * hr[1]))
    psca2 = math.exp(c.cln * (hp[1] - hd[1])) / hg[1]
    rtest2 = math.exp(c.cln * hr[1])
    if psca2 <= rtest2:
        cvel = svel[1]
        psca = (psca2 * rtest2) ** 0.5
        ovrtrn_mod.taucz = psca / cvel
        return
    for k in range(2, m + 1):
        psca1 = psca2
        rtest1 = rtest2
        if not lrot:
            hg[k] = hs1[k] * math.exp(c.cln * (c.cgl - 2.0 * hr[k]))
        psca2 = math.exp(c.cln * (hp[k] - hd[k])) / hg[k]
        rtest2 = math.exp(c.cln * hr[k])
        if psca2 <= rtest2:
            fx = (rtest1 - psca1) / ((psca2 - rtest2) - (psca1 - rtest1))
            cvel = svel[k - 1] + fx * (svel[k] - svel[k - 1])
            psca = psca1 + fx * (psca2 - psca1)
            ovrtrn_mod.taucz = psca / cvel
            return  # to label 95
    cvel = svel[m]
    psca = psca2
    ovrtrn_mod.taucz = psca / cvel
