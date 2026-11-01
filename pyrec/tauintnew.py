"""Python version of tauintnew.f (YREC 09c3b9b): the convective overturn
time by G. Somers' method (3/17)."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import ovrtrn_mod
from pyrec.kspline import kspline
from pyrec.ksplint import ksplint

# TAUINTNEW has SAVE. Of its locals only JCORE can be read before it is set
# in a call, so it alone is kept between calls here.
jcore = 0


def _psca(hp, hd, hg, k):
    # The pressure scale height P/(rho g) at shell k.
    cln = constants_mod.cln
    return math.exp(cln * (hp[k] - hd[k])) / hg[k]


def tauintnew(hs1, lcz, hr, hp, hd, hg, m, mint, svel, del1, del2, rbcz):
    """Compute the convective overturn time TAUCZ (in ovrtrn_mod) and the
    radius of the base of the convection zone: find the zone from the
    centre up and from the surface down; for a single zone, H_P/v one
    pressure scale height above its base; otherwise the lowest cell's width
    over the velocity at its middle.

    Translates: tauintnew.f#TAUINTNEW (09c3b9b)

    Arguments: the profile of m points (mint of them interior shells), index
    0 unused: hs1 (mass, g), lcz (convective or not), hr, hp, hd (log r, P,
    rho), hg (gravity), svel (convective velocity), del1, del2 (the two
    gradients compared to find the zone's edges); rbcz is output only.

    Returns rbcz, the base's radius in solar radii (0 for a radiative
    envelope).
    """
    global jcore
    cln = constants_mod.cln
    xdel = np.zeros(5)
    xrad = np.zeros(5)
    yrad = np.zeros(5)
    yvel = np.zeros(5)
    ypsc = np.zeros(5)
    yspl = np.zeros(5)

    lallcz = False
    ibeg = 1
    while True:  # label 50
        # Find the top of the main internal convection zone.
        # Question (see the page): jcore keeps its value from an earlier pass
        # or call when the search starts in a radiative shell.
        # Question (see the page): nothing stops ibeg from passing m.
        lradzone = not lcz[ibeg]
        found = False
        for i in range(ibeg, m + 1):  # loop 60
            if not lradzone and not lcz[i]:
                jcore = i - 1
                lradzone = True
            if lradzone and lcz[i]:
                found = True
                break  # to label 70
        if not found:
            # Fully convective, or a radiative envelope (maybe with a
            # convective core).
            if lradzone:
                lallcz = False
                lczsur = False
                imax = m
                jmax = m
                jcore = m
                to_100 = True
            else:
                lallcz = True
                lczsur = True
                i = 1
                j = 0
                jcore = 1
                to_100 = False
        else:
            # label 70: find the base of the uppermost convection zone.
            # LCZ(J+1) traps the case where the last points are already
            # radiative (JvS 02/2026).
            j = 0  # the Fortran's J after the loop runs out
            for jj in range(m - 1, 0, -1):  # loop 80
                if lcz[jj + 1] and not lcz[jj]:
                    j = jj
                    break  # to label 90
            to_100 = False

        if not to_100:
            # label 90
            imax = i
            jmax = j + 1
            ijmax = jmax - imax
            # Skip convective shells just above the convective core.
            # The check on imax and jcore (KC 2025-05-30).
            # Question (see the page): for a fully convective star with a
            # stitched envelope, this restarts the search for ever.
            if imax > 0 and jcore > 0:
                if m != mint and hp[jcore] - hp[imax] < 1:
                    ibeg = imax + 1
                    continue  # GOTO 50
            # At least 3 shells thick, or search again from above.
            # Question (see the page): lcz[imax + 2] can be beyond m.
            if not (lcz[imax] and lcz[imax + 1] and lcz[imax + 2]):
                if not lcz[imax + 1]:
                    ibeg = imax + 1
                else:
                    ibeg = imax + 2
                continue  # GOTO 50
            # Interior model only: use jmax.
            if ijmax != 0 and m == mint:
                imax = jmax
                ijmax = 0
            # Masses at the top and bottom of the zone (unused below).
            hstop = hs1[m]
            if jmax > 1:
                hsbot = 0.5 * (hs1[jmax] + hs1[jmax - 1])
            else:
                hsbot = 0.0
            lczsur = (hstop - hsbot) / constants_mod.cmsun > 0.0

        # label 100: taucz from the top down; 0 for a radiative envelope.
        if lczsur:
            if lallcz:
                # Question (see the page): exp(cln * 1) is 10 cm, not 1 cm.
                rbcz = 1.0
            else:
                # The base where del1 - del2 changes sign, linear in log r.
                dd2 = del1[imax - 1] - del2[imax - 1]
                dd1 = del1[imax] - del2[imax]
                fx = dd2 / (dd2 - dd1)
                rbcz = hr[imax - 1] + fx * (hr[imax] - hr[imax - 1])
            rbcz = math.exp(cln * rbcz)
            to_110 = ijmax != 0  # several cells: the averaging method
            if not to_110:
                psca2 = _psca(hp, hd, hg, imax)
                rtest2 = math.exp(cln * hr[imax]) - rbcz
                valtest = rtest2 - psca2
                if valtest > 0.0:
                    # The innermost point's scale height reaches below.
                    k = imax
                    cvel = svel[imax]
                    psca = psca2
                    ovrtrn_mod.taucz = psca / cvel
                    # Possible bug (see the page): with the ENDIF where it now
                    # is, this branch falls through to label 110, which
                    # overwrites taucz (the ENDIF moved, KC 2025-05-31).
                    to_110 = True
                else:
                    to_110 = True
                    for k in range(imax + 1, m):
                        psca1 = psca2
                        rtest1 = rtest2
                        psca2 = _psca(hp, hd, hg, k)
                        rtest2 = math.exp(cln * hr[k]) - rbcz
                        valtest = rtest2 - psca2
                        # Where r - r_b = H_P.
                        if valtest > 0.0:
                            kk = k
                            if k == imax + 1:
                                kk = k + 1
                            # Question (see the page): kk + 1 can reach m + 1.
                            for n in range(1, 5):
                                kn = kk - 3 + n
                                xdel[n] = (math.exp(cln * hr[kn]) - rbcz
                                           - _psca(hp, hd, hg, kn))
                            for n in range(1, 5):
                                yvel[n] = svel[kk - 3 + n]
                            kspline(xdel, yvel, yspl)
                            cvel = ksplint(xdel, yvel, yspl, 0.0, 0.0)
                            for n in range(1, 5):
                                ypsc[n] = _psca(hp, hd, hg, kk - 3 + n)
                            kspline(xdel, ypsc, yspl)
                            psca = ksplint(xdel, ypsc, yspl, 0.0, 0.0)
                            for n in range(1, 5):
                                yrad[n] = hr[kk - 3 + n]
                            kspline(xdel, yrad, yspl)
                            # Question (see the page): envrcl is not used.
                            envrcl = ksplint(xdel, yrad, yspl, 0.0, 0.0)
                            ovrtrn_mod.taucz = psca / cvel
                            to_110 = False
                            break  # to label 140
            if to_110:
                # label 110: r_b + H_P is outside the star, or several cells
                # lie in the radiative envelope; use the lowest cell's
                # middle. First find its upper edge.
                k = m  # the Fortran's K after the loop runs out
                for kk2 in range(imax, m):
                    if not lcz[kk2]:
                        k = kk2
                        break  # to label 120
                itop = k - 1  # label 120
                dd2 = del1[itop] - del2[itop]
                dd1 = del1[itop + 1] - del2[itop + 1]
                fx = dd2 / (dd2 - dd1)
                rtcz = hr[itop] + fx * (hr[itop + 1] - hr[itop])
                rtcz = math.exp(cln * rtcz)
                delr = rtcz - rbcz
                # The shells around the cell's middle.
                rtest = math.log10(rbcz + 0.5 * delr)
                k = itop + 2  # the Fortran's K after the loop runs out
                for kk2 in range(imax, itop + 2):
                    if hr[kk2] > rtest:
                        k = kk2
                        break  # to label 130
                rtest = math.exp(cln * rtest)  # label 130
                kk = k
                if k == imax + 1:
                    kk = k + 1
                for n in range(1, 5):
                    xrad[n] = math.exp(cln * hr[kk - 3 + n])
                    yvel[n] = svel[kk - 3 + n]
                kspline(xrad, yvel, yspl)
                cvel = ksplint(xrad, yvel, yspl, rtest, 0.0)
                ovrtrn_mod.taucz = delr / cvel
            # label 140: the base's radius in solar units.
            rbcz = rbcz / constants_mod.crsun
        else:
            ovrtrn_mod.taucz = 0.0
            rbcz = 0.0
        # A taucz found deep in the interior: search again above.
        if ovrtrn_mod.taucz > 1.0e20:
            ibeg = imax + 1
            continue  # GOTO 50
        return rbcz
