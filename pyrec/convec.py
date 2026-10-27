"""Python version of convec.f (YREC 09c3b9b): the convection zones (without
and with overshoot), the radiative regions, JCORE and JENV."""

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.oversh import oversh


def convec(hcomp, hd, hp, hr, hs, ht, lc, m, mrzone, mxzone, mxzon0,
           jcore, jenv, nrzone, nzone, nzone0):
    """Find the convection zones, without overshoot (mxzon0, nzone0) and
    with it (mxzone, nzone, merging zones that overlap), the radiative
    regions between them (mrzone, nrzone), the edge of the convective
    core jcore and the base of the convective envelope jenv (both before
    overshoot). Overshoot is overmixing: the overshoot shells are mixed
    with the zone but keep the radiative gradient.

    Translates: convec.f#CONVEC (09c3b9b)

    Arguments: the model hcomp, hd, hp, hr, hs, ht (used only by oversh);
    lc[1..m], convective or not (lc[m+1] is set to False here); m, the
    number of shells; the zone arrays mrzone[k, 1..2], mxzone[k, 1..2],
    mxzon0[k, 1..2] (changed in place); jcore, jenv, nrzone, nzone,
    nzone0, passed in (with any value).

    Returns (jcore, jenv, nrzone, nzone, nzone0).
    """
    n = nmlpar_mod

    # The standard convection zones; a run of a single shell is ignored.
    # Question (see the page): one-shell zones treated as radiative.
    j = 1
    lcc = False
    lc[m + 1] = False
    for i in range(1, m + 2):
        if lc[i]:
            if not lcc:
                # The start of a convection zone.
                lcc = True
                i1 = i
            continue
        # label 10
        if not lcc:
            continue
        # The end of a convection zone.
        lcc = False
        if i1 != i - 1:
            mxzone[j, 1] = i1
            mxzone[j, 2] = i - 1
            j = j + 1
        if j < 12:
            continue
        iounits_mod.write(iounits_mod.ishort, ' -----TOO MANY MIX ZONES')
        break  # GOTO 12
    # label 12

    # MHP 5/91: no convection zone.
    if j == 1:
        jcore = 1
        jenv = m
        nzone0 = 0
        nzone = 0
        nrzone = 1
        mrzone[1, 1] = 1
        mrzone[1, 2] = m
        return jcore, jenv, nrzone, nzone, nzone0  # GOTO 9999
    nzone = j - 1
    for i in range(1, nzone + 1):
        mxzon0[i, 1] = mxzone[i, 1]
        mxzon0[i, 2] = mxzone[i, 2]
    nzone0 = nzone

    # The outer edge of the convective core (JCORE) and the inner edge of
    # the convective envelope (JENV), before overshoot.
    if mxzone[1, 1] == 1:
        if mxzone[1, 2] == m:
            # A fully convective star (MRZONE is left as it was).
            jcore = 1
            jenv = 1
            nrzone = 0
            return jcore, jenv, nrzone, nzone, nzone0  # GOTO 9999
        else:
            jcore = mxzone[1, 2]
    else:
        jcore = 1
    if mxzone[nzone, 2] == m:
        jenv = mxzone[nzone, 1]
    else:
        jenv = m

    # Overshoot, if any switch is on: oversh moves the edges in mxzone;
    # then merge zones whose edges now meet.
    if n.lovstc or n.lovste or n.lovstm:
        oversh(hcomp, hd, hp, hr, hs, ht, m, mxzone, mxzon0, nzone)
        if nzone != 1:
            j = 1
            while True:  # label 85
                if mxzone[j, 2] >= mxzone[j + 1, 1]:
                    # Two convection zones have merged.
                    iounits_mod.write(
                        iounits_mod.ishort,
                        '  CONVECTION ZONES MERGED DUE TO OVERSHOOT\n'
                        f'  OLD[{mxzone[j, 1]:3d}-{mxzone[j, 2]:3d}]'
                        f'[{mxzone[j + 1, 1]:3d}-{mxzone[j + 1, 2]:3d}]'
                        f' NEW[{mxzone[j, 1]:3d}-{mxzone[j + 1, 2]:3d}]')
                    mxzone[j + 1, 1] = mxzone[j, 1]
                    for k in range(j, nzone):
                        for j2 in range(1, 3):
                            mxzone[k, j2] = mxzone[k + 1, j2]
                    nzone = nzone - 1
                    if j <= nzone - 1:
                        continue  # GOTO 85
                    break  # GOTO 100
                j = j + 1
                if j <= nzone - 1:
                    continue  # GOTO 85
                break
    # label 100: the radiative regions.
    # Below the first convection zone.
    if mxzone[1, 1] > 1:
        nrzone = 1
        mrzone[1, 1] = 1
        mrzone[1, 2] = mxzone[1, 1] - 1
    else:
        nrzone = 0
    # Between convection zones.
    if nzone > 1:
        for i in range(1, nzone):
            nrzone = nrzone + 1
            mrzone[nrzone, 1] = mxzone[i, 2] + 1
            mrzone[nrzone, 2] = mxzone[i + 1, 1] - 1
    # Above the last convection zone.
    if mxzone[nzone, 2] < m:
        nrzone = nrzone + 1
        mrzone[nrzone, 1] = mxzone[nzone, 2] + 1
        # Question (see the page): MRZONE(1,2), not MRZONE(NRZONE,2).
        mrzone[1, 2] = m
    # label 9999
    return jcore, jenv, nrzone, nzone, nzone0
