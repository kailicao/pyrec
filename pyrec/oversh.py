"""Python version of oversh.f (YREC 09c3b9b): convective overshoot at the
edges of each convection zone."""

import math

from pyrec import constants_mod
from pyrec import dpmix_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.hsubp import hsubp


def _zone_lines(label, mx, nzone):
    """FORMAT(1X,label,4('[',I4,'-',I4,' ]')): four zones a line; further
    lines restart at the repeated group, without the label (question 2)."""
    groups = [f"[{mx[i, 1]:4d}-{mx[i, 2]:4d} ]" for i in range(1, nzone + 1)]
    lines = [" " + label + "".join(groups[:4])]
    for n in range(4, len(groups), 4):
        lines.append("".join(groups[n:n + 4]))
    return "\n".join(lines)


def oversh(hcomp, hd, hp, hr, hs, ht, m, mxzone, mxzon0, nzone):
    """Extend each convection zone by overshoot: by alpha_c H_P above a
    core (with LOVMAX, at most beta_c r), alpha_e H_P below the envelope,
    alpha_m H_P both ways for an intermediate zone. H_P is the pressure
    scale height at the zone's formal edge (HSUBP). Also sets IOV1, IOV2
    in DPMIX_MOD for the last zone overshooting downward.

    Translates: oversh.f#OVERSH (09c3b9b)

    Arguments: the model hcomp, hd, hp, hr, hs, ht (index 0 unused); m, the
    number of shells; mxzone[12+1, 2+1], the zones' edges, changed in place;
    mxzon0, the edges without overshoot; nzone. Returns nothing.
    """
    nml = nmlpar_mod
    cln = constants_mod.cln
    ishort = iounits_mod.ishort
    pscalu = 0.0
    pscald = 0.0
    # Question (see the page): IOV1, IOV2 keep only the last zone.
    dpmix_mod.iov1 = -1
    dpmix_mod.iov2 = -1
    for i in range(1, nzone + 1):
        if mxzone[i, 1] == 1:
            if mxzone[i, 2] == m:
                iounits_mod.write(ishort,
                                  " FULLY CONVECTIVE MODEL - NO OVERSHOOT")
                # Question (see the page): returns before the second line.
                iounits_mod.write(ishort, _zone_lines(
                    "MIXED REGIONS WITHOUT OVERSHOOT", mxzon0, 1))
                return
            if not nml.lovstc:
                continue
            lflagu = True
            lflagd = False
            jedge = mxzone[i, 2]
            pscalu = hsubp(hcomp, hd, hp, hr, hs, ht, jedge, pscalu)
            if nml.lovmax:
                pscalu = min(pscalu * nml.alphac,
                             nml.betac * math.exp(cln * hr[jedge]))
            else:
                pscalu = pscalu * nml.alphac
        elif mxzone[i, 2] == m:
            if not nml.lovste:
                continue
            lflagu = False
            lflagd = True
            jedge = mxzone[i, 1]
            pscald = hsubp(hcomp, hd, hp, hr, hs, ht, jedge, pscald)
            pscald = pscald * nml.alphae
        else:
            if not nml.lovstm:
                continue
            lflagu = True
            lflagd = True
            jedge = mxzone[i, 1]
            pscald = hsubp(hcomp, hd, hp, hr, hs, ht, jedge, pscald)
            pscald = pscald * nml.alpham
            jedge = mxzone[i, 2]
            pscalu = hsubp(hcomp, hd, hp, hr, hs, ht, jedge, pscalu)
            pscalu = pscalu * nml.alpham
        if lflagd:
            # Downward: the innermost point with r >= r_edge - d.
            id_ = mxzone[i, 1]
            rcz = math.exp(cln * hr[id_])
            ros = rcz - pscald
            for j in range(id_ - 1, 0, -1):
                r = math.exp(cln * hr[j])
                if r < ros:
                    break       # goto 20
            else:
                j = 0
            # label 20
            mxzone[i, 1] = j + 1
            if nml.lrot and nml.linstb and mxzone[i, 1] == mxzon0[i, 1]:
                mxzone[i, 1] = mxzone[i, 1] - 1
            dpmix_mod.iov2 = id_
            dpmix_mod.iov1 = mxzone[i, 1]
        if lflagu:
            # Upward: the outermost point with r <= r_edge + d.
            id_ = mxzone[i, 2]
            rcz = math.exp(cln * hr[id_])
            ros = rcz + pscalu
            for j in range(id_ + 1, m + 1):
                r = math.exp(cln * hr[j])
                if r > ros:
                    break       # goto 40
            else:
                j = m + 1
            # label 40
            mxzone[i, 2] = j - 1
            if nml.lrot and nml.linstb and mxzone[i, 2] == mxzon0[i, 2]:
                mxzone[i, 2] = mxzone[i, 2] + 1
    # label 100 ends the loop
    iounits_mod.write(ishort, _zone_lines(
        "MIXED REGIONS WITHOUT OVERSHOOT", mxzon0, nzone))
    iounits_mod.write(ishort, _zone_lines(
        "MIXED REGIONS WITH OVERSHOOT   ", mxzone, nzone))
