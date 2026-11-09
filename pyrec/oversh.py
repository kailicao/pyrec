"""Python version of oversh.f (YREC 7be4133): convective overshoot at the
edges of each convection zone."""

import math
import numpy as np

from pyrec import constants_mod
from pyrec import dpmix_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.hsubp import hsubp


def _zone_lines(label: str, mx: np.ndarray, nzone: int) -> str:
    """FORMAT(1X,label,4('[',I4,'-',I4,' ]')): four zones a line; further
    lines restart at the repeated group, without the label (Fortran's
    format reversion; Fortran lines 170–176)."""

    groups = [f"[{mx[i, 1]:4d}-{mx[i, 2]:4d} ]" for i in range(1, nzone + 1)]
    lines = [" " + label + "".join(groups[:4])]
    for n in range(4, len(groups), 4):
        lines.append("".join(groups[n:n + 4]))
    return "\n".join(lines)


def oversh(hcomp: np.ndarray, hd: np.ndarray, hp: np.ndarray, hr: np.ndarray, hs: np.ndarray, ht: np.ndarray, m: int,
           mxzone: np.ndarray, mxzon0: np.ndarray, nzone: int) -> None:
    """Extend each convection zone by overshoot: by alpha_c H_P above a
    core (with LOVMAX, at most beta_c r), alpha_e H_P below the envelope,
    alpha_m H_P both ways for an intermediate zone. H_P is the pressure
    scale height at the zone's formal edge (HSUBP). Also sets IOV1, IOV2
    in DPMIX_MOD for the last zone overshooting downward.

    Translates: oversh.f#OVERSH (7be4133)

    Arguments: the model hcomp, hd, hp, hr, hs, ht (index 0 unused); m, the
    number of shells; mxzone[12+1, 2+1], the zones' edges, changed in place;
    mxzon0, the edges without overshoot; nzone. Returns nothing.

    """

    nml = nmlpar_mod
    cln = constants_mod.cln
    ishort = iounits_mod.ishort
    pscalu = 0.0
    pscald = 0.0

    # Fortran lines 61–62.
    # Question (see the page): IOV1, IOV2 keep only the last zone.
    dpmix_mod.iov1 = -1
    dpmix_mod.iov2 = -1

    for i in range(1, nzone + 1):
        # The zone's type and the overshoot distances (Fortran lines 64–119).
        if mxzone[i, 1] == 1:
            if mxzone[i, 2] == m:
                iounits_mod.write(ishort,
                                  " FULLY CONVECTIVE MODEL - NO OVERSHOOT")
                # Never reached: CONVEC returns before calling OVERSH for a
                # fully convective star (see the page).
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

        # Downward: the innermost point with r >= r_edge - d (Fortran lines
        # 120–146).
        if lflagd:
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

        # Upward: the outermost point with r <= r_edge + d (Fortran lines
        # 147–168).
        if lflagu:
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

    # label 100 ends the loop. The zones without and with overshoot
    # (Fortran lines 170–176).
    iounits_mod.write(ishort, _zone_lines(
        "MIXED REGIONS WITHOUT OVERSHOOT", mxzon0, nzone))
    iounits_mod.write(ishort, _zone_lines(
        "MIXED REGIONS WITH OVERSHOOT   ", mxzone, nzone))
