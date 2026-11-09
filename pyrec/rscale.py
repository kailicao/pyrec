"""Python version of rscale.f (YREC 7be4133): rescaling a starting model
as a kind card asks (X, Z, one species, the total mass, the core mass of a
horizontal-branch star), and the ZRAMP option."""

import math
import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import label_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod
from pyrec.findsh import findsh


def _error_stop(text_lines: list[str]) -> None:
    """Write an error message to ISHORT and stop the run."""

    for line in text_lines:
        iounits_mod.write(iounits_mod.ishort, line)
    raise SystemExit(1)


def rscale(hl: np.ndarray, hcomp: np.ndarray, hs: np.ndarray, hstot: float, m: int, nk: int, smass: float,
           lc: np.ndarray) -> tuple[float, float]:
    """Rescale the model for kind card nk: X, Z, one species (main
    sequence), the total mass (main sequence: every mass point scaled;
    horizontal branch: the mass coordinate stretched above the hydrogen
    shell or its base) and the core mass (horizontal branch); then reset Y
    and apply ZRAMP. The targets are runloop_mod.rescal[1..4, nk], the
    species nmlpar_mod.xnewcp, inewcp, and the ramp nmlpar_mod.rsclzc,
    rsclzm1, rsclzm2.

    Translates: rscale.f#RSCALE (7be4133)

    Arguments: hl, hcomp, hs (changed in place), hstot (log total mass),
    m, nk, smass (total mass in solar masses), lc. Also sets comp_mod.xnew,
    znew, zenv, stotal and label_mod.xenv0, zenv0. Returns (hstot, smass).

    """

    c = constants_mod
    comp = comp_mod
    nml = nmlpar_mod
    rescal = runloop_mod.rescal
    w = iounits_mod.write
    ishort = iounits_mod.ishort

    # X: every shell scaled by X_new / X_surf, capped at X_new (Fortran lines 91–117).
    if rescal[2, nk] > 0.0:
        if rescal[2, nk] <= 1.0:
            xfact = rescal[2, nk] / max(comp.xnew, 1.0e-20)
            for i in range(1, m + 1):
                hcomp[1, i] = min(rescal[2, nk], hcomp[1, i] * xfact)
            comp.xnew = rescal[2, nk]
            label_mod.xenv0 = rescal[2, nk]
        else:
            i1 = 2
            w(ishort, " ERROR IN SUBROUTINE RSCALE")  # format 1000
            w(ishort, f" RESCALING OF X IN KIND CARD #{nk:3d} FAILED - "
                      f"DESIRED COMP{rescal[i1, nk]:9.6f}")
            w(ishort, " GREATER THAN UNITY.  X NOT RESCALED")

    # Z: the same amount added in every shell; species 5-11 scaled (Fortran lines 122–149).
    if rescal[3, nk] > 0.0:
        if rescal[3, nk] <= 1.0:
            delz = rescal[3, nk] - comp.znew
            for j in range(1, m + 1):
                zfact = max(0.0, hcomp[3, j] + delz) / (hcomp[3, j] + 1.0e-30)
                hcomp[3, j] = max(0.0, hcomp[3, j] + delz)
                # Question (see the page): He3 (species 4) is not scaled.
                for i in range(5, 12):
                    hcomp[i, j] = hcomp[i, j] * zfact
            comp.znew = rescal[3, nk]
            label_mod.zenv0 = rescal[3, nk]
            comp.zenv = label_mod.zenv0
        else:
            i1 = 3
            w(ishort, " ERROR IN SUBROUTINE RSCALE")  # format 1002
            w(ishort, f" RESCALING OF Z IN KIND CARD #{nk:3d} FAILED - "
                      f"DESIRED COMP{rescal[i1, nk]:9.6f}")
            w(ishort, " GREATER THAN UNITY.  Z NOT RESCALED")

    # Main sequence: one species set everywhere (Fortran lines 150–180).
    if hcomp[1, 1] > 1.0e-10:
        # Question (see the page): one xnewcp for all kind cards.
        if nml.lnewcp and nml.xnewcp >= 0.0:
            # A relative target is converted into a local, xtarg.
            xtarg = nml.xnewcp
            if nml.lrel:
                xtarg = math.exp(c.cln * (nml.xnewcp - 12.0)) * hcomp[1, m]
            if xtarg < 1.0:
                for i in range(1, m + 1):
                    hcomp[nml.inewcp, i] = xtarg
            else:
                # Suspected bug (see the page): ICOMP is never set (0 here).
                icomp = 0
                w(ishort, " ERROR IN SUBROUTINE RSCALE")  # format 1004
                w(ishort, f" RESCALING OF SPECIES {icomp:3d} IN KIND CARD # "
                          f"{nk:3d} FAILED - DESIRED COMP{xtarg:9.6f}")
                w(ishort, " GREATER THAN UNITY.  Z NOT RESCALED")

        # The total mass: every mass point scaled (Fortran lines 181–192).
        if rescal[1, nk] > 0.0:
            hstot1 = math.log10(rescal[1, nk] / smass)
            hstot = hstot + hstot1
            smass = rescal[1, nk]
            comp.stotal = hstot
            for i in range(1, m + 1):
                hs[i] = hs[i] + hstot1
    else:
        # Horizontal branch: the old masses (Fortran lines 193–226).
        if rescal[1, nk] > 0.0 or rescal[4, nk] > 0.0:
            jcore, jenv, jxbeg, jxend, jxmid, lshell = findsh(
                hcomp, hl, lc, m, 0, 0, 0, 0, 0, False)
            corold = math.exp(c.cln * hs[jxbeg - 1]) / c.cmsun
            envold = (math.exp(c.cln * hstot)
                      - math.exp(c.cln * hs[jxend])) / c.cmsun
            envold2 = (math.exp(c.cln * hs[m])
                       - math.exp(c.cln * hs[jxend])) / c.cmsun
            shlold = (math.exp(c.cln * hs[jxend])
                      - math.exp(c.cln * hs[jxbeg - 1])) / c.cmsun
            denv = 0.0
            dcor = 0.0
            w(ishort, f" HB-OLD MASSES: CORE {corold:9.6f} SHELL "
                      f"{shlold:9.6f} ENV {envold:9.6f} TOTAL{smass:9.6f}")

        # The total mass: the mass coordinate stretched above m_0 (Fortran lines 228–237).
        if rescal[1, nk] > 0.0:
            denv = rescal[1, nk] - smass
            envnew = envold + denv
            # Question (see the page): the 0.01 Msun switch.
            if ((10.0 ** hstot) / c.cmsun) - corold >= 0.01:
                # The old method: m_0 at the top of the shell (Fortran lines 240–283).
                # Question (see the page): debug output.
                print(" Entering old method ",
                      (10.0 ** hstot - 10.0 ** hs[m])
                      / (10.0 ** hs[m] - 10.0 ** hs[jxbeg - 1]))
                print(" SENV ", 10.0 ** hstot - 10.0 ** hs[m])
                print(" Envelope ", 10.0 ** hs[m] - 10.0 ** hs[jxbeg - 1])
                if envnew <= 0.0:
                    _error_stop([  # format 69
                        " ERROR IN SUBROUTINE RSCALE",
                        " DESIRED NEW ENVELOPE MASS LESS THAN ZERO",
                        f" OLD ENVELOPE MASS {envold:9.2E} NEW ENVELOPE "
                        f"{envold + denv:9.2E} NEW AND OLD TOTAL MASS "
                        f"{rescal[1, nk]:10.2E}{smass:10.2E}"])

                # The scale factor and the new masses (Fortran lines 260–276).
                hse = hs[jxend]
                sfactor = ((rescal[1, nk] - math.exp(c.cln * hse) / c.cmsun)
                           / (smass - math.exp(c.cln * hse) / c.cmsun))
                hstot = math.log10(10.0 ** hse
                                   + sfactor * (10.0 ** hstot - 10.0 ** hse))
                smass = rescal[1, nk]
                comp.stotal = hstot
                for i in range(jxend + 1, m + 1):
                    hs[i] = math.log10(10.0 ** hs[jxend]
                                       + sfactor * (10.0 ** hs[i]
                                                    - 10.0 ** hs[jxend]))
                envold = envnew
            else:
                # The new method: m_0 at the base of the shell (Fortran lines 285–331).
                # Question (see the page): debug output.
                print(" Entering new method ",
                      (10.0 ** hstot - 10.0 ** hs[m])
                      / (10.0 ** hs[m] - 10.0 ** hs[jxbeg - 1]))
                print(" SENV ", 10.0 ** hstot - 10.0 ** hs[m])
                print(" Envelope ", 10.0 ** hs[m] - 10.0 ** hs[jxbeg - 1])

                # Suspected bug (see the page): the new total mass is not checked
                # against the core.
                chkenv = ((10.0 ** hstot) / c.cmsun) - corold
                if chkenv <= 0.0:
                    _error_stop([  # format 69
                        " ERROR IN SUBROUTINE RSCALE",
                        " DESIRED NEW ENVELOPE MASS LESS THAN ZERO",
                        f" OLD ENVELOPE MASS {envold:9.2E} NEW ENVELOPE "
                        f"{envold + denv:9.2E} NEW AND OLD TOTAL MASS "
                        f"{rescal[1, nk]:10.2E}{smass:10.2E}"])

                # The scale factor and the new masses (Fortran lines 307–324).
                hsb = hs[jxbeg - 1]
                sfactor = ((rescal[1, nk] - math.exp(c.cln * hsb) / c.cmsun)
                           / (smass - math.exp(c.cln * hsb) / c.cmsun))
                hstot = math.log10(10.0 ** hsb
                                   + sfactor * (10.0 ** hstot - 10.0 ** hsb))
                smass = rescal[1, nk]
                comp.stotal = hstot
                # HS(JXBEG-1) itself is not in the loop, so hsb stays valid.
                for i in range(jxbeg, m + 1):
                    hs[i] = math.log10(10.0 ** hsb + sfactor
                                       * (10.0 ** hs[i] - 10.0 ** hsb))
                envold = (math.exp(c.cln * hstot)
                          - math.exp(c.cln * hs[jxend])) / c.cmsun
                envtotal = (math.exp(c.cln * hs[m])  # noqa: F841  (unused in the Fortran too)
                            - math.exp(c.cln * hs[jxbeg - 1])) / c.cmsun

        # The core mass: core scaled, shell shifted, envelope squeezed (Fortran lines 337–380).
        if rescal[4, nk] > 0.0:
            shlp = math.exp(c.cln * hs[jxend])
            dcor = rescal[4, nk] - corold
            denv = -dcor
            envnew = envold2 + denv
            if envnew <= 0.0:
                _error_stop([  # format 71
                    " ERROR IN SUBROUTINE RSCALE",
                    " NEW ENVELOPE MASS LESS THAN ZERO BECAUSE OF CORE "
                    "RESCALING",
                    f" OLD ENVELOPE MASS {envold:10.2E} NEW ENVELOPE "
                    f"{envnew:10.2E} NEW AND OLD CORE MASS "
                    f"{rescal[4, nk]:10.2E}{corold:10.2E}"])
            hstot1 = math.log10(rescal[4, nk] / corold)
            for i in range(1, jxbeg):
                hs[i] = hs[i] + hstot1
            shift = dcor * c.cmsun
            for i in range(jxbeg, jxend + 1):
                temp = math.exp(c.cln * hs[i])
                hs[i] = math.log10(temp + shift)
            sfactor = ((math.exp(c.cln * hs[m])
                        - math.exp(c.cln * hs[jxend]))
                       / (math.exp(c.cln * hs[m]) - shlp))
            for i in range(jxend + 1, m + 1):
                hs[i] = math.log10(10.0 ** hs[m]
                                   - sfactor * (10.0 ** hs[m] - 10.0 ** hs[i]))

        # The new masses printed (Fortran lines 382–390).
        if rescal[1, nk] > 0.0 or rescal[4, nk] > 0.0:
            cornew = math.exp(c.cln * hs[jxbeg - 1]) / c.cmsun
            envnew = (math.exp(c.cln * hstot)
                      - math.exp(c.cln * hs[jxend])) / c.cmsun
            shlnew = (math.exp(c.cln * hs[jxend])
                      - math.exp(c.cln * hs[jxbeg - 1])) / c.cmsun
            w(ishort, f" HB-RESCALED MASSES: CORE {cornew:9.6f} SHELL "
                      f"{shlnew:9.6f} ENV {envnew:9.6f} TOTAL {smass:9.6f}")

    # Y = 1 - X - Z - X(He3) everywhere (Fortran lines 392–395).
    for i in range(1, m + 1):
        hcomp[2, i] = 1.0 - hcomp[1, i] - hcomp[3, i] - hcomp[4, i]

    # ZRAMP: Z = Z_c inside RSCLZM1, rising linearly to the surface Z at
    # RSCLZM2 (in fractions of the total mass) (Fortran lines 397–426).
    if (nml.lzramp and nml.rsclzc[nk] > 0.0 and nml.rsclzm1[nk] > 0.0
            and nml.rsclzm2[nk] > 0.0):
        slope = (comp.znew - nml.rsclzc[nk]) / (nml.rsclzm2[nk]
                                                 - nml.rsclzm1[nk])
        for i in range(1, m + 1):
            fs = 10.0 ** hs[i] / (smass * c.cmsun)
            if fs < nml.rsclzm1[nk]:
                zp = nml.rsclzc[nk]
            elif fs < nml.rsclzm2[nk]:
                zp = (fs - nml.rsclzm1[nk]) * slope + nml.rsclzc[nk]
            else:
                continue
            zfact = zp / (hcomp[3, i] + 1.0e-30)
            hcomp[3, i] = zp
            hcomp[1, i] = 1.0 - hcomp[3, i] - hcomp[4, i] - hcomp[2, i]
            for j in range(5, 12):
                hcomp[j, i] = hcomp[j, i] * zfact

    return hstot, smass
