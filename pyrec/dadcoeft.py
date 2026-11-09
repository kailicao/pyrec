"""Python version of dadcoeft.f (YREC 7be4133): angular-momentum transport
by diffusion and advection (meridional circulation) over a step, as a band
system with four unknowns per shell."""

import numpy as np

from pyrec import difad3_mod
from pyrec import difad_mod
from pyrec import dimensions_mod
from pyrec import egridder_mod
from pyrec import nmlpar_mod
from pyrec import ovrtrn_mod
from pyrec.bandw import bandw

JSON = dimensions_mod.json
NMAX = 8000

# SAVE: the local arrays. The first five were COMMON blocks declared only
# here (/EGRIDCHI/, /DIFADDT/, /DIFAD4/) and are set nowhere in the code,
# so they hold zeros (suspected bug, see the page). ETHVN, ETHVP and QCHIRC
# are declared and never used, and are left out.
qchire = np.zeros(JSON + 1)
ewst = np.zeros(JSON + 1)
eqwrst = np.zeros(JSON + 1)
req = np.zeros(JSON + 1)
fgeomix = np.zeros(JSON + 1)
vesn = np.zeros(JSON + 1)
a = np.zeros((NMAX + 1, 10 + 1))
aa = np.zeros((NMAX + 1, 10 + 1))
b = np.zeros(NMAX + 1)
c = np.zeros(NMAX + 1)
d = np.zeros(NMAX + 1)
ewprev = np.zeros(JSON + 1)
dwit = np.zeros(50 + 1)
iwit = np.zeros(50 + 1, dtype=int)
fx1 = np.zeros(JSON + 1)
qwr = np.zeros(JSON + 1)
wm = np.zeros(JSON + 1)
ewit = np.zeros(JSON + 1)
wm0 = np.zeros(JSON + 1)
ewit2 = np.zeros(JSON + 1)
wmp = np.zeros(JSON + 1)
qwrp = np.zeros(JSON + 1)
ewt = np.zeros(JSON + 1)
wminit = np.zeros(JSON + 1)
qwrinit = np.zeros(JSON + 1)
qqwcc = np.zeros(JSON + 1)


def dadcoeft(dr: float, dt: float, ei: np.ndarray, ew: np.ndarray, ntot: int, wind1: float, wind2: float,
             dj: np.ndarray, ecod2: np.ndarray, sumdj: float, lfix: bool, lokad: bool) -> tuple[float, float, float,
             bool]:
    """Transport angular momentum on the even grid over a step with
    diffusion and advection: four unknowns per shell (Omega, Omega'',
    Omega', Omega''') in a band system solved by BANDW, iterated on the
    Omega-dependent coefficients, rescaled to conserve angular momentum
    with the wind, in substeps that are halved on failure; then the
    composition mixing coefficient.

    Translates: dadcoeft.f#DADCOEFT (7be4133)

    Arguments: dr, the grid spacing; dt, the step; ei, the moments of
    inertia; ew, Omega at the start; ntot, the number of points; wind1,
    wind2, the wind's change of Omega (start, end); dj, ecod2 (arrays,
    set in place); sumdj (set only); lfix, a fixed surface Omega; lokad
    (set only).

    Returns (dt, wind2, sumdj, lokad): dt as passed (changed during the
    substeps, restored at the end); wind2 as last rescaled; the total
    change of angular momentum; and convergence (True on return, since
    failure stops the run).

    """

    d3 = difad3_mod
    ecod3 = difad_mod.ecod3
    ecod4 = difad_mod.ecod4
    eqqcod = egridder_mod.eqqcod
    eqqqcod = egridder_mod.eqqqcod
    nml = nmlpar_mod
    ess = d3.ess
    egsf = d3.egsf

    # Initial values (Fortran lines 144–148).
    ncut = 0
    nstep = 1
    dt0 = dt
    wind20 = wind2
    tiny = 1.0e-30

    # Start-of-step gradients and averaged Omega (Fortran lines 149–154).
    for i in range(2, ntot + 1):
        qwr[i] = qchire[i] * (ew[i] - ew[i - 1]) / dr
        wm0[i] = 0.5 * (ew[i] + ew[i - 1])

    # Suspected bug (see the page): nm must not exceed NMAX; nothing checks
    # it (Fortran line 155).
    nm = 4 * ntot - 2

    while True:  # label 5: restart with more substeps
        # Start-of-step values of Omega and J (Fortran lines 156–172).
        sumdt = 0.0
        wind2 = wind20
        sumj = 0.0
        for i in range(1, ntot + 1):
            ewprev[i] = ew[i]
            ewit[i] = ew[i]
            ewit2[i] = ew[i]
            ewt[i] = ew[i]
            sumj = sumj + ew[i] * ei[i]
        for i in range(2, ntot + 1):
            wmp[i] = wm0[i]
            qwrp[i] = qwr[i]

        # The substeps; goto950 stands for GOTO 950 (Fortran lines
        # 173–176).
        goto950 = False
        for nnnn in range(1, nstep + 1):  # substeps
            dt = dt0 / float(nstep)
            sumdt = sumdt + dt
            fxx = dt / dt0

            # Outer iteration, meant for the d(theta)/dt term (now zero).
            # Question (see the page): can this loop go? (Fortran lines
            # 177–183.)
            for nnn in range(1, nml.itdif2 + 1):
                lokad = False

                # Rescale the wind to the current surface Omega (Fortran
                # lines 184–211).
                if nnn > 1 and wind2 > 0.0:
                    if nml.lrossby:
                        wsat = nml.wmax * nml.pmmsoltau / ovrtrn_mod.taucz
                    else:
                        wsat = nml.wmax
                    wp = min(ew[ntot], wsat)
                    wn = min(ewprev[ntot], wsat)
                    wind2 = (wind20 * (wn / wp)**(nml.exw - 1.0)
                             * (ewprev[ntot] / ew[ntot]))

                # Inner iteration on the coefficients. In Fortran, nn ends
                # at itdif2 + 1 when the loop runs out (used in WRITE 100)
                # (Fortran line 212).
                nn = 0
                while nn < nml.itdif2:
                    nn = nn + 1

                    # Start values for the coefficients; unused? (read only
                    # in commented-out lines) (Fortran lines 213–225).
                    if nnnn == 1:
                        for i in range(2, ntot + 1):
                            wminit[i] = 0.5 * (ewst[i] + ewst[i - 1])
                            qwrinit[i] = eqwrst[i]
                    else:
                        for i in range(2, ntot + 1):
                            wminit[i] = wmp[i]
                            qwrinit[i] = qwrp[i]

                    # Second differences of Omega (Fortran lines 226–231).
                    qqwcc[1] = 0.0
                    for i in range(2, ntot):
                        qqwcc[i] = ewprev[i + 1] - 2.0 * ewprev[i] + ewprev[i - 1]
                    qqwcc[ntot] = 0.0

                    # Coefficients from the current Omega. Suspected bug
                    # (see the page): ev0..ev2b, fgeom, ess, egsf and qchire
                    # are set by no routine (Fortran lines 232–271).
                    for i in range(2, ntot + 1):
                        wm[i] = (0.5 * (ewprev[i] + ewprev[i - 1])
                                 - 0.125 * (qqwcc[i] - qqwcc[i - 1]))
                        wmi = 0.5 * (ewit2[i] + ewit2[i - 1])  # noqa: F841  (unused in the Fortran too)
                        qwr2 = qchire[i] * (ewprev[i] - ewprev[i - 1]) / dr
                        qwri = qchire[i] * (ewit2[i] - ewit2[i - 1]) / dr  # noqa: F841  (unused in the Fortran too)
                        if abs(qwr[i]) > tiny:
                            ess[i] = ess[i] * (qwr2 / qwr[i])**2
                            egsf[i] = egsf[i] * (qwr2 / qwr[i])**2
                        v1 = nml.fw * d3.ev0[i] * wm[i]**2 * (d3.ev1a[i] + wm[i]**2 * d3.ev1b[i])
                        vthn = 0.0
                        vthp = 0.0
                        if nnn == 1:
                            vthn = 0.0
                            vthp = 0.0
                        v2 = nml.fw * d3.ev0[i] * wm[i]**2 * (d3.ev2a[i] + d3.ev2b[i]) + vthn
                        v1 = v1 + vthp
                        # Question (see the page): the factors 0.2 here and
                        # 5 in vesn below.
                        ecod3[i] = 0.2 * v1 * d3.fgeom[i]
                        ecod4[i] = (0.2 * v2 + ess[i] + egsf[i]) * d3.fgeom[i]
                        ecod4[i] = max(0.0, ecod4[i])

                    # Correction factors for Omega (Fortran lines 272–277).
                    for i in range(2, ntot + 1):
                        fx1[i] = 2.0 * wm[i] / (ew[i] + ew[i - 1])

                    # Right-hand side: Omega, then three zero rows per shell
                    # (Fortran lines 278–289).
                    for i in range(1, ntot + 1):
                        ii = 1 + 4 * (i - 1)
                        b[ii] = ewt[i]
                        b[ii + 1] = 0.0
                        b[ii + 2] = 0.0
                        b[ii + 3] = 0.0
                        d[ii] = ewt[i]
                        d[ii + 1] = 0.0
                        d[ii + 2] = 0.0
                        d[ii + 3] = 0.0

                    # The wind's loss of angular momentum in the surface row
                    # (Fortran lines 290–293).
                    delj = -0.5 * (wind1 + wind2) * ei[ntot] * fxx
                    b[4 * ntot - 3] = b[4 * ntot - 3] * (1.0 + delj / ei[ntot] / ewt[ntot])
                    d[4 * ntot - 3] = b[4 * ntot - 3]

                    # The band matrix, column 5 the diagonal: factors, zeros
                    # (Fortran lines 294–311).
                    fact0 = dt / dr
                    fact = fact0 / ei[1]
                    facta = 0.5 * dt / ei[1] / dr
                    for i in range(1, nm + 1):
                        for j in range(1, 11):
                            a[i, j] = 0.0

                    # The first shell's rows (Fortran lines 312–336).
                    a[1, 5] = 1.0 - facta * ecod3[2]
                    a[1, 9] = -facta * ecod3[2]
                    a[1, 7] = -fact * qchire[2] * ecod4[2]
                    a[1, 10] = -facta * qchire[2] * d3.ecod5[2] * fx1[2]
                    a[1, 8] = fact * qchire[2] * d3.ecod6[2] * fx1[2]
                    a[2, 5] = 1.0
                    a[3, 3] = 1.0 / dr
                    a[3, 5] = 1.0
                    a[3, 7] = -1.0 / dr
                    a[4, 5] = 1.0
                    a[4, 9] = -1.0 / 3.0

                    # The interior shells' rows (Fortran lines 337–374).
                    for ii in range(2, ntot):
                        fact = fact0 / ei[ii]
                        facta = 0.5 * dt / ei[ii] / dr
                        i = 1 + (ii - 1) * 4
                        a[i, 1] = facta * ecod3[ii]
                        a[i, 5] = 1.0 + facta * (ecod3[ii] - ecod3[ii + 1])
                        a[i, 9] = -facta * ecod3[ii + 1]
                        a[i, 3] = fact * qchire[ii] * ecod4[ii]
                        a[i, 7] = -fact * qchire[ii + 1] * ecod4[ii + 1]
                        a[i, 2] = facta * qchire[ii] * d3.ecod5[ii] * fx1[ii]
                        a[i, 6] = facta * (d3.ecod5[ii] * qchire[ii] * fx1[ii]
                                           - d3.ecod5[ii + 1] * qchire[ii + 1] * fx1[ii + 1])
                        a[i, 10] = -facta * d3.ecod5[ii + 1] * qchire[ii + 1] * fx1[ii + 1]
                        a[i, 4] = fact * d3.ecod6[ii] * qchire[ii] * fx1[ii]
                        a[i, 8] = -fact * d3.ecod6[ii + 1] * qchire[ii + 1] * fx1[ii + 1]
                        a[i + 1, 2] = eqqcod[ii] * qchire[ii] / dr * fx1[ii]
                        a[i + 1, 5] = 1.0
                        a[i + 1, 6] = -eqqcod[ii + 1] * qchire[ii + 1] / dr * fx1[ii + 1]
                        a[i + 2, 3] = 1.0 / dr
                        a[i + 2, 5] = 1.0
                        a[i + 2, 7] = -1.0 / dr
                        a[i + 3, 3] = qchire[ii] * eqqqcod[ii] / dr
                        a[i + 3, 5] = 1.0
                        a[i + 3, 7] = -qchire[ii] * eqqqcod[ii + 1] / dr

                    # The last shell's rows (Fortran lines 375–406).
                    fact = fact0 / ei[ntot]
                    facta = 0.5 * dt / ei[ntot] / dr
                    i = 4 * ntot - 3
                    a[i - 4, 10] = 0.0
                    a[5, 2] = 0.0
                    if not lfix:
                        a[i, 1] = facta * ecod3[ntot]
                        a[i, 5] = 1.0 + facta * ecod3[ntot]
                        a[i, 3] = fact * qchire[ntot] * ecod4[ntot]
                        a[i, 2] = facta * qchire[ntot] * d3.ecod5[ntot] * fx1[ntot]
                        a[i, 4] = fact * qchire[ntot] * d3.ecod6[ntot]
                    else:
                        a[i, 5] = 1.0
                    a[i + 1, 5] = 1.0

                    # Keep the matrix, and solve (Fortran lines 407–427).
                    for i in range(1, 11):
                        for j in range(1, nm + 1):
                            aa[j, i] = a[j, i]
                    bandw(a, nm, b)

                    # Residual check, computed and never used; unused?
                    # (Fortran lines 428–463).
                    for i in range(1, nm + 1):
                        c[i] = 0.0
                    for i in range(5, 11):
                        c[1] = c[1] + aa[1, i] * b[i - 4]
                    for i in range(4, 11):
                        c[2] = c[2] + aa[2, i] * b[i - 3]
                    for i in range(3, 11):
                        c[3] = c[3] + aa[3, i] * b[i - 2]
                    for i in range(2, 11):
                        c[4] = c[4] + aa[4, i] * b[i - 1]
                    for j in range(5, nm - 4):
                        for i in range(1, 11):
                            c[j] = c[j] + aa[j, i] * b[i - 5 + j]
                    for i in range(1, 10):
                        c[nm - 4] = c[nm - 4] + aa[nm - 4, i] * b[i + nm - 9]
                    for i in range(1, 9):
                        c[nm - 3] = c[nm - 3] + aa[nm - 3, i] * b[i + nm - 8]
                    for i in range(1, 8):
                        c[nm - 2] = c[nm - 2] + aa[nm - 2, i] * b[i + nm - 7]
                    for i in range(1, 7):
                        c[nm - 1] = c[nm - 1] + aa[nm - 1, i] * b[i + nm - 6]
                    for i in range(1, 6):
                        c[nm] = c[nm] + aa[nm, i] * b[i + nm - 5]

                    # Largest relative changes: from the last iterate, from
                    # the start of the step, from the last outer iterate
                    # (Fortran lines 464–498).
                    sumdj = 0.0
                    sumdj2 = 0.0
                    dwmax = (b[1] - ewprev[1]) / ewprev[1]
                    iwmax = 1
                    dww = (b[1] - ew[1]) / ew[1]
                    iww = 1
                    dwi = (b[1] - ewit[1]) / ewit[1]
                    iwi = 1  # noqa: F841  (unused in the Fortran too)
                    for i in range(1, ntot + 1):
                        ii = 1 + 4 * (i - 1)
                        dj[i] = (b[ii] - ew[i]) * ei[i]
                        sumdj2 = sumdj2 + (b[ii] - ewt[i]) * ei[i]
                        sumdj = sumdj + dj[i]
                        dw = (b[ii] - ewprev[i]) / ewprev[i]
                        dw2 = (b[ii] - ew[i]) / ew[i]
                        dwi2 = (b[ii] - ewit[i]) / ewit[i]
                        if abs(dw) > abs(dwmax):
                            iwmax = i
                            dwmax = dw
                        if abs(dw2) > abs(dww):
                            iww = i
                            dww = dw2
                        if abs(dwi2) > abs(dwi):
                            iwi = i  # noqa: F841  (unused in the Fortran too)
                            dwi = dwi2

                    # Corrections over 10%: cut the step; the damping FX is
                    # always 1 (Fortran lines 499–512). Question (see the
                    # page): the screen output.
                    if abs(dwmax) > 0.1:
                        print(nnnn, nnn, nn)
                        print(" CORRECTIONS TOO LARGE ", dwmax, iwmax)
                        lokad = False
                        goto950 = True
                        break
                    elif abs(dwmax) < 1.0e-2:
                        fx = 1.0
                    else:
                        fx = 1.0

                    # Update, then rescale to conserve angular momentum
                    # (Fortran lines 513–529).
                    sumj2 = 0.0
                    for i in range(1, ntot + 1):
                        ii = 1 + 4 * (i - 1)
                        ewprev[i] = ewprev[i] + fx * (b[ii] - ewprev[i])
                        dj[i] = dj[i] * fx
                        sumj2 = sumj2 + ewprev[i] * ei[i]
                    corr = (sumj + delj) / sumj2
                    print(sumj, delj, sumj2, corr)
                    sumdj = 0.0
                    for i in range(1, ntot + 1):
                        ewprev[i] = ewprev[i] * corr
                        dj[i] = (ewprev[i] - ew[i]) * ei[i]
                        sumdj = sumdj + dj[i]
                    print(sumdj, delj)

                    # History, and inner convergence (Fortran lines
                    # 534–540).
                    dwit[nn] = dwmax
                    iwit[nn] = iwmax
                    if abs(dwmax) < nml.djok:
                        lokad = True
                        break  # GOTO 900
                else:
                    nn = nn + 1  # Fortran's DO leaves nn = itdif2 + 1
                if goto950:
                    break

                # label 900. Suspected bug (see the page): dwit[itdif2 + 1]
                # is printed unset (Fortran lines 542–546).
                print(nnnn, nnn, nn)
                print(f" MAX D(OMEGA)/OMEGA{dww:12.3E} AT PT.{iww:5d} BY ITERATION")
                items = [f" {dwit[j]:11.3E}{iwit[j]:4d}" for j in range(1, nn + 1)]
                for k in range(0, len(items), 5):
                    print("".join(items[k:k + 5]))

                # The outer iterate (Fortran lines 547–557).
                if nnn <= 2:
                    for k in range(1, ntot + 1):
                        ewit[k] = ewprev[k]
                        ewit2[k] = ewprev[k]
                else:
                    for k in range(1, ntot + 1):
                        ewit[k] = ewprev[k]
                        ewit2[k] = 0.5 * ewit[k] + 0.5 * ewprev[k]

                # Outer convergence. Suspected bugs (see the page): never
                # with itdif2 = 1; a failed substep is not cut (Fortran
                # lines 558–562).
                if abs(dwi) <= nml.djok and nnn >= 2:
                    break  # GOTO 9999
                else:
                    lokad = False
            if goto950:
                break

            # label 9999: end of the substep (Fortran lines 564–572).
            for k in range(1, ntot + 1):
                ewt[k] = ewprev[k]
            sumj = sumj + delj
            for k in range(2, ntot + 1):
                wmp[k] = 0.5 * (ewprev[k] + ewprev[k - 1])
                qwrp[k] = qchire[k] * (ewprev[k] - ewprev[k - 1]) / dr

        # label 950: cut the step, at most 5 times (Fortran lines 574–585).
        if not lokad:
            ncut = ncut + 1
            nstep = 2 * nstep
            if ncut > 5:
                print(" TIMESTEP CUT MORE THAN 5 TIMES - RUN STOPPED")
                raise SystemExit("111")
            else:
                print(" TIMESTEP CUT #", ncut, " IN DADCOEFT")
            continue  # GOTO 5
        break

    # Fortran lines 586–587.
    print(f"{b[1]:12.3E}{ew[1]:12.3E}{b[2]:12.3E}{b[3]:12.3E}{b[4]:12.3E}")

    # Composition mixing from the final solution (Fortran lines 588–610).
    for i in range(2, ntot + 1):
        ii = 1 + 4 * (i - 1)
        wmm = 0.5 * (b[ii] + b[ii - 4])
        w2m = 0.5 * (b[ii + 1] + b[ii - 3])
        v0 = ecod3[i] / d3.fgeom[i]
        v1 = ecod4[i] / d3.fgeom[i] * b[ii - 2] * qchire[i] / wmm
        v1e = 0.2 * nml.fw * d3.ev0[i] * wmm * qchire[i] * (d3.ev2a[i] + d3.ev2b[i]) * b[ii - 2]
        if ecod4[i] <= 0.0:
            v1e = 0.0
        v2 = d3.ecod5[i] / d3.fgeom[i] * qchire[i] * w2m / wmm
        v3 = d3.ecod6[i] / d3.fgeom[i] * b[ii - 1] * qchire[i] / wmm
        vtot = v0 + v1 + v2 + v3
        vtot2 = v0 + v1e + v2 + v3
        if i <= 10 or (ntot - i) <= 10:
            print("".join(f"{v:12.3E}" for v in (
                b[ii], ew[i], v0, v1, v2, v3, vtot, b[ii + 1], b[ii + 2], b[ii + 3])))
        vesn[i] = 5.0 * abs(nml.fc * nml.fesc * vtot2 * req[i])
        dcmix = vesn[i] + nml.fc * nml.fssc * ess[i] + nml.fc * nml.fgsfc * egsf[i]
        ecod2[i] = dcmix * fgeomix[i]

    # Fortran lines 611–613.
    dt = dt0
    delj0 = delj / fxx
    print(sumdj, delj0)
    return dt, wind2, sumdj, lokad
