"""Python version of condopacp.f (YREC 09c3b9b): the electron thermal
conductivity interpolated in A. Y. Potekhin's table (condall.d)."""

import math

import numpy as np

from pyrec import iounits_mod

# The table's dimensions; they must match the file.
MAXT = 19
MAXR = 64
MAXZ = 15

# Locals kept between calls (SAVE, DATA): the table, read at the first
# call, and the last indices found, HUNT's first guesses at the next call.
at = np.zeros(MAXT + 1)
ar = np.zeros(MAXR + 1)
az = np.zeros(MAXZ + 1)
akap = np.zeros((MAXT + 1, MAXR + 1, MAXZ + 1))
krun = -1
iz = 0
it = 0
ir = 0


def _read_values(unit, n):
    """A list-directed READ of n numbers: records are read until n values
    are found."""
    vals = []
    while len(vals) < n:
        rec = iounits_mod.read(unit)
        vals.extend(float(v.replace('D', 'E').replace('d', 'e'))
                    for v in rec.replace(',', ' ').split())
    return vals[:n]


def condopacp(zion, tlg, rlg, ck, drk, dtk):
    """Interpolate log10 of the electron thermal conductivity, and its
    derivatives in log rho and log T, in Potekhin's table: cubic in
    log rho and log T, linear in log Z_ion. Reads the table at the first
    call.

    Translates: condopacp.f#CONDOPACP (09c3b9b)

    Arguments: zion, the ion charge; tlg, log T (K); rlg, log rho
    (g/cm^3); ck, drk, dtk, passed in (with any value).

    Returns (ck, drk, dtk): log10 K (cgs), d log K / d log rho and
    d log K / d log T.
    """
    global krun, iz, it, ir
    if krun != 12345:
        # Read the table (the file is opened by PARMIN).
        ip = iounits_mod.icondopacp
        iounits_mod.read(ip)  # skip the first line
        for iz in range(1, MAXZ + 1):
            vals = _read_values(ip, 1 + MAXT)
            z = vals[0]
            for it in range(1, MAXT + 1):
                at[it] = vals[it]
            az[iz] = math.log10(z)
            for ir in range(1, MAXR + 1):
                vals = _read_values(ip, 1 + MAXT)
                ar[ir] = vals[0]
                for it in range(1, MAXT + 1):
                    akap[it, ir, iz] = vals[it]
        iounits_mod.close(ip)
        krun = 12345
        # The middle of the table: INT(MAXZ/2.)+1, the same as MAXZ/2+1 (KC 2025-05-30).
        iz = int(MAXZ / 2.0) + 1
        it = int(MAXT / 2.0) + 1
        ir = int(MAXR / 2.0) + 1
    zlg = math.log10(zion)
    iz = hunt(az, MAXZ, zlg, iz)
    if iz == 0 or iz == MAXZ:
        raise SystemExit('CONINTER: Z out of range')
    it = hunt(at, MAXT, tlg, it)
    # Question (see the page): a point exactly at the last grid value
    # (log T = 9, log rho = 9.75) also stops here.
    if it == 0 or it == MAXT:
        # A debugging leftover: the Z grid is printed before stopping.
        print(*az[1:])
        print(MAXZ, tlg, iz)
        raise SystemExit('CONINTER: T out of range')
    ir = hunt(ar, MAXR, rlg, ir)
    if ir == 0 or ir == MAXR:
        raise SystemExit('CONINTER: rho out of range')
    itm = max(1, it - 1)
    itp = min(MAXT, it + 2)
    irm = max(1, ir - 1)
    irp = min(MAXR, ir + 2)

    # Cubic interpolation in log rho, at the 4 temperatures ITM, IT,
    # IT+1, ITP and the 2 charges IZ, IZ+1.
    def in_rho(jt, jz):
        vf, df, d2, xr = cinterp3(ar[irm], ar[ir], ar[ir + 1], ar[irp],
                                  rlg, ir, MAXR, akap[jt, irm, jz],
                                  akap[jt, ir, jz], akap[jt, ir + 1, jz],
                                  akap[jt, irp, jz], 0.0, 0.0, 0.0, 0.0)
        return vf, df, d2

    cktmz0, drktmz0, dr2ktmz0 = in_rho(itm, iz)
    ckt0z0, drkt0z0, dr2kt0z0 = in_rho(it, iz)
    ckt1z0, drkt1z0, dr2kt1z0 = in_rho(it + 1, iz)
    cktpz0, drktpz0, dr2ktpz0 = in_rho(itp, iz)
    cktmz1, drktmz1, dr2ktmz1 = in_rho(itm, iz + 1)
    ckt0z1, drkt0z1, dr2kt0z1 = in_rho(it, iz + 1)
    ckt1z1, drkt1z1, dr2kt1z1 = in_rho(it + 1, iz + 1)
    cktpz1, drktpz1, dr2ktpz1 = in_rho(itp, iz + 1)

    # Linear interpolation in log Z.
    xz1 = (zlg - az[iz]) / (az[iz + 1] - az[iz])
    xz0 = 1.0 - xz1
    cktm = xz0 * cktmz0 + xz1 * cktmz1
    drktm = xz0 * drktmz0 + xz1 * drktmz1
    dr2ktm = xz0 * dr2ktmz0 + xz1 * dr2ktmz1
    ckt0 = xz0 * ckt0z0 + xz1 * ckt0z1
    drkt0 = xz0 * drkt0z0 + xz1 * drkt0z1
    dr2kt0 = xz0 * dr2kt0z0 + xz1 * dr2kt0z1
    ckt1 = xz0 * ckt1z0 + xz1 * ckt1z1
    drkt1 = xz0 * drkt1z0 + xz1 * drkt1z1
    dr2kt1 = xz0 * dr2kt1z0 + xz1 * dr2kt1z1
    cktp = xz0 * cktpz0 + xz1 * cktpz1
    drktp = xz0 * drktpz0 + xz1 * drktpz1
    dr2ktp = xz0 * dr2ktpz0 + xz1 * dr2ktpz1

    # Cubic interpolation in log T: of log K (giving d log K / d log T),
    # and of d log K / d log rho.
    ck, dtk, dt2k, xt = cinterp3(at[itm], at[it], at[it + 1], at[itp], tlg,
                                 it, MAXT, cktm, ckt0, ckt1, cktp,
                                 0.0, 0.0, 0.0, 0.0)
    drk, drtk, drt2k, xt = cinterp3(at[itm], at[it], at[it + 1], at[itp],
                                    tlg, it, MAXT, drktm, drkt0, drkt1,
                                    drktp, 0.0, 0.0, 0.0, 0.0)
    return ck, drk, dtk


def cinterp3(zm, z0, z1, zp, z, n0, mxnv, vm, v0, v1, vp, vf, df, d2, xh):
    """Interpolate at z on [z0, z1] from 4 tabulated points: a cubic
    Hermite polynomial with end slopes from the neighboring secants, and
    quadratic on the table's first and last intervals.

    Translates: condopacp.f#CINTERP3 (09c3b9b)

    Arguments: zm, z0, z1, zp and vm, v0, v1, vp, the 4 points; z; n0, the
    index of z0 in the table; mxnv, the table's length; vf, df, d2, xh,
    passed in (with any value).

    Returns (vf, df, d2, xh): the value, the first and the second
    derivative, and the fraction of the way from z0 to z1.
    """
    if n0 <= 0 or n0 >= mxnv:
        raise SystemExit('CINTERP: N0 out of range')
    x = z - z0
    h = z1 - z0  # the basic interval
    xh = x / h
    if n0 > 1:
        hm = z0 - zm  # the left adjoining interval
        v01 = (((v1 - v0) / h ** 2 + (v0 - vm) / hm ** 2)
               / (1.0 / h + 1.0 / hm))  # the left derivative
    if n0 < mxnv - 1:
        hp = zp - z1  # the right adjoining interval
        v11 = (((v1 - v0) / h ** 2 + (vp - v1) / hp ** 2)
               / (1.0 / h + 1.0 / hp))  # the right derivative
    if n0 > 1 and n0 < mxnv - 1:
        # Cubic interpolation.
        c2 = 3.0 * (v1 - v0) - h * (v11 + 2.0 * v01)
        c3 = h * (v01 + v11) - 2.0 * (v1 - v0)
        vf = v0 + v01 * x + c2 * xh ** 2 + c3 * xh ** 3
        df = v01 + (2.0 * c2 * xh + 3.0 * c3 * xh ** 2) / h
        d2 = (2.0 * c2 + 6.0 * c3 * xh) / h ** 2
        return vf, df, d2, xh  # GOTO 10
    if n0 == 1:
        # Quadratic interpolation, first interval.
        c2 = v0 - v1 + v11 * h
        vf = v1 - v11 * (h - x) + c2 * (1.0 - xh) ** 2
        df = v11 - 2.0 * c2 * (1.0 - xh) / h
        d2 = 2.0 * c2 / h ** 2
    else:
        # Quadratic interpolation, last interval (N0 = MXNV-1).
        c2 = v1 - v0 - v01 * h
        vf = v0 + v01 * x + c2 * xh ** 2
        df = v01 + 2.0 * c2 * xh / h
        d2 = 2.0 * c2 / h ** 2
    return vf, df, d2, xh


def hunt(xx, n, x, jlo):
    """Numerical Recipes' HUNT: find jlo with x between xx[jlo] and
    xx[jlo+1] in a monotonic table, starting from the guess jlo and
    doubling the step before bisecting.

    Translates: condopacp.f#HUNT (09c3b9b)

    Arguments: xx[1..n], the table; n; x; jlo, the first guess.

    Returns jlo: 0 or n outside the table (n also for x equal to xx[n]).
    """
    ascnd = xx[n] > xx[1]  # True if ascending
    if jlo <= 0 or jlo > n:
        # The guess is not useful: go to bisection.
        jlo = 0
        jhi = n + 1
    else:
        inc = 1
        if (x >= xx[jlo]) == ascnd:
            # Hunt up (label 1).
            while True:
                jhi = jlo + inc
                if jhi > n:
                    jhi = n + 1  # off the end of the table
                    break
                elif (x >= xx[jhi]) == ascnd:
                    jlo = jhi
                    inc = inc + inc
                    continue
                break
        else:
            # Hunt down (label 2).
            jhi = jlo
            while True:
                jlo = jhi - inc
                if jlo < 1:
                    jlo = 0  # off the end of the table
                    break
                elif (x < xx[jlo]) == ascnd:
                    jhi = jlo
                    inc = inc + inc
                    continue
                break
    # label 3: bisection.
    while jhi - jlo != 1:
        jm = (jhi + jlo) // 2
        if (x >= xx[jm]) == ascnd:
            jlo = jm
        else:
            jhi = jm
    return jlo
