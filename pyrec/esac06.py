"""Python version of esac06.f (YREC development, 09c3b9b): interpolation
of the OPAL 2006 equation-of-state tables in X, T6 and rho, and their
one-time setup."""

import numpy as np

from pyrec import aeos06_mod
from pyrec import eeeos06_mod
from pyrec import eos06misc_mod
from pyrec import iounits_mod
from pyrec.gmass06 import gmass06
from pyrec.quadeos06 import quadeos06
from pyrec.radsub06 import radsub06
from pyrec.readcoeos06 import readcoeos06
from pyrec.t6rinteos06 import t6rinteos06

MX = 5                    # X tables
MV = 10                   # quantities
NR = 169                  # densities
NT = 197                  # temperatures
APROP = 83.14511          # gas constant, 1e12 erg / mole / 1e6 K (DATA)
ID = "OPALEOS/ESAC06:"

# Locals kept between calls (SAVE).
frac = np.zeros(8)
moles = 0.0


def esac06(xh, t6, r, iorder, irad):
    """Interpolate the OPAL 2006 EOS tables in X (overlapping quadratics),
    then T6 and rho (t6rinteos06), and scale the results; the ten
    quantities go to eos06misc_mod.eos.

    Translates: esac06.f#ESAC06 (09c3b9b)

    Arguments: xh (X), t6 (T in 1e6 K), r (density), iorder (how many
    quantities), irad (1: add radiation through radsub06). None is set.

    Returns alt1, True when the Fortran takes its alternate return (*):
    out of the table (label 62) or in its empty region (label 65).
    """
    global moles
    a = aeos06_mod
    e = eeeos06_mod
    m = eos06misc_mod
    io = iounits_mod

    if iorder > 9:
        # Question (see the page): the table has ten quantities.
        io.write(io.ishort, f"{ID} iorder cannot exceed 9")
    if irad != 0 and irad != 1:
        io.write(io.ishort, f"{ID} Irad must be 0 or 1")
        raise SystemExit(f"{ID} Irad must be 0 or 1")
    xxi = xh                                       # unused, as in the Fortran
    ri = r                                         # unused, as in the Fortran
    slt = t6
    slr = r
    # Z from the table, every call (the tables are read by seteos06).
    z = m.zz[1]

    # Is the point within the grid? (t6a decreases with its index.)
    if slt > a.t6a[1] or slt < a.t6a[NT] or slr < a.rho[1] or slr > a.rho[NR]:
        return _out_of_table(slt, slr)             # go to 62

    # Bisection for the X tables mf, mg, mh, mi.
    ilo = 2
    ihi = MX
    while ihi - ilo > 1:                           # label 8
        imd = (ihi + ilo) // 2
        if xh <= a.xa[imd] + 1.0e-7:
            ihi = imd
        else:
            ilo = imd
    i = ihi
    a.mf = i - 2
    mg = i - 1
    mh = i
    mi = i + 1
    mf2 = mi
    if xh < 1.0e-6:
        a.mf = 1
        mg = 1
        mh = 1
        mi = 2
        mf2 = 1
    if xh <= a.xa[2] + 1.0e-7 or xh >= a.xa[MX - 2] - 1.0e-7:
        mf2 = mh
    mf = a.mf

    # Bisection in rho: l1 ... l4.
    ilo = 2
    ihi = NR
    while ihi - ilo > 1:                           # label 12
        imd = (ihi + ilo) // 2
        if slr == a.rho[imd]:
            ihi = imd
            break                                  # go to 13
        if slr <= a.rho[imd]:
            ihi = imd
        else:
            ilo = imd
    i = ihi                                        # label 13
    m.l1 = i - 2
    m.l2 = i - 1
    m.l3 = i
    m.l4 = m.l3 + 1
    iqu = 3
    if m.l4 > NR:
        iqu = 2

    # Bisection in T6 (along t6list(1, .)): k1 ... k4.
    ilo = NT
    ihi = 2
    while ilo - ihi > 1:                           # label 11
        imd = (ihi + ilo) // 2
        if t6 == a.t6list[1, imd]:
            ilo = imd
            break                                  # go to 14
        if t6 <= a.t6list[1, imd]:
            ihi = imd
        else:
            ilo = imd
    i = ilo                                        # label 14
    m.k1 = i - 2
    m.k2 = i - 1
    m.k3 = i
    m.k4 = m.k3 + 1
    ipu = 3
    if m.k4 > NT:
        ipu = 2
    if m.k3 == 0:
        io.write(io.ishort, f"{ID} ihi,ilo,imd{ihi:5d}{ilo:5d}{imd:5d}")

    # Does the stencil avoid the empty region (values above 1e30)?
    # Choose the largest stencil allowed.
    sum1 = 0.0
    sum2 = 0.0
    sum23 = 0.0
    sum33 = 0.0
    l1, k1 = m.l1, m.k1
    for midx in range(mf, mf2 + 1):
        a.midx = midx
        for ir in range(l1, l1 + 2):
            for it in range(k1, k1 + 2):
                sum1 = sum1 + a.xz[midx, 1, it, ir]
        for ir in range(l1, l1 + 3):
            for it in range(k1, k1 + 3):
                sum2 = sum2 + a.xz[midx, 1, it, ir]
        if ipu == 3:
            for ir in range(l1, l1 + 3):
                for it in range(k1, k1 + ipu + 1):
                    sum23 = sum23 + a.xz[midx, 1, it, ir]
        else:
            sum23 = 2.0e30
        if iqu == 3:
            for ir in range(l1, l1 + 4):
                for it in range(k1, k1 + ipu + 1):
                    sum33 = sum33 + a.xz[midx, 1, it, ir]
        else:
            sum33 = 2.0e30
    m.iq = 2
    m.ip = 2
    shifted = False
    if sum2 > 1.0e30:
        if sum1 < 1.0e25:
            # Move the stencil back by one in T6 and rho, used as 3x3
            # without a new check.
            m.k1 = m.k3 - 3
            m.k2 = m.k1 + 1
            m.k3 = m.k2 + 1
            m.l1 = m.l3 - 3
            m.l2 = m.l1 + 1
            m.l3 = m.l2 + 1
            shifted = True                         # go to 15
        else:
            return _empty_region(xh, t6, r)        # go to 65
    if not shifted:
        if sum23 < 1.0e30:
            m.ip = 3
        if sum33 < 1.0e30:
            m.iq = 3
        if t6 >= a.t6list[1, 2] + 1.0e-7:
            m.ip = 2
        if slr <= a.rho[2] + 1.0e-15:
            m.iq = 2
        if m.l3 == NR or m.k3 == NT:
            m.iq = 2
            m.ip = 2

    # label 15: for each quantity, interpolate in X, then in T6 and rho.
    l1, k1, ip, iq = m.l1, m.k1, m.ip, m.iq
    for iv in range(1, iorder + 1):
        for midx in range(mf, mf2 + 1):
            a.midx = midx
            is_ = 0
            for ir in range(l1, l1 + iq + 1):
                for it in range(k1, k1 + ip + 1):
                    e.epl[midx, it, ir] = a.xz[midx, iv, it, ir]
                    is_ = 1
        if m.zz[mg] != m.zz[mf] or m.zz[mh] != m.zz[mf]:
            io.write(io.ishort,
                     f"{ID}Z does not match Z in OEOS06 files you are using")
            raise SystemExit(
                f"{ID}Z does not match Z in OEOS06 files you are using")
        if z != m.zz[mf]:
            _z_mismatch(mf, xh, t6, r, z)          # go to 66 (stops)
        # A quadratic through the tables mf, mg, mh.
        is_ = 0
        iw = 1
        for ir in range(l1, l1 + iq + 1):
            for it in range(k1, k1 + ip + 1):
                if mf2 == 1:
                    a.esk[it, ir] = e.epl[mf, it, ir]
                    continue                       # go to 46
                a.esk[it, ir] = quadeos06(
                    is_, iw, xh, e.epl[mf, it, ir], e.epl[mg, it, ir],
                    e.epl[mh, it, ir], e.xx[mf], e.xx[mg], e.xx[mh])
                if a.esk[it, ir] > 1.0e20:
                    io.write(io.ishort, f"{ID} problem it,ir,l3,k3,iq,ip="
                             f"{it:5d}{ir:5d}{m.l3:5d}{m.k3:5d}"
                             f"{iq:5d}{ip:5d}")
                    io.write(io.ishort, "".join(
                        f"{e.epl[ms, it, ir]:12.4e}"
                        for ms in range(mf, mf + 3)))
                is_ = 1
        if mi == mf2:
            # A second quadratic through mg, mh, mi, blended linearly.
            is_ = 0
            iw = 1
            dixr = (e.xx[mh] - xh) * a.dfsx[mh]
            for ir in range(l1, l1 + iq + 1):          # loop 47
                for it in range(k1, k1 + ip + 1):
                    a.esk2[it, ir] = quadeos06(
                        is_, iw, xh, e.epl[mg, it, ir], e.epl[mh, it, ir],
                        e.epl[mi, it, ir], e.xx[mg], e.xx[mh], e.xx[mi])
                    # Question (see the page): tests esk, not esk2.
                    if a.esk[it, ir] > 1.0e20:
                        io.write(io.ishort, f"{ID} problem it,ir,l3,k3,iq,ip="
                                 f"{it:5d}{ir:5d}{m.l3:5d}{m.k3:5d}"
                                 f"{iq:5d}{ip:5d}")
                        io.write(io.ishort, "".join(
                            f"{e.epl[ms, it, ir]:12.4e}"
                            for ms in range(mg, mg + 3)))
                    a.esk[it, ir] = (a.esk[it, ir] * dixr
                                     + a.esk2[it, ir] * (1.0 - dixr))
                    is_ = 1
        is_ = 0
        # Interpolate in T6 and rho (overlapping quadratics on 4x4).
        t6rinteos06(slr, slt)
        m.eos[iv] = m.esact

    # Unscale: P was tabulated as P/(T6 rho), E as E/T6.
    p0 = t6 * r
    m.eos[m.iri[1]] = m.eos[m.iri[1]] * p0
    m.eos[m.iri[2]] = m.eos[m.iri[2]] * t6
    tmass, moles, eground, fracz = gmass06(xh, z, moles, 0.0, 0.0, frac)
    if irad == 1:
        # Question (see the page): irad is always 1 here.
        radsub06(irad, t6, r, moles, tmass)
    else:
        m.eos[m.iri[5]] = m.eos[m.iri[5]] * moles * APROP / tmass
    return False


def _out_of_table(slt, slr):
    """Label 62: write the message and take the alternate return."""
    a = aeos06_mod
    io = iounits_mod
    io.write(io.ishort, f" {ID} T6/LogR outside of table range (62)")
    io.write(io.ishort, f" slt, t6a(1),t6a(nt): {slt} {a.t6a[1]} {a.t6a[NT]}")
    io.write(io.ishort, f" slr, rho(1), rho(nr): {slr} {a.rho[1]} {a.rho[NR]}")
    return True


def _empty_region(xh, t6, r):
    """Label 65: write the message and take the alternate return."""
    io = iounits_mod
    io.write(io.ishort, f" {ID}T6/log rho in empty region of table (65)")
    io.write(io.ishort, f"xh,t6,r={xh:12.4e}{t6:12.4e}{r:12.4e}")
    return True


def _z_mismatch(mf, xh, t6, r, z):
    """Label 66: the tables' Z differs from Z; write and stop."""
    m = eos06misc_mod
    io = iounits_mod
    io.write(io.ishort, f" {ID} Z does not match Z in EOSdata* files you "
             "are using (66)")
    io.write(io.ishort, f"mf,zz(mf)={mf:5d}{m.zz[mf]:12.4e}")
    # Question (see the page): ip, iq printed under the heading iq, ip.
    io.write(io.ishort, f"  iq,ip,k3,l3,xh,t6,r,z= {m.ip:5d}{m.iq:5d}"
             f"{m.k3:5d}{m.l3:5d}{xh:12.4e}{t6:12.4e}{r:12.4e}{z:12.4e}")
    raise SystemExit("ESAC06: Z does not match Z in EOSdata* files (66)")


def seteos06():
    """One-time setup of the OPAL 2006 EOS, called by setups when lopale06
    is set: map index into iri, copy the X values into xx, and read the
    table (the file is opened by parmin).

    Translates: esac06.f#SETEOS06 (09c3b9b)
    """
    m = eos06misc_mod
    for i in range(1, 11):
        for j in range(1, 11):
            if m.index[i] == j:
                m.iri[i] = j
    for i in range(1, MX + 1):
        eeeos06_mod.xx[i] = aeos06_mod.xa[i]
    readcoeos06()
