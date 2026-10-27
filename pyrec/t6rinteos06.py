"""Python version of t6rinteos06.f (YREC 09c3b9b): interpolation of one
OPAL 2006 equation-of-state quantity in temperature and density."""

from pyrec import aeos06_mod
from pyrec import eos06misc_mod
from pyrec import iounits_mod
from pyrec.quadeos06 import quadeos06


def t6rinteos06(slr, slt):
    """Interpolate one EOS quantity at density slr and temperature slt on
    the 4 x 4 grid ESAC06 has chosen, by blending overlapping quadratics;
    the result goes to eos06misc_mod.esact.

    Translates: t6rinteos06.f#T6RINTEOS06 (09c3b9b)

    Arguments: slr, slt. Reads the grid indices K1-K4, L1-L4 and the
    switches IP, IQ from eos06misc_mod, the tables from aeos06_mod. Uses
    eos06misc_mod.h and .q as work arrays. None is returned.
    """
    a = aeos06_mod
    e = eos06misc_mod
    l1, l2, l3, l4 = e.l1, e.l2, e.l3, e.l4
    k1, k2, k3, k4 = e.k1, e.k2, e.k3, e.k4
    ip, iq = e.ip, e.iq
    iu = 0
    is_ = 0

    # Quadratics in density at each temperature: h through rows l1-l3, and,
    # if IQ = 3, q through l2-l4.
    for kx in range(k1, k1 + ip + 1):
        iw = 1
        iu = iu + 1
        e.h[iu] = quadeos06(is_, iw, slr, a.esk[kx, l1], a.esk[kx, l2],
                            a.esk[kx, l3], a.rho[l1], a.rho[l2], a.rho[l3])
        if iq == 3:
            iw = 2
            e.q[iu] = quadeos06(is_, iw, slr, a.esk[kx, l2], a.esk[kx, l3],
                                a.esk[kx, l4], a.rho[l2], a.rho[l3],
                                a.rho[l4])
        is_ = 1

    # Quadratics in temperature, blended if IP = 3.
    is_ = 0
    iw = 1
    e.esact = quadeos06(is_, iw, slt, e.h[1], e.h[2], e.h[3], a.t6a[k1],
                        a.t6a[k2], a.t6a[k3])
    if iq == 3:
        esactq = quadeos06(is_, iw, slt, e.q[1], e.q[2], e.q[3], a.t6a[k1],
                           a.t6a[k2], a.t6a[k3])
    if ip == 3:
        esact2 = quadeos06(is_, iw, slt, e.h[2], e.h[3], e.h[4], a.t6a[k2],
                           a.t6a[k3], a.t6a[k4])
        dix = (a.t6a[k3] - slt) * a.dfs[k3]
        e.esact = e.esact * dix + esact2 * (1.0 - dix)
        if iq == 3:
            esactq2 = quadeos06(is_, iw, slt, e.q[2], e.q[3], e.q[4],
                                a.t6a[k2], a.t6a[k3], a.t6a[k4])
            esactq = esactq * dix + esactq2 * (1.0 - dix)
    # The blend in density.
    if iq == 3:
        dix2 = (a.rho[l3] - slr) * a.dfsr[l3]
        if ip == 3:
            e.esact = e.esact * dix2 + esactq * (1 - dix2)

    if e.esact > 1.0e15:
        iounits_mod.write(iounits_mod.ishort,
                          "T6RINTEOS06: Interpolation indices out of "
                          "range;please report conditions.")
        raise SystemExit()
