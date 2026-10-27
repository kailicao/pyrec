"""Python version of eqbound06.f (YREC 09c3b9b): whether (T, rho) lies in
the OPAL 2006 equation-of-state table, and the blending weight near its
edges."""

from pyrec import aeos06_mod
from pyrec import eos06misc_mod

# PARAMETER (mx=5, mv=10, nr=169, nt=197): the table's dimensions.
MX = 5
MV = 10
NR = 169
NT = 197


def eqbound06(t, dl, fac, ltab, lramp):
    """Test whether (T, log rho) lies inside the OPAL 2006 table, following
    its irregular edge (nta, nra), and give the weight fac = f_rho * f_T for
    blending OPAL with the Yale/SCV equation of state near the edges.

    Translates: eqbound06.f#EQBOUND06 (09c3b9b)

    Arguments: t (K), dl (log10 rho); fac, ltab, lramp (set only).

    Returns (fac, ltab, lramp): the weight (0 out of table), whether the
    point is in the table, and whether ramping is on (True out of table).
    """
    rho = aeos06_mod.rho
    t6a = aeos06_mod.t6a
    nta = eos06misc_mod.nta
    nra = eos06misc_mod.nra
    t6 = t * 1.0e-6
    d = 10.0**dl

    while True:  # one pass; break goes to the out-of-table exit, label 9999
        if d < rho[1] or d >= rho[NR]:
            break  # out of table in density
        if t6 > t6a[1] or t6 <= t6a[NT]:
            break  # out of table in temperature

        drampf = 1.0
        trampf = 1.0
        ltab = True
        lramp = False

        # Ramps at the easy borders, linear in T6 and rho.
        # Question (see the page): linear, where the 1995 and 2001
        # versions ramp in log rho.
        if t6 > t6a[2]:
            trampf = (t6 - t6a[1]) / (t6a[2] - t6a[1])
            lramp = True
        if d < rho[2]:
            drampf = (d - rho[1]) / (rho[2] - rho[1])
            lramp = True

        # Find k, l by a linear search from ESAC06's last indices.
        # Question (see the page): k3, l2 are left by ESAC06's last call.
        k = eos06misc_mod.k3
        while t6 <= t6a[k]:  # label 10
            k = k + 1
        while t6 > t6a[k - 1]:  # label 20
            k = k - 1
        l = eos06misc_mod.l2
        while d > rho[l]:  # label 30
            l = l + 1
        while d <= rho[l - 1]:  # label 40
            l = l - 1

        # The irregular edge.
        if l > nra[k]:
            break  # out of table in density
        if l == nra[k]:
            lramp = True
            drampf = (rho[l] - d) / (rho[l] - rho[l - 1])
        if k > nta[l]:
            break  # out of table in temperature
        if k == nta[l]:
            lramp = True
            # Question (see the page): replaces, not combines, the
            # high-temperature ramp above.
            trampf = (t6 - t6a[k]) / (t6a[k - 1] - t6a[k])

        fac = drampf * trampf
        return fac, ltab, lramp

    # label 9999: out of table
    ltab = False
    lramp = True
    fac = 0.0
    return fac, ltab, lramp
