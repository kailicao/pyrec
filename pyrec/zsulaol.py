"""Python version of zsulaol.f (YREC 09c3b9b): splines of the pure-Z LAOL89
opacity table in density."""

import math

import numpy as np

from pyrec import zlaol_mod
from pyrec import zslaol_mod
from pyrec.cspline import cspline

# The Fortran's local work arrays (SAVE).
top = np.zeros(105)
tr = np.zeros(105)
tdor = np.zeros(105)


def zsulaol():
    """Convert the pure-Z table's temperatures to log10 T in place, and
    for each temperature fit a natural cubic spline of log kappa against
    log rho through the nonzero entries, if there are at least four
    (otherwise numrsz is 0 and gtpurz skips that temperature).

    Translates: zsulaol.f#ZSULAOL (09c3b9b)

    No arguments; reads zlaol_mod (from rdzlaol), sets zslaol_mod. None is
    returned.
    """
    z = zlaol_mod
    zs = zslaol_mod
    # Question (see the page): the logarithm is taken in place.
    for it in range(1, z.numtz + 1):
        z.zot[it] = math.log10(z.zot[it])
    for it in range(1, z.numtz + 1):
        nums = 0
        for ir in range(1, z.numrhoz + 1):
            zs.zslaoll[ir, it] = 0.0
            zs.zsrhol[ir, it] = 0.0
            zs.zsdorl[ir, it] = 0.0
            if z.zolaol[ir, it] != 0.0:
                nums = nums + 1
                top[nums] = math.log10(z.zolaol[ir, it])
                tr[nums] = math.log10(z.zorho[ir])
        if nums >= 4:
            zs.numrsz[it] = nums
            cspline(tr, top, nums, 1.0e30, 1.0e30, tdor)
            for ir in range(1, nums + 1):
                zs.zslaoll[ir, it] = top[ir]
                zs.zsrhol[ir, it] = tr[ir]
                zs.zsdorl[ir, it] = tdor[ir]
        else:
            zs.numrsz[it] = 0
