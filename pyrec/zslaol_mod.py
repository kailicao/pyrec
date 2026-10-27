"""Python version of zslaol_mod.f (YREC 09c3b9b): the density splines of the pure-metal LAOL89 opacity table, formerly COMMON/ZSLAOL/.

Translates: zslaol_mod.f#ZSLAOL_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ZSLAOL_MOD import this module and write ``zslaol_mod.<name>``.
Set by zsulaol, read by gtpurz. Index 0 unused in every dimension.
"""

import numpy as np

zslaoll = np.zeros((105, 53))   # log kappa at the nonzero entries
zsrhol = np.zeros((105, 53))    # log rho there
zsdorl = np.zeros((105, 53))    # the spline's second derivatives (cspline)
numrsz = np.zeros(53, dtype=np.int32)   # entries per temperature (0: no spline)
