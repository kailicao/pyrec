"""Python version of op95interp_mod.f (YREC 09c3b9b): the interpolation stencil and weights for the OPAL 1995 opacity tables, formerly COMMON/OP95INDX/, /OP95FACT/ and /OP95EXT/.

Translates: op95interp_mod.f#OP95INTERP_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE OP95INTERP_MOD import this module and write ``op95interp_mod.<name>``.
All are set by getopal95 and read by op952d, op953d, op954d.
"""

import numpy as np

# First table index of each four-point stencil (initializers, formerly
# ll95tbl's DATA): in Z, in X for each Z, in T, in R for each T.
jz = 1
jx = np.ones((5, 5), dtype=np.int32)    # index 0 unused in both dimensions
jt = 1
jd = np.ones(5, dtype=np.int32)
# Four-point weights (bz, bx from intrp2; bt, bd from interp) and the
# derivative weights in T and R (btt, bdd).
bz = np.zeros(5)
bx = np.zeros((5, 5))
bt = np.zeros(5)
btt = np.zeros(5)
bd = np.zeros((5, 5))
bdd = np.zeros((5, 5))
rho_t63l = 0.0          # log R = log rho - 3 (log T - 6), formerly RL
rho_t63ll = 0.0         # lowest log R of the table, formerly RLL
rlh = np.zeros(5)       # table edge in log R at each T (extrapolation)
lextl = False           # extrapolation below the table in R
lexth = False           # extrapolation above it (high-density corner)
lexh = np.zeros(5, dtype=bool)   # the same, for each of the four T
