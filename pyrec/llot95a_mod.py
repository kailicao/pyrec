"""Python version of llot95a_mod.f (YREC 09c3b9b): all the OPAL95 opacity
tables and their grids, formerly COMMON/LLOT95A/.

Translates: llot95a_mod.f#LLOT95A_MOD (09c3b9b)

Routines that USE LLOT95A_MOD import this module and write
``llot95a_mod.<name>``. Arrays keep Fortran's index order, index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_NUMD = 19
_NUMT = 70
_NUMX = 10
_NUMXZ = 126
_NUMZ = 13

tgr = np.zeros(_NUMT + 1)          # log T grid (ll95tbl)
# X values of the tables (initializer).
xxg = np.array([0.0, 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0])
r0gr = np.zeros(_NUMD + 1)         # log R grid (ll95tbl)
# Z values of the tables (initializer).
zzg = np.array([0.0, 0.0, 0.0001, 0.0003, 0.001, 0.002, 0.004, 0.01, 0.02,
                0.03, 0.04, 0.06, 0.08, 0.10])
cappa2 = np.zeros((_NUMXZ + 1, _NUMT + 1, _NUMD + 1))   # log kappa, every table (ll95tbl)
# Number of X tables per Z (initializer; read nowhere).
numxm = np.array([0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 9, 9, 8])
# Offset of each Z's tables in cappa2's first index (initializer).
nz = np.array([0, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 109, 118])
