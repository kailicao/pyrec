"""Python version of galot06_mod.f (YREC 09c3b9b): the grids of the
Alexander et al. 2006 opacity tables, formerly COMMON/GALOT06/.

Translates: galot06_mod.f#GALOT06_MOD (09c3b9b)

Routines that USE GALOT06_MOD import this module and write
``galot06_mod.<name>``. Arrays have index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_NUMD = 19
_NUMT = 85
_NUMX = 9
_NUMZ = 16

# log T grid (initializer).
tgr06 = np.array([0.0] + [
    2.70, 2.75, 2.80, 2.85, 2.90, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97,
    2.98, 2.99, 3.00, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09,
    3.10, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.20, 3.21,
    3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.30, 3.31, 3.32, 3.33,
    3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.40, 3.41, 3.42, 3.43, 3.44, 3.45,
    3.46, 3.47, 3.48, 3.49, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85,
    3.90, 3.95, 4.00, 4.05, 4.10, 4.15, 4.20, 4.25, 4.30, 4.35, 4.40, 4.45,
    4.50])
# X values of the tables.
xxg06 = np.array([0.0, 0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 1.0])
# log R grid.
r0gr06 = np.array([0.0, -8.0, -7.5, -7.0, -6.5, -6.0, -5.5, -5.0, -4.5,
                   -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5,
                   1.0])
# Z values of the tables.
zzg06 = np.array([0.0, 0.0, 0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.002,
                  0.004, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10])

assert len(tgr06) == _NUMT + 1 and len(xxg06) == _NUMX + 1
assert len(r0gr06) == _NUMD + 1 and len(zzg06) == _NUMZ + 1
