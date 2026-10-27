"""Python version of aeos06_mod.f (YREC 09c3b9b): the OPAL 2006 equation-of-state tables and grids, formerly COMMON/AEOS06/.

Translates: aeos06_mod.f#AEOS06_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE AEOS06_MOD import this module and write ``aeos06_mod.<name>``.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_MX = 5
_MV = 10
_NR = 169
_NT = 197

xz = np.zeros((_MX + 1, _MV + 1, _NT + 1, _NR + 1))   # the tables (readcoeos06)
t6list = np.zeros((_NR + 1, _NT + 1))
rho = np.zeros(_NR + 1)            # density grid (readcoeos06)
t6a = np.zeros(_NT + 1)            # T6 grid (readcoeos06)
esk = np.zeros((_NT + 1, _NR + 1))   # values interpolated in X (esac06)
esk2 = np.zeros((_NT + 1, _NR + 1))
dfsx = np.zeros(_MX + 1)           # inverse spacings in X, T6, rho
dfs = np.zeros(_NT + 1)
dfsr = np.zeros(_NR + 1)
midx = 0                           # formerly M (readcoeos06)
mf = 0                             # the X index (esac06)
# The X values of the tables, index 0 unused (DATA, from esac06; the
# Fortran's single-precision literals, as stored in REAL*8).
xa = np.array([0.0, 0.0, np.float32(0.2), np.float32(0.4),
               np.float32(0.6), np.float32(0.8)], dtype=float)
