"""Python version of tridi_mod.f (YREC 7be4133): a tridiagonal system and its solution, formerly COMMON/TRIDI/.

Translates: tridi_mod.f#TRIDI_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE TRIDI_MOD import this module and write ``tridi_mod.<name>``.
The system is a[i] u[i-1] + b[i] u[i] + c[i] u[i+1] = d[i]. Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–33.
a = np.zeros(_JSON + 1)      # sub-diagonal (dcoeft, ccoeft, getqua)
b = np.zeros(_JSON + 1)      # diagonal
c = np.zeros(_JSON + 1)      # super-diagonal
d = np.zeros(_JSON + 1)      # right-hand side
u = np.zeros(_JSON + 1)      # solution (tridia, ctridi, getqua)
gama = np.zeros(_JSON + 1)   # Thomas algorithm's work array
