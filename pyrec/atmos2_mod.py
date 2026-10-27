"""Python version of atmos2_mod.f (YREC 09c3b9b): the surface-pressure table for KTTAU = 3 or 4, formerly COMMON/ATMOS2/.

Translates: atmos2_mod.f#ATMOS2_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ATMOS2_MOD import this module and write ``atmos2_mod.<name>``.
"""

import numpy as np

_NG = 11                           # PRIVATE in the Fortran
_NT = 57

# All read by setups from the file FATM.
atmpl = np.zeros((_NT + 1, _NG + 1))   # log P
atmtl = np.zeros(_NT + 1)          # log Teff grid
atmgl = np.zeros(_NG + 1)          # log g grid
atmz = 0.0                         # the table's Z (read by no routine)
