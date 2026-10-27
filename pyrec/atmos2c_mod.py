"""Python version of atmos2c_mod.f (YREC 09c3b9b): the surface-pressure table for KTTAU = 5, formerly COMMON/ATMOS2C/.

Translates: atmos2c_mod.f#ATMOS2C_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ATMOS2C_MOD import this module and write ``atmos2c_mod.<name>``.
"""

import numpy as np

_NGC = 11                          # PRIVATE in the Fortran
_NTC = 76

# All read by setups from the file FATM.
atmplc = np.zeros((_NTC + 1, _NGC + 1))   # log P
atmtlc = np.zeros(_NTC + 1)        # log Teff grid
atmglc = np.zeros(_NGC + 1)        # log g grid
