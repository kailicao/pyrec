"""Python version of burn_mod.f (YREC 09c3b9b): a copy of the composition passed between the mixing routines, formerly COMMON/BURN/.

Translates: burn_mod.f#BURN_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE BURN_MOD import this module and write ``burn_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

# hcompm[k, j] for HCOMPM(K, J): set by mix (and hpoint), read by rotmix.
hcompm = np.zeros((16, dimensions_mod.json + 1))
