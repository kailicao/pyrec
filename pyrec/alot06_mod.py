"""Python version of alot06_mod.f (YREC 7be4133): the low-temperature opacity
table for the current composition, formerly COMMON/ALOT06/.

Translates: alot06_mod.f#ALOT06_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ALOT06_MOD import this module and write ``alot06_mod.<name>``.
"""

import numpy as np

_NUMD = 19                         # PRIVATE in the Fortran (lines 34–35)
_NUMT = 85

# The table for the current composition (Fortran line 37).
capalex06 = np.zeros((_NUMT + 1, _NUMD + 1))   # log kappa (alex06tab)

# Initializers, from readalex06's DATA (Fortran lines 38–42).
xe = 0.0                           # composition of the current table
ze = 0.0
ix = 4                             # first X table used (alex06tab)
it = 43                            # starting indices of getalex06's search
id = 10
