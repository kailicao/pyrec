"""Python version of alot06all_mod.f (YREC 09c3b9b): all the low-temperature opacity tables, formerly COMMON/ALOT06ALL/.

Translates: alot06all_mod.f#ALOT06ALL_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ALOT06ALL_MOD import this module and write ``alot06all_mod.<name>``.
"""

import numpy as np

_NUMD = 19                         # PRIVATE in the Fortran
_NUMT = 85
_NUMXZ = 143

# The 143 compositions' tables (readalex06).
captalex06 = np.zeros((_NUMXZ + 1, _NUMT + 1, _NUMD + 1))
