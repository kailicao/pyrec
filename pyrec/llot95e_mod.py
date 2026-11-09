"""Python version of llot95e_mod.f (YREC 7be4133): the OPAL95 opacity table
at the surface X and Z, formerly COMMON/LLOT95E/ (XTAB now XOP95).

Translates: llot95e_mod.f#LLOT95E_MOD (7be4133)

Routines that USE LLOT95E_MOD import this module and write
``llot95e_mod.<name>``. Index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran) (Fortran lines 28–29).
_NUMD = 19
_NUMT = 70

# The table at X = xop95 (Fortran lines 31–32).
cappae = np.zeros((_NUMT + 1, _NUMD + 1))   # log kappa at X = xop95, Z = ztab (op95xtab)
xop95 = 0.0                                 # the X of cappae (op95xtab)
