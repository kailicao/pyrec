"""Python version of llot95_mod.f (YREC 09c3b9b): the OPAL95 opacity tables
at the current Z, formerly COMMON/LLOT95/.

Translates: llot95_mod.f#LLOT95_MOD (09c3b9b)

Routines that USE LLOT95_MOD import this module and write
``llot95_mod.<name>``. Arrays keep Fortran's index order, index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_NUMD = 19
_NUMT = 70
_NUMX = 10

cappa = np.zeros((_NUMX + 1, _NUMT + 1, _NUMD + 1))   # log kappa at Z = ztab (op95ztab)
ztab = 0.0                                          # the Z of cappa (op95ztab)
