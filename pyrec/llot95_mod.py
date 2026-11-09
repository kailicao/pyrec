"""Python version of llot95_mod.f (YREC 7be4133): the OPAL95 opacity tables
at the current Z, formerly COMMON/LLOT95/.

Translates: llot95_mod.f#LLOT95_MOD (7be4133)

Routines that USE LLOT95_MOD import this module and write
``llot95_mod.<name>``. Arrays keep Fortran's index order, index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran) (Fortran lines 28–30).
_NUMD = 19
_NUMT = 70
_NUMX = 10

# The tables at Z = ztab (Fortran lines 32–33).
cappa = np.zeros((_NUMX + 1, _NUMT + 1, _NUMD + 1))   # log kappa at Z = ztab (op95ztab)
ztab = 0.0                                          # the Z of cappa (op95ztab)
