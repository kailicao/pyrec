"""Python version of debhu_mod.f (YREC 7be4133): the inputs of the
Debye-Hueckel correction, formerly COMMON/DEBHU/.

Translates: debhu_mod.f#DEBHU_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case (REAL*8).
Routines that USE DEBHU_MOD import this module and write
``debhu_mod.<name>``.
"""

import numpy as np

# Fortran lines 24–29.
cdh = 0.0               # Debye-Hueckel constant (setups)
zdh = np.zeros(19)      # abundances of the 18 heavy elements, index 0 unused
xxdh = 0.0              # local X (copied by the callers of the EOS)
yydh = 0.0              # local Y
zzdh = 0.0              # local Z
dhnue = np.zeros(19)    # charges of the 18 heavy elements (setups), index 0 unused
