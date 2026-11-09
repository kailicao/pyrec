"""Python version of fac_mod.f (YREC 7be4133): the valid gravity ranges of
the surface-pressure tables, formerly COMMON/FAC/.

Translates: fac_mod.f#FAC_MOD (7be4133)

Routines that USE FAC_MOD import this module and write ``fac_mod.<name>``.
All INTEGER*4; arrays with index 0 unused.
"""

import numpy as np

# Table dimensions, PRIVATE in the Fortran (Fortran lines 28–29).
_NT = 57
_NTC = 76

# Fortran lines 31–36.
imin = np.zeros(_NT + 1, dtype=int)         # lowest valid gravity per T (setups)
iminmax = np.zeros(_NT + 1, dtype=int)      # highest valid gravity per T (setups)
jp = 0                                      # last T index (surfp, kcsurfp); read nowhere
kp = np.zeros(5, dtype=int)                 # last gravity indices; read nowhere
imin2 = np.zeros(_NTC + 1, dtype=int)       # the same for the 76-T table (kcsurfp)
iminmax2 = np.zeros(_NTC + 1, dtype=int)
