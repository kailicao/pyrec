"""Python version of deuter_mod.f (YREC 7be4133): the deuterium burning
rates, the accreted mass and the base of the convection zone, formerly
COMMON/DEUTER/.

Translates: deuter_mod.f#DEUTER_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE DEUTER_MOD import this module and write ``deuter_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–31. FMASSACC has no initial value in the Fortran (0
# here); MDOT sets it only when mass is gained (see dburn's page).
drate = np.zeros(_JSON + 1)     # deuterium burning rate per shell (deutrate)
drate0 = np.zeros(_JSON + 1)    # the same at the start of the step
fmassacc = 0.0                  # mass accreted in the step, g (mdot)
jcz = 0                         # shell at the base of the surface CZ (main, massloss)
