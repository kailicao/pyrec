"""Python version of deuter_mod.f (YREC 09c3b9b): the deuterium burning
rates, the accreted mass and the base of the convection zone, formerly
COMMON/DEUTER/.

Translates: deuter_mod.f#DEUTER_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE DEUTER_MOD import this module and write ``deuter_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

drate = np.zeros(_JSON + 1)     # deuterium burning rate per shell (deutrate)
drate0 = np.zeros(_JSON + 1)    # the same at the start of the step
fmassacc = 0.0                  # mass accreted in the step (mdot)
jcz = 0                         # shell at the base of the surface CZ (main, massloss)
