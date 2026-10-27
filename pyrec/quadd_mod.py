"""Python version of quadd_mod.f (YREC 09c3b9b): terms of the rotating star's potential and the circulation's quadrupole factor, formerly COMMON/QUADD/.

Translates: quadd_mod.f#QUADD_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE QUADD_MOD import this module and write ``quadd_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Index 0 unused.
# Question (see the page): set by fpft and read by no routine.
phisp = np.zeros(_JSON + 1)     # GM/r
phirot = np.zeros(_JSON + 1)    # Omega^2 r^2
phidis = np.zeros(_JSON + 1)    # distortion term
rat = np.zeros(_JSON + 1)       # quadrupole factor of the circulation (vcirc)
