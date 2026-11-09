"""Python version of quadd_mod.f (YREC 7be4133): terms of the rotating
star's potential and the circulation's quadrupole factor, formerly
COMMON/QUADD/.

Translates: quadd_mod.f#QUADD_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE QUADD_MOD import this module and write ``quadd_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Index 0 unused (Fortran lines 28–31).
# phisp, phirot and phidis are set by fpft and read by no routine (unused?).
phisp = np.zeros(_JSON + 1)     # GM/r
phirot = np.zeros(_JSON + 1)    # Omega^2 r^2
phidis = np.zeros(_JSON + 1)    # distortion term
rat = np.zeros(_JSON + 1)       # quadrupole factor of the circulation (vcirc)
