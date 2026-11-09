"""Python version of oldrot_mod.f (YREC 7be4133): the rotation of the model
at the start of the step, formerly COMMON/OLDROT/.

Translates: oldrot_mod.f#OLDROT_MOD (7be4133)

Routines that USE OLDROT_MOD import this module and write
``oldrot_mod.<name>``. Saved by hpoint (wold also by starin).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–33.
wold = np.zeros(_JSON + 1)      # Omega
hjx = np.zeros(_JSON + 1)       # J/M
hio = np.zeros(_JSON + 1)       # moment of inertia; read nowhere (unused?)
hgo = np.zeros(_JSON + 1)       # gravity (read by midmod)
r0x = np.zeros(_JSON + 1)       # distortion term r0
eta2x = np.zeros(_JSON + 1)     # distortion term eta2
