"""Python version of oldrot_mod.f (YREC 09c3b9b): the rotation of the model
at the start of the step, formerly COMMON/OLDROT/.

Translates: oldrot_mod.f#OLDROT_MOD (09c3b9b)

Routines that USE OLDROT_MOD import this module and write
``oldrot_mod.<name>``. Saved by hpoint (wold also by starin).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

wold = np.zeros(_JSON + 1)      # Omega
hjx = np.zeros(_JSON + 1)       # J/M
hio = np.zeros(_JSON + 1)       # moment of inertia; read nowhere
hgo = np.zeros(_JSON + 1)       # gravity (read by midmod)
r0x = np.zeros(_JSON + 1)       # distortion term r0
eta2x = np.zeros(_JSON + 1)     # distortion term eta2
