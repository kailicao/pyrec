"""Python version of oldmod_mod.f (YREC 7be4133): the model at the start of
the step, formerly COMMON/OLDMOD/.

Translates: oldmod_mod.f#OLDMOD_MOD (7be4133)

Routines that USE OLDMOD_MOD import this module and write
``oldmod_mod.<name>``. Saved by main, hpoint and starin (hpoint also uses
hpo, hto, hdo as scratch arrays).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# The old model, its composition and flags (Fortran lines 28–41).
hpo = np.zeros(_JSON + 1)               # log P
hto = np.zeros(_JSON + 1)               # log T
hro = np.zeros(_JSON + 1)               # log r
hlo = np.zeros(_JSON + 1)               # L
hdo = np.zeros(_JSON + 1)               # log rho
hcompp = np.zeros((16, _JSON + 1))      # composition, hcompp[k, j]
hss = np.zeros(_JSON + 1)               # mass coordinates
lco = np.zeros(_JSON + 1, dtype=bool)   # convective flags (read by midmod)
lczo = np.zeros(_JSON + 1, dtype=bool)  # read by crrect (for getrot)
tefflo = 0.0                            # log Teff
hlsuro = 0.0                            # surface L (for entime)
