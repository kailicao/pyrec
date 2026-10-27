"""Python version of oldmod_mod.f (YREC 09c3b9b): the model at the start of
the step, formerly COMMON/OLDMOD/.

Translates: oldmod_mod.f#OLDMOD_MOD (09c3b9b)

Routines that USE OLDMOD_MOD import this module and write
``oldmod_mod.<name>``. Saved by main, hpoint and starin (hpoint also uses
hpo, hto, hdo as scratch arrays).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

hpo = np.zeros(_JSON + 1)               # log P
hto = np.zeros(_JSON + 1)               # log T
hro = np.zeros(_JSON + 1)               # log r
hlo = np.zeros(_JSON + 1)               # L
hdo = np.zeros(_JSON + 1)               # log rho
hcompp = np.zeros((16, _JSON + 1))      # composition, hcompp[k, j]
hss = np.zeros(_JSON + 1)               # mass coordinates
lco = np.zeros(_JSON + 1, dtype=bool)   # convective flags (read by midmod)
lczo = np.zeros(_JSON + 1, dtype=bool)  # read nowhere
tefflo = 0.0                            # log Teff
mo = 0                                  # number of shells
