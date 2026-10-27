"""Python version of envstruct_mod.f (YREC 09c3b9b): the envelope's
structure, point by point, formerly COMMON/ENVSTRUCT/.

Translates: envstruct_mod.f#ENVSTRUCT_MOD (09c3b9b)

Routines that USE ENVSTRUCT_MOD import this module and write
``envstruct_mod.<name>``. Filled by envint (point 1 at the bottom after it
reverses them); shifted by getnewenv and starin.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

envp = np.zeros(_JSON + 1)              # log P
envt = np.zeros(_JSON + 1)              # log T
envs = np.zeros(_JSON + 1)              # log(m/M)
envd = np.zeros(_JSON + 1)              # log rho
envr = np.zeros(_JSON + 1)              # log r
envx = np.zeros(_JSON + 1)              # X
envz = np.zeros(_JSON + 1)              # Z
lcenv = np.zeros(_JSON + 1, dtype=bool) # convective point
edels = np.zeros((4, _JSON + 1))        # the three gradients, edels[k, j]
evels = np.zeros(_JSON + 1)             # convective velocity
ebetas = np.zeros(_JSON + 1)            # beta
egam1 = np.zeros(_JSON + 1)             # Gamma_1
eqcp = np.zeros(_JSON + 1)              # thermodynamic derivative QQCP
efxions = np.zeros((4, _JSON + 1))      # ionization fractions, efxions[k, j]
envo = np.zeros(_JSON + 1)              # opacity
envl = np.zeros(_JSON + 1)              # luminosity
eqdt = np.zeros(_JSON + 1)              # thermodynamic derivative QQDT
numenv = 0                              # number of stored points
