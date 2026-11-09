"""Python version of egridder_mod.f (YREC 7be4133): the higher-derivative
coefficients of the diffusion-plus-advection treatment, formerly
COMMON/EGRIDDER/.

Translates: egridder_mod.f#EGRIDDER_MOD (7be4133)

Routines that USE EGRIDDER_MOD import this module and write
``egridder_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–31.
eqqcod = np.zeros(_JSON + 1)    # qqcod * Omega (vcirc); read by dadcoeft
# Question (see the page): eqqqcod is read by dadcoeft but set nowhere.
eqqqcod = np.zeros(_JSON + 1)
qqcod = np.zeros(_JSON + 1)     # setupv; read by vcirc
qqqcod = np.zeros(_JSON + 1)    # setupv; read nowhere
