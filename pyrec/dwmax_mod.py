"""Python version of dwmax_mod.f (YREC 09c3b9b): the critical rotation
gradient of the shear instability, formerly COMMON/DWMAX/.

Translates: dwmax_mod.f#DWMAX_MOD (09c3b9b)

Routines that USE DWMAX_MOD import this module and write
``dwmax_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

qwrmax = np.zeros(_JSON + 1)    # critical dOmega/dln r per interface (physic)
qwrmx0 = np.zeros(_JSON + 1)    # its start-of-step copy (hpoint; read nowhere)
