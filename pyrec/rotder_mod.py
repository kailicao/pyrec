"""Python version of rotder_mod.f (YREC 09c3b9b): derivatives at each shell for the rotational mixing, formerly COMMON/ROTDER/.

Translates: rotder_mod.f#ROTDER_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ROTDER_MOD import this module and write ``rotder_mod.<name>``.
Set by coefft (when lrot), read by setupv.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

sqot = np.zeros(_JSON + 1)      # dln kappa / dln T, index 0 unused
sqet = np.zeros(_JSON + 1)      # dln epsilon / dln T
fracnuc = np.zeros(_JSON + 1)   # nuclear share of the energy generation
