"""Python version of pualpha_mod.f (YREC 09c3b9b): mixing-length quantities for the pulsation output, formerly COMMON/PUALPHA/.

Translates: pualpha_mod.f#PUALPHA_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PUALPHA_MOD import this module and write ``pualpha_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# At the current point: set to 0 by tpgrad and computed nowhere (see
# TPGRAD's page); written by envint.
alfmlt = 0.0
phmlt = 0.0
cmxmlt = 0.0
# At each shell (coefft, when lpulse); written by wrtmod. Index 0 unused.
valfmlt = np.zeros(_JSON + 1)
vphmlt = np.zeros(_JSON + 1)
vcmxmlt = np.zeros(_JSON + 1)
