"""Python version of splin_mod.f (YREC 09c3b9b): work arrays for interpolation with osplin, formerly COMMON/SPLIN/.

Translates: splin_mod.f#SPLIN_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE SPLIN_MOD import this module and write ``splin_mod.<name>``.
Scratch space of eq2mod, mixcom, mixgrid and rotgrid. Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

xval = np.zeros(_JSON + 1)    # abscissas of evaluation
yval = np.zeros(_JSON + 1)    # values returned by osplin
xtab = np.zeros(_JSON + 1)    # the table's abscissas
ytab = np.zeros(_JSON + 1)    # the table's values
