"""Python version of difad_mod.f (YREC 7be4133): the advective and diffusive
coefficients of angular momentum transport on the even grid, formerly
COMMON/DIFAD/.

Translates: difad_mod.f#DIFAD_MOD (7be4133)

Routines that USE DIFAD_MOD import this module and write
``difad_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–29.
ecod3 = np.zeros(_JSON + 1)     # advective coefficient (rotgrid, dadcoeft)
ecod4 = np.zeros(_JSON + 1)     # diffusive coefficient (rotgrid, dadcoeft)
