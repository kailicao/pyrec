"""Python version of gravez_mod.f (YREC 09c3b9b): the metals and their
coefficients on grsett's equally spaced grid, formerly COMMON/GRAVEZ/.

Translates: gravez_mod.f#GRAVEZ_MOD (09c3b9b)

Routines that USE GRAVEZ_MOD import this module and write
``gravez_mod.<name>``. "_h" marks the grid's midpoints.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

ecod1z = np.zeros(_JSON + 1)        # settling coefficient at the points (model_to_equal)
ecod1z_h = np.zeros(_JSON + 1)      # the same at the midpoints (model_to_equal, grsett)
ecod2z_h = np.zeros(_JSON + 1)      # diffusion coefficient (model_to_equal); read nowhere
eqcod1z_h = np.zeros(_JSON + 1)     # abundance-dependent part of ecod1z_h
eqcod2z_h = np.zeros(_JSON + 1)     # the same for ecod2z_h; read nowhere
ez = np.zeros(_JSON + 1)            # Z at the points, then its change (grsett)
ez_h = np.zeros(_JSON + 1)          # Z at the midpoints (model_to_equal, lax_wendrof1)
