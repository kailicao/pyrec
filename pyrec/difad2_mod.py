"""Python version of difad2_mod.f (YREC 09c3b9b): the advective and
diffusive parts of the circulation velocity, formerly COMMON/DIFAD2/.

Translates: difad2_mod.f#DIFAD2_MOD (09c3b9b)

Routines that USE DIFAD2_MOD import this module and write
``difad2_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

vesa = np.zeros(_JSON + 1)      # advective part (vcirc)
vesa0 = np.zeros(_JSON + 1)     # its previous value (vcirc)
vesd = np.zeros(_JSON + 1)      # diffusive part (vcirc)
vesd0 = np.zeros(_JSON + 1)     # its previous value (vcirc)
