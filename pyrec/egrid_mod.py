"""Python version of egrid_mod.f (YREC 7be4133): the evenly spaced grid of
the rotation solvers, formerly COMMON/EGRID/.

Translates: egrid_mod.f#EGRID_MOD (7be4133)

Routines that USE EGRID_MOD import this module and write
``egrid_mod.<name>``. All set by getgrid.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–32.
chi = np.zeros(_JSON + 1)       # grid variable at the model's shells
echi = np.zeros(_JSON + 1)      # grid variable at the even points
es1 = np.zeros(_JSON + 1)       # mass (g) at the even points
dchi = 0.0                      # spacing of the even points in chi
ntot = 0                        # number of even points (at least 3)
