"""Python version of temp2_mod.f (YREC 09c3b9b): the velocities of rotational mixing at each interface, formerly COMMON/TEMP2/.

Translates: temp2_mod.f#TEMP2_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE TEMP2_MOD import this module and write ``temp2_mod.<name>``.
Set by vcirc (hle by codiff). Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

ves = np.zeros(_JSON + 1)     # Eddington-Sweet circulation velocity
ves0 = np.zeros(_JSON + 1)    # its value at the previous iteration
vss = np.zeros(_JSON + 1)     # secular-shear velocity
vss0 = np.zeros(_JSON + 1)
hle = np.zeros(_JSON + 1)     # length scale, always 0 (codiff)
vgsf = np.zeros(_JSON + 1)    # Goldreich-Schubert-Fricke velocity
vgsf0 = np.zeros(_JSON + 1)
vmu = np.zeros(_JSON + 1)     # mu-current velocity
