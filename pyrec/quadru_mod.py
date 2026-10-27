"""Python version of quadru_mod.f (YREC 09c3b9b): the quadrupole moment of the potential and the gravity at each shell, formerly COMMON/QUADRU/.

Translates: quadru_mod.f#QUADRU_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE QUADRU_MOD import this module and write ``quadru_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

quad = np.zeros(_JSON + 1)    # quadrupole moment (getqua), index 0 unused
gg = np.zeros(_JSON + 1)      # local gravity used by getqua
