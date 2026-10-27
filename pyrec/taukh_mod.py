"""Python version of taukh_mod.f (YREC 09c3b9b): factors of the circulation velocities at each interface, formerly COMMON/TAUKH/.

Translates: taukh_mod.f#TAUKH_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE TAUKH_MOD import this module and write ``taukh_mod.<name>``.
Set by setupv, read by vcirc. Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

fact6 = np.zeros(_JSON + 1)   # Eddington-Sweet term in d(Omega^2 r^2)
fes1 = np.zeros(_JSON + 1)    # Zahn-type terms
fes2 = np.zeros(_JSON + 1)
fgsfj = np.zeros(_JSON + 1)   # GSF factor of an earlier form (read nowhere)
fgsfk = np.zeros(_JSON + 1)   # GSF factor used by vcirc
fes3 = np.zeros(_JSON + 1)
