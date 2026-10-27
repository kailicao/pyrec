"""Python version of vfact_mod.f (YREC 09c3b9b): factors of the circulation velocities and stability criteria at each interface, formerly COMMON/VFACT/.

Translates: vfact_mod.f#VFACT_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE VFACT_MOD import this module and write ``vfact_mod.<name>``.
Set by setupv, read by vcirc. Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

fact1 = np.zeros(_JSON + 1)   # Eddington-Sweet velocity: ves = fact1*fact2*omega**2
fact2 = np.zeros(_JSON + 1)
fact3 = np.zeros(_JSON + 1)   # |H_T r| for the GSF velocity (read nowhere)
fact4 = np.zeros(_JSON + 1)   # mu currents
fact5 = np.zeros(_JSON + 1)   # Richardson criterion with mu gradients
fact7 = np.zeros(_JSON + 1)   # d omega/dr terms (ldifad only)
fact8 = np.zeros(_JSON + 1)   # stays 0
