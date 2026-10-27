"""Python version of runloop_mod.f (YREC 09c3b9b): per-run state of the run loop that is not in the namelist, formerly COMMON/CKIND/ and COMMON/SETT/.

Translates: runloop_mod.f#RUNLOOP_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE RUNLOOP_MOD import this module and write ``runloop_mod.<name>``.
Each array is indexed by the run (1 to 50; index 0 unused).
"""

import numpy as np

# Rescaling targets: rescal[1..4, nk] = mass, X, Z, core mass (parmin;
# changed by the calibrations setcal, setscal, chkcal, chkscal).
rescal = np.zeros((5, 51))
iresca = np.zeros(51, dtype=np.int32)   # run kind: 1 evolve, 2 rescale, 3 both
lendag = np.zeros(51, dtype=bool)       # evolve until a stop condition
lsetdt = np.zeros(51, dtype=bool)       # fixed time step
