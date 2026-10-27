"""Python version of atmprt_mod.f (YREC 09c3b9b): the state at the last point of the atmosphere, formerly COMMON/ATMPRT/.

Translates: atmprt_mod.f#ATMPRT_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ATMPRT_MOD import this module and write ``atmprt_mod.<name>``.
"""

import numpy as np

taul = 0.0             # log tau (qatm; formerly TAU)
ap = 0.0               # log P (qatm, surfp, kcsurfp)
at = 0.0               # log T (qatm, surfp, kcsurfp, alsurfp)
ad = 0.0               # log rho (qatm, when the model is stored)
ao = 0.0               # opacity (qatm, likewise)
afxion = np.zeros(4)   # ionization fractions, index 0 unused (qatm)
