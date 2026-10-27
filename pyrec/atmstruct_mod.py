"""Python version of atmstruct_mod.f (YREC 09c3b9b): the atmosphere's structure, point by point, for the stitched output, formerly COMMON/ATMSTRUCT/.

Translates: atmstruct_mod.f#ATMSTRUCT_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ATMSTRUCT_MOD import this module and write ``atmstruct_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# All set by envint, read by stitch; index 0 unused.
atmop = np.zeros(_JSON + 1)        # log P
atmot = np.zeros(_JSON + 1)        # log T
atmod = np.zeros(_JSON + 1)        # log rho
atmor = np.zeros(_JSON + 1)        # log r
adels = np.zeros((4, _JSON + 1))   # the three gradients
abetas = np.zeros(_JSON + 1)       # beta
agam1 = np.zeros(_JSON + 1)        # Gamma_1
aqdt = np.zeros(_JSON + 1)
afxions = np.zeros((4, _JSON + 1))  # ionization fractions
atmoo = np.zeros(_JSON + 1)        # opacity
atmocp = np.zeros(_JSON + 1)       # c_P
numatm = 0                         # the number of points
