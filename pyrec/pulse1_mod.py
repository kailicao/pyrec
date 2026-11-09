"""Python version of pulse1_mod.f (YREC 7be4133): derivatives at every shell
for the pulsation output, formerly COMMON/PULSE1/.

Translates: pulse1_mod.f#PULSE1_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PULSE1_MOD import this module and write ``pulse1_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# At each shell, index 0 unused (coefft, when lpulse; pqdp and pqdt always)
# (Fortran lines 28–36).
pqdp = np.zeros(_JSON + 1)     # dln rho / dln P
pqed = np.zeros(_JSON + 1)     # energy generation's derivative in rho
pqet = np.zeros(_JSON + 1)     # ... in T
pqod = np.zeros(_JSON + 1)     # opacity's derivative in rho
pqot = np.zeros(_JSON + 1)     # ... in T
pqcp = np.zeros(_JSON + 1)     # c_P
prmu = np.zeros(_JSON + 1)
pqdt = np.zeros(_JSON + 1)     # dln rho / dln T
pemu = np.zeros(_JSON + 1)

# The envelope switch (Fortran line 37).
lpumod = False                 # envelope pulsation output on (envint)
