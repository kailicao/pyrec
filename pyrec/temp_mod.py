"""Python version of temp_mod.f (YREC 09c3b9b): physics at each shell for rotational mixing, formerly COMMON/TEMP/.

Translates: temp_mod.f#TEMP_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE TEMP_MOD import this module and write ``temp_mod.<name>``.
Index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

cp = np.zeros(_JSON + 1)      # c_P (physic)
hamu = np.zeros(_JSON + 1)    # mean molecular weight, fully ionized (physic)
sqdt = np.zeros(_JSON + 1)    # dln rho / dln T (physic)
thdif = np.zeros(_JSON + 1)   # radiative diffusivity (viscos)
visc = np.zeros(_JSON + 1)    # viscosity, radiative plus molecular (viscos)
