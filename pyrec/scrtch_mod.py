"""Python version of scrtch_mod.f (YREC 09c3b9b): the interior's local physics at each shell, formerly COMMON/SCRTCH/.

Translates: scrtch_mod.f#SCRTCH_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE SCRTCH_MOD import this module and write ``scrtch_mod.<name>``.
Set by coefft at each shell (so, sdel, svel also by physic). Index 0 unused
in every dimension.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

sesum = np.zeros(_JSON + 1)          # nuclear energy generation, EG(1)..EG(5)
seg = np.zeros((8, _JSON + 1))       # 1-5: fractions of sesum; 6: neutrinos; 7: gravitational
sbeta = np.zeros(_JSON + 1)          # beta, gas over total pressure
seta = np.zeros(_JSON + 1)           # degeneracy parameter eta
locons = np.zeros(_JSON + 1, dtype=bool)   # conduction dominates (read nowhere)
so = np.zeros(_JSON + 1)             # opacity
sdel = np.zeros((4, _JSON + 1))      # 1: nabla_rad, 2: nabla, 3: nabla_ad
sfxion = np.zeros((4, _JSON + 1))    # ionization fractions
svel = np.zeros(_JSON + 1)           # convective velocity, cm/s
scp = np.zeros(_JSON + 1)            # c_P
