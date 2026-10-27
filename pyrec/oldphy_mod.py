"""Python version of oldphy_mod.f (YREC 09c3b9b): the physics of the model
at the start of the step, formerly COMMON/OLDPHY/.

Translates: oldphy_mod.f#OLDPHY_MOD (09c3b9b)

Routines that USE OLDPHY_MOD import this module and write
``oldphy_mod.<name>``. Saved by hpoint, read by midmod (amuo also getw).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

delo = np.zeros(_JSON + 1)      # actual temperature gradient
delao = np.zeros(_JSON + 1)     # adiabatic gradient
amuo = np.zeros(_JSON + 1)      # mean molecular weight
oo = np.zeros(_JSON + 1)        # opacity
cpo = np.zeros(_JSON + 1)       # c_P
qdto = np.zeros(_JSON + 1)      # thermodynamic derivative
velo = np.zeros(_JSON + 1)      # convective velocity
visco = np.zeros(_JSON + 1)     # viscosity
thdifo = np.zeros(_JSON + 1)    # thermal diffusivity
esumo = np.zeros(_JSON + 1)     # energy generation
delro = np.zeros(_JSON + 1)     # radiative gradient
epso = np.zeros(_JSON + 1)      # energy generation plus seg(6), seg(7)
