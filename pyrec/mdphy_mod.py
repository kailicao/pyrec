"""Python version of mdphy_mod.f (YREC 09c3b9b): the physics at each shell
at an intermediate time of the model step, formerly COMMON/MDPHY/.

Translates: mdphy_mod.f#MDPHY_MOD (09c3b9b)

Routines that USE MDPHY_MOD import this module and write
``mdphy_mod.<name>``. Set by midmod (amum also by checkc and getw).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

amum = np.zeros(_JSON + 1)      # mean molecular weight
cpm = np.zeros(_JSON + 1)       # c_P
delm = np.zeros(_JSON + 1)      # actual temperature gradient
delam = np.zeros(_JSON + 1)     # adiabatic gradient
delrm = np.zeros(_JSON + 1)     # radiative gradient
esumm = np.zeros(_JSON + 1)     # energy generation
om = np.zeros(_JSON + 1)        # opacity
qdtm = np.zeros(_JSON + 1)      # thermodynamic derivative
thdifm = np.zeros(_JSON + 1)    # thermal diffusivity
velm = np.zeros(_JSON + 1)      # convective velocity; read nowhere
viscm = np.zeros(_JSON + 1)     # viscosity
epsm = np.zeros(_JSON + 1)      # energy generation plus seg(6), seg(7)
