"""Python version of intvar_mod.f (YREC 09c3b9b): the structure at the
interfaces between shells, formerly COMMON/INTVAR/.

Translates: intvar_mod.f#INTVAR_MOD (09c3b9b)

Routines that USE INTVAR_MOD import this module and write
``intvar_mod.<name>``. All set by setupv; index i is the interface between
shells i-1 and i.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

alm = np.zeros(_JSON + 1)       # luminosity (cgs)
delami = np.zeros(_JSON + 1)    # adiabatic gradient
delmi = np.zeros(_JSON + 1)     # radiative gradient
dm = np.zeros(_JSON + 1)        # density
epsilm = np.zeros(_JSON + 1)    # energy generation rate
hgm = np.zeros(_JSON + 1)       # gravity
hs3 = np.zeros(_JSON + 1)       # mass
pm = np.zeros(_JSON + 1)        # pressure
qdtmi = np.zeros(_JSON + 1)     # thermodynamic derivative
rm = np.zeros(_JSON + 1)        # radius (mean of the two shells')
tm = np.zeros(_JSON + 1)        # temperature
