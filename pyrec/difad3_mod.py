"""Python version of difad3_mod.f (YREC 09c3b9b): the Zahn-type factors of
the diffusion-plus-advection treatment of angular momentum, formerly
COMMON/DIFAD3/.

Translates: difad3_mod.f#DIFAD3_MOD (09c3b9b)

Routines that USE DIFAD3_MOD import this module and write
``difad3_mod.<name>``. All arrays have JSON + 1 elements, index 0 unused.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

facd2 = np.zeros(_JSON + 1)     # setupv; read by vcirc
facd3 = np.zeros(_JSON + 1)     # setupv; read by vcirc
vesd2 = np.zeros(_JSON + 1)     # vcirc; read nowhere else
vesd3 = np.zeros(_JSON + 1)     # vcirc; read nowhere else
# Question (see the page): ecod5, ecod6, fgeom, ev0-ev2b, ess and egsf
# are read by dadcoeft but set by no routine (zero here).
ecod5 = np.zeros(_JSON + 1)
ecod6 = np.zeros(_JSON + 1)
fgeom = np.zeros(_JSON + 1)
fv0 = np.zeros(_JSON + 1)       # setupv; rescaled by vcirc
fv1a = np.zeros(_JSON + 1)      # setupv; read nowhere
fv1b = np.zeros(_JSON + 1)      # setupv; read nowhere
fv2a = np.zeros(_JSON + 1)      # setupv; read nowhere
fv2b = np.zeros(_JSON + 1)      # setupv; read nowhere
ev0 = np.zeros(_JSON + 1)       # read by dadcoeft, set nowhere
ev1a = np.zeros(_JSON + 1)
ev1b = np.zeros(_JSON + 1)
ev2a = np.zeros(_JSON + 1)
ev2b = np.zeros(_JSON + 1)
dss = np.zeros(_JSON + 1)       # r |v_SS| (vcirc); read nowhere
dgsf = np.zeros(_JSON + 1)      # r |v_GSF| (vcirc); read nowhere
ess = np.zeros(_JSON + 1)       # read (and rescaled) by dadcoeft, set nowhere
egsf = np.zeros(_JSON + 1)      # read (and rescaled) by dadcoeft, set nowhere
