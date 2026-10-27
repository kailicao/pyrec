"""Python version of newrat_mod.f (YREC 09c3b9b): the light-element burning
rates of the current model, formerly COMMON/NEWRAT/.

Translates: newrat_mod.f#NEWRAT_MOD (09c3b9b)

Routines that USE NEWRAT_MOD import this module and write
``newrat_mod.<name>``. Set by lirate88 (j = 1).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

rli6 = np.zeros(_JSON + 1)      # Li6 + p rate
rli7 = np.zeros(_JSON + 1)      # Li7 + p rate
rbe9 = np.zeros(_JSON + 1)      # Be9 + p rate
