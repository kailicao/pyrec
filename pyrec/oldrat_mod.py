"""Python version of oldrat_mod.f (YREC 09c3b9b): the light-element burning
rates at the start of the step, formerly COMMON/OLDRAT/.

Translates: oldrat_mod.f#OLDRAT_MOD (09c3b9b)

Routines that USE OLDRAT_MOD import this module and write
``oldrat_mod.<name>``. Set by lirate88 (j = 2), getw and bursmix.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

rli60 = np.zeros(_JSON + 1)     # Li6 + p rate
rli70 = np.zeros(_JSON + 1)     # Li7 + p rate
rbe90 = np.zeros(_JSON + 1)     # Be9 + p rate
