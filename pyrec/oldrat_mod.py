"""Python version of oldrat_mod.f (YREC 7be4133): the light-element burning
rates at the start of the step, formerly COMMON/OLDRAT/.

Translates: oldrat_mod.f#OLDRAT_MOD (7be4133)

Routines that USE OLDRAT_MOD import this module and write
``oldrat_mod.<name>``. Set by lirate88 (j = 2), getw and bursmix.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–30.
rli60 = np.zeros(_JSON + 1)     # Li6 + p rate
rli70 = np.zeros(_JSON + 1)     # Li7 + p rate
rbe90 = np.zeros(_JSON + 1)     # Be9 + p rate
