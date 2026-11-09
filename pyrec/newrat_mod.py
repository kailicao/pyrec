"""Python version of newrat_mod.f (YREC 7be4133): the light-element burning
rates of the current model, formerly COMMON/NEWRAT/.

Translates: newrat_mod.f#NEWRAT_MOD (7be4133)

Routines that USE NEWRAT_MOD import this module and write
``newrat_mod.<name>``. Set by lirate88 (j = 1).
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran lines 28–30.
rli6 = np.zeros(_JSON + 1)      # Li6 + p rate
rli7 = np.zeros(_JSON + 1)      # Li7 + p rate
rbe9 = np.zeros(_JSON + 1)      # Be9 + p rate
