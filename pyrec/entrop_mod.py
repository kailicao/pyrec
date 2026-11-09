"""Python version of entrop_mod.f (YREC 7be4133): the rates of change of
the structure in the last model, formerly COMMON/ENTROP/.

Translates: entrop_mod.f#ENTROP_MOD (7be4133)

Routines that USE ENTROP_MOD import this module and write
``entrop_mod.<name>``. Set by main, remapped by hpoint.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# The rates of change of the structure (Fortran lines 28–31).
htt0 = np.zeros(_JSON + 1)      # d log T / dt in the last step
hpp0 = np.zeros(_JSON + 1)      # d log P / dt
hll0 = np.zeros(_JSON + 1)      # relative dL / dt
hrr0 = np.zeros(_JSON + 1)      # d log r / dt
