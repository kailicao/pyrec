"""Python version of sound_mod.f (YREC 7be4133): Gamma_1 at each shell and the
sound-speed output switch, formerly COMMON/SOUND/.

Translates: sound_mod.f#SOUND_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE SOUND_MOD import this module and write ``sound_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Gamma_1 and the switch (Fortran lines 28–29).
gam1 = np.zeros(_JSON + 1)   # Gamma_1 at each shell (coefft), index 0 unused
lsound = False               # sound-speed output at the run's last model (main)
