"""Python version of roten_mod.f (YREC 7be4133): the change of rotational
kinetic energy at each shell, formerly COMMON/ROTEN/.

Translates: roten_mod.f#ROTEN_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ROTEN_MOD import this module and write ``roten_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Fortran line 28.
derot = np.zeros(_JSON + 1)   # dt (E_rot - E_rot,old)/dm (coefft), index 0 unused
