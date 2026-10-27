"""Python version of alatm02_mod.f (YREC 09c3b9b): the valid range of the Allard table, formerly COMMON/ALATM02/.

Translates: alatm02_mod.f#ALATM02_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ALATM02_MOD import this module and write ``alatm02_mod.<name>``.
"""

import numpy as np

_NTA = 250                         # PRIVATE in the Fortran

# Set by altabinit (tefflmin, tefflmax also by alprint).
glmin = np.zeros(_NTA + 1)         # per log Teff row: allowed log g range
glmax = np.zeros(_NTA + 1)
iglmin = np.zeros(_NTA + 1, dtype=int)   # first and last valid column
iglmax = np.zeros(_NTA + 1, dtype=int)
tefflmin = 0.0                     # log Teff range
tefflmax = 0.0
glxmin = 0.0                       # extremes of log g over all rows
glxmax = 0.0
