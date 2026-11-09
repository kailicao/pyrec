"""Python version of zlaol_mod.f (YREC 7be4133): the LAOL89 opacity table for pure metals,
as read, formerly COMMON/ZLAOL/.

Translates: zlaol_mod.f#ZLAOL_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ZLAOL_MOD import this module and write ``zlaol_mod.<name>``.
Read from the file by rdzlaol. Index 0 unused in every dimension.
"""

import numpy as np

# The table and its grids (Fortran lines 24–28).
zolaol = np.zeros((105, 53))   # opacity at (density, temperature); 0 = no entry
zot = np.zeros(53)             # temperature grid (zsulaol takes its log10 in place)
zorho = np.zeros(105)          # density grid
numrhoz = 0                    # number of densities
numtz = 0                      # number of temperatures
