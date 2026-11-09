"""Python version of alatm01_mod.f (YREC 7be4133): the Allard model-atmosphere table, formerly COMMON/ALATM01/.

Translates: alatm01_mod.f#ALATM01_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ALATM01_MOD import this module and write ``alatm01_mod.<name>``.
"""

import numpy as np

# Sizes, PRIVATE in the Fortran (Fortran lines 28–29).
_NGA = 25
_NTA = 250

# All set by alfilein (Fortran lines 31–42).
teffls = np.zeros(_NTA + 1)        # log Teff grid
gls = np.zeros(_NGA + 1)           # log g grid
fehs = np.zeros(_NGA + 1)          # [Fe/H] values in the file
alphas = np.zeros(_NGA + 1)        # [alpha/Fe] values in the file
pls = np.zeros((_NTA + 1, _NGA + 1))      # log P where T = Teff
p100ls = np.zeros((_NTA + 1, _NGA + 1))   # log P at tau = 100
t100ls = np.zeros((_NTA + 1, _NGA + 1))   # log T at tau = 100
loldnextgen = False                # the file's format
nteff = 0
ngl = 0
nfeh = 0
nalpha = 0
