"""Python version of scveos_mod.f (YREC 09c3b9b): the tables of the SCV equation of state, formerly COMMON/SCVEOS/.

Translates: scveos_mod.f#SCVEOS_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE SCVEOS_MOD import this module and write ``scveos_mod.<name>``.
(LSCV, formerly in the block, is in nmlpar_mod.) Index 0 unused in every
dimension.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_NPS = 76      # pressures per temperature
_NTS = 63      # temperatures

tlogx = np.zeros(_NTS + 1)                      # log T of the tables (setups)
tablex = np.zeros((_NTS + 1, _NPS + 1, 13))     # hydrogen (setups; column 12 setscv)
tabley = np.zeros((_NTS + 1, _NPS + 1, 13))     # helium
smix = np.zeros((_NTS + 1, _NPS + 1))           # entropy of mixing (setscv)
tablez = np.zeros((_NTS + 1, _NPS + 1, 14))     # metals (setups)
tablenv = np.zeros((_NTS + 1, _NPS + 1, 13))    # envelope mixture (setscv)
nptsx = np.zeros(_NTS + 1, dtype=np.int32)      # pressure points at each T
idtscv = 0     # stencil start in T, kept between calls (formerly IDT)
idp = 0        # stencil start in P
