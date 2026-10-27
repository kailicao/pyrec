"""Python version of acdpth_mod.f (YREC 09c3b9b): the acoustic-depth outputs of calcad and their switches, formerly COMMON/ACDPTH/.

Translates: acdpth_mod.f#ACDPTH_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ACDPTH_MOD import this module and write ``acdpth_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

# Results of calcad (written by wrtout).
tauczn = 0.0                       # t_CZ / t_norm (wrtout resets it to 0)
deladj = np.zeros(_JSON + 1)       # nabla_ad at each point (calcad)
tauhe = 0.0
tnorm = 0.0
tcz = 0.0
whe = 0.0
acatmr = np.zeros(_JSON + 1)       # unused
acatmd = np.zeros(_JSON + 1)       # unused
acatmp = np.zeros(_JSON + 1)       # unused
acatmt = np.zeros(_JSON + 1)       # unused
tatmos = 0.0                       # always 0 (calcad)
# Ages (Gyr) of the outputs, index 0 unused (initializer, from parmin).
ageout = np.array([0.0, 0.5, 1.0, 5.0, 10.0, 20.0])
# Switches (set by main from lacout and the age).
lclcd = False
ljlast = False
ljwrt = False
ladon = False
laoly = False                      # unused
# Question (see the page): in no namelist and set nowhere, so always False.
lacout = False
