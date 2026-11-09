"""Python version of fluxes_mod.f (YREC 7be4133): the solar neutrino
fluxes and capture rates, formerly COMMON/FLUXES/.

Translates: fluxes_mod.f#FLUXES_MOD (7be4133)

Routines that USE FLUXES_MOD import this module and write
``fluxes_mod.<name>``.
"""

import numpy as np

# Fortran lines 24–27.
flux = np.zeros(11)         # fluxes per gram at one shell (engeb), index 0 unused
fluxtot = np.zeros(11)      # mass-weighted sums over the star (coefft, main)
clsnu = 0.0                 # Cl37 capture rate, SNU (wrtout)
gasnu = 0.0                 # Ga71 capture rate, SNU (wrtout)
