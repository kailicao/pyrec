"""Python version of envprt_mod.f (YREC 09c3b9b): the state at the current
point of the envelope integration, formerly COMMON/ENVPRT/.

Translates: envprt_mod.f#ENVPRT_MOD (09c3b9b)

Routines that USE ENVPRT_MOD import this module and write
``envprt_mod.<name>``. Set by qenv, read by envint.
"""

import numpy as np

ep = 0.0                # log P
et = 0.0                # log T
er = 0.0                # log r
es = 0.0                # log(m/M)
ed = 0.0                # log rho
eo = 0.0                # opacity
ebeta = 0.0             # beta
edel = np.zeros(4)      # radiative, adiabatic, actual gradient; index 0 unused
efxion = np.zeros(4)    # ionization fractions; index 0 unused
evel = 0.0              # convective velocity
