"""Python version of oldrot2_mod.f (YREC 09c3b9b): the d(theta)/dt terms of
the Maeder-Zahn treatment of angular momentum, formerly COMMON/OLDROT2/.

Translates: oldrot2_mod.f#OLDROT2_MOD (09c3b9b)

Routines that USE OLDROT2_MOD import this module and write
``oldrot2_mod.<name>``. Index i is the interface between shells i-1 and i.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

tho = np.zeros(_JSON + 1)       # theta of the old model (hpoint)
thn = np.zeros(_JSON + 1)       # d theta / dt term (vcirc)
thm = np.zeros(_JSON + 1)       # r^4 / (sqrt(m) G) (setupv)
ddelm = np.zeros(_JSON + 1)     # nabla_ad - nabla (hpoint, then setupv)
festime = np.zeros(_JSON + 1)   # P / (g dnabla rho) (setupv, rescaled by vcirc)
thp = np.zeros(_JSON + 1)       # theta at the start of the substep (getw)
qwrst = np.zeros(_JSON + 1)     # dOmega/dr of the old model (hpoint)
wmst = np.zeros(_JSON + 1)      # Omega at the start of the substep (getw)
qwrmst = np.zeros(_JSON + 1)    # dOmega/dr at the start (getw); read nowhere
