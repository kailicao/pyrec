"""Python version of comp_mod.f (YREC 09c3b9b): the envelope composition
and the star's total and envelope masses, formerly COMMON/COMP/.

Translates: comp_mod.f#COMP_MOD (09c3b9b)

Module-level names, as in the Fortran module (all REAL*8). Routines that
USE COMP_MOD import this module and write ``comp_mod.xenv``.
"""

import numpy as np

xenv = 0.0              # envelope X (set in starin)
zenv = 0.0              # envelope Z (starin; rscale, parmin may reset it)
zenvm = 0.0             # envelope Z other than C, N, O (starin; read by mixcz)
amuenv = 0.0            # sum of X_i/A_i over the 12 species (starin)
fxenv = np.zeros(13)    # number fractions of the 12 species, index 0 unused
xnew = 0.0              # surface X of the current surface opacity table
znew = 0.0              # surface Z of the current surface opacity table
stotal = 0.0            # log10 of the total mass
senv = 0.0              # log10(m/M) at the envelope fitting point
