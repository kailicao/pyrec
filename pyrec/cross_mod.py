"""Python version of cross_mod.f (YREC 7be4133): the cross-section
adjustment factors of the nuclear rates, formerly COMMON/CROSS/.

Translates: cross_mod.f#CROSS_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case (REAL*8).
Routines that USE CROSS_MOD import this module and write
``cross_mod.<name>``.
"""

import numpy as np

# Fortran lines 24–27.
qs0e = np.zeros(9)      # S'/S factors of the 8 pp/CNO rates (remap), index 0 unused
qqs0ee = np.zeros(9)    # S''/S factors (remap), index 0 unused
fo16 = 0.0              # factor on N15+p -> O16+gamma (remap)
fc12 = 0.0              # factor on N15+p -> C12+alpha (remap)
