"""Python version of eeeos06_mod.f (YREC 7be4133): workspace of the OPAL
2006 EOS interpolation, formerly COMMON/EEEOS06/.

Translates: eeeos06_mod.f#EEEOS06_MOD (7be4133)

Routines that USE EEEOS06_MOD import this module and write
``eeeos06_mod.<name>``. Arrays keep Fortran's index order, with index 0
unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran) (Fortran lines 27–29).
_MX = 5
_NR = 169
_NT = 197

# Fortran lines 31–32.
epl = np.zeros((_MX + 1, _NT + 1, _NR + 1))   # one quantity's tables (esac06)
xx = np.zeros(_MX + 1)                        # X of the tables, from xa (seteos06)
