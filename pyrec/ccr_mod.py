"""Python version of ccr_mod.f (YREC 7be4133): the Fermi-Dirac electron table, formerly COMMON/CCR/.

Translates: ccr_mod.f#CCR_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE CCR_MOD import this module and write ``ccr_mod.<name>``.
"""

import numpy as np

# All set by setups (read from unit IFERMI; iden built there); read by eqrelv
# (Fortran lines 24–27).
dval = np.zeros(44)                # density abscissae
etat8 = np.zeros(44)               # degeneracies
dtable = np.zeros((6, 44, 21))     # dtable[k, i, j] for DTABLE(K, I, J)
iden = np.zeros(262, dtype=int)    # density bin (width 0.05) -> table row
