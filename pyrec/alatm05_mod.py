"""Python version of alatm05_mod.f (YREC 7be4133): the Allard table's first and last log Teff, formerly COMMON/ALATM05/.

Translates: alatm05_mod.f#ALATM05_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE ALATM05_MOD import this module and write ``alatm05_mod.<name>``.
"""


# Fortran lines 24–25.
al_tefflmin = 0.0                  # set by altabinit, read nowhere
al_tefflmax = 0.0                  # set by altabinit, read by surfbc
