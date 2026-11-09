"""Python version of pualpha_mod.f (YREC 7be4133): mixing-length quantities
for the pulsation output, formerly COMMON/PUALPHA/.

Translates: pualpha_mod.f#PUALPHA_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PUALPHA_MOD import this module and write ``pualpha_mod.<name>``.
"""

# At the current point (Fortran lines 28–30): set to 0 by tpgrad and
# computed nowhere (see TPGRAD's page); written by envint.
alfmlt = 0.0
phmlt = 0.0
cmxmlt = 0.0
