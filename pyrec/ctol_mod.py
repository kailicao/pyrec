"""Python version of ctol_mod.f (YREC 7be4133): the fraction of the Henyey
correction applied, formerly part of COMMON/CTOL/.

Translates: ctol_mod.f#CTOL_MOD (7be4133)

Routines that USE CTOL_MOD import this module and write ``ctol_mod.fcorr``.
"""

# Fortran line 24.
fcorr = 0.0     # set by main each model, raised by crrect each iteration
