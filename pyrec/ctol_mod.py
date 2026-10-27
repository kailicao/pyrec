"""Python version of ctol_mod.f (YREC 09c3b9b): the fraction of the Henyey
correction applied, formerly part of COMMON/CTOL/.

Translates: ctol_mod.f#CTOL_MOD (09c3b9b)

Routines that USE CTOL_MOD import this module and write ``ctol_mod.fcorr``.
"""

fcorr = 0.0     # set by main each model, raised by crrect each iteration
