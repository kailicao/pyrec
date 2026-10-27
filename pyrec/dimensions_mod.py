"""Python version of dimensions_mod.f (YREC 09c3b9b): the maximum number of
shells, formerly PARAMETER(JSON=5000) in 108 files.

Translates: dimensions_mod.f#DIMENSIONS_MOD (09c3b9b)

Routines import this module and write ``dimensions_mod.json``; arrays of
shells get json + 1 elements, so that index 0 is unused.
"""

json = 5000     # maximum number of shells (PARAMETER)
