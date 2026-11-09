"""Python version of dimensions_mod.f (YREC 7be4133): the maximum number of
shells, formerly PARAMETER(JSON=5000) in 108 files.

Translates: dimensions_mod.f#DIMENSIONS_MOD (7be4133)

Routines import this module and write ``dimensions_mod.json``; arrays of
shells get json + 1 elements, so that index 0 is unused.
"""

# Fortran line 42.
json = 5000     # maximum number of shells (PARAMETER)
