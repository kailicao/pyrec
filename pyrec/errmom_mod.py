"""Python version of errmom_mod.f (YREC 7be4133): the tolerance of the
moment-of-inertia iteration, formerly COMMON/ERRMOM/.

Translates: errmom_mod.f#ERRMOM_MOD (7be4133)

Routines that USE ERRMOM_MOD import this module and write
``errmom_mod.toleri``.
"""

# Given a value in the module; checkj only relaxes it temporarily (Fortran
# line 28).
toleri = 1.0e-12
