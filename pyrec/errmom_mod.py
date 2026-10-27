"""Python version of errmom_mod.f (YREC 09c3b9b): the tolerance of the
moment-of-inertia iteration, formerly COMMON/ERRMOM/.

Translates: errmom_mod.f#ERRMOM_MOD (09c3b9b)

Routines that USE ERRMOM_MOD import this module and write
``errmom_mod.toleri``.
"""

# Question (see the page): no routine gives toleri a value; checkj only
# relaxes it temporarily. Zero here, as the Fortran's static storage.
toleri = 0.0
