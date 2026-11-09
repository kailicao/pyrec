"""Python version of disk_mod.f (YREC 7be4133): the age used for disk
locking, formerly part of COMMON/DISK/.

Translates: disk_mod.f#DISK_MOD (7be4133)

Routines that USE DISK_MOD import this module and write ``disk_mod.sage``.
"""

# Fortran line 24.
sage = 0.0      # total age in Gyr (main, after each model)
