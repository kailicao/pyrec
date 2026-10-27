"""Python version of dpmix_mod.f (YREC 09c3b9b): the overshoot region and
the current shell, read by tpgrad, formerly part of COMMON/DPMIX/.

Translates: dpmix_mod.f#DPMIX_MOD (09c3b9b)

Routines that USE DPMIX_MOD import this module and write
``dpmix_mod.<name>``.
"""

iov1 = 0        # first shell of the overshoot region, -1 if none (oversh)
iov2 = 0        # last shell of the overshoot region, -1 if none (oversh)
iovim = 0       # shell whose gradient tpgrad computes (set by its callers)
