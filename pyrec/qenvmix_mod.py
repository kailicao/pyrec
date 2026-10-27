"""Python version of qenvmix_mod.f (YREC 09c3b9b): the mixing-length parameter relayed from envint to qenv.

Translates: qenvmix_mod.f#QENVMIX_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE QENVMIX_MOD import this module and write ``qenvmix_mod.<name>``.
"""

cmixlq = 0.0      # mixing-length parameter (envint sets it to cmixl; qenv reads it)
