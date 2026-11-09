"""Python version of qenvmix_mod.f (YREC 7be4133): the mixing-length parameter relayed from envint to qenv.

Translates: qenvmix_mod.f#QENVMIX_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE QENVMIX_MOD import this module and write ``qenvmix_mod.<name>``.
"""

# Fortran line 24.
cmixlq = 0.0      # mixing-length parameter (envint sets it to cmixl; qenv reads it)
