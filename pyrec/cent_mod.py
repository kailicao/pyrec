"""Python version of cent_mod.f (YREC 09c3b9b): the central conditions and the surface convection zone's base, formerly COMMON/CENT/.

Translates: cent_mod.f#CENT_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE CENT_MOD import this module and write ``cent_mod.<name>``.
"""

# All set by wrtout; read by main for the neutrino-flux output.
tcenter = 0.0          # central log T
pcenter = 0.0          # central log P
dcenter = 0.0          # central log rho
czsurf_m = 0.0         # mass of the surface convection zone (Msun)
czsurf_r = 0.0         # radius of its base (Rsun)
