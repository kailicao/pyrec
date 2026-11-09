"""Python version of confac_mod.f (YREC 7be4133): the unit conversions of
the diffusion routines, formerly COMMON/CONFAC/.

Translates: confac_mod.f#CONFAC_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE CONFAC_MOD import this module and write ``confac_mod.<name>``.
"""

# Set by microdiff_setup and setup_grsett (Fortran lines 24–27).
con_rad = 0.0          # 1/Rsun
con_mass = 0.0         # 1e-2 * con_rad**3
con_temp = 0.0         # 1e-7
con_time = 0.0         # 2.7e13 yr in seconds
