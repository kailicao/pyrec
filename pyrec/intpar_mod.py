"""Python version of intpar_mod.f (YREC 7be4133): the settings of the
Bulirsch-Stoer integration, formerly COMMON/INTPAR/.

Translates: intpar_mod.f#INTPAR_MOD (7be4133)

Routines that USE INTPAR_MOD import this module and write
``intpar_mod.<name>``. All three are PHYSICS namelist variables.
"""

# Fortran lines 36–38.
stolr0 = 1.0e-3     # integration tolerance (envint)
imax = 11           # largest number of step subdivisions (bsstep)
nuse = 7            # points used in the extrapolation (bsstep, ratext)
