"""Python version of liov_mod.f (YREC 7be4133): the overshoot distance below
the surface convection zone, for light-element burning, formerly
COMMON/LIOV/.

Translates: liov_mod.f#LIOV_MOD (7be4133)

Routines that USE LIOV_MOD import this module and write ``liov_mod.<name>``.
"""

# Fortran lines 24–25.
pscal0 = 0.0    # overshoot distance (cm) at the start of the step
pscal1 = 0.0    # overshoot distance (cm) at the end of the step (getw, main)
