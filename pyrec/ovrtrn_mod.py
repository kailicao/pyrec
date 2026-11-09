"""Python version of ovrtrn_mod.f (YREC 7be4133): the convective
overturn time and photospheric pressure at the start and end of a step,
formerly COMMON/OVRTRN/.

Translates: ovrtrn_mod.f#OVRTRN_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE OVRTRN_MOD import this module and write ``ovrtrn_mod.<name>``.
(LNEWTCZ and LCALCENV, formerly in the block, are in nmlpar_mod.)
"""

# Fortran lines 24–28.
taucz = 0.0       # convective overturn time, s (taucal, tauint, tauintnew)
taucz0 = 0.0      # taucz at the start of the step (starin, wrtout)
pphot = 0.0       # photospheric pressure (envint)
pphot0 = 0.0      # pphot at the start of the step (starin, wrtout)
fracstep = 0.0    # fraction of the step for the wind (starin, getw)
