"""Python version of prevcz_mod.f (YREC 7be4133): the convection zone and
light-element burning rates carried from one step to the next, formerly
COMMON/PREVCZ/.

Translates: prevcz_mod.f#PREVCZ_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PREVCZ_MOD import this module and write ``prevcz_mod.<name>``.
Set and read by liburn and liburn2; main clears lbliold when lexcom is set.
"""

# Fortran lines 24–33.
rbase0 = 0.0      # radius of the convection zone's base at the step's start
bli6old = 0.0     # CZ-averaged 6Li, 7Li, 9Be burning rates at the step's end
bli7old = 0.0
bbe9old = 0.0
lbliold = False   # True when the three hold liburn's end-of-step rates
jenvos = 0        # shell of the envelope's base at the step's end
