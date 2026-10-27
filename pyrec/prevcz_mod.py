"""Python version of prevcz_mod.f (YREC 09c3b9b): the convection zone and light-element burning rates carried from one step to the next, formerly COMMON/PREVCZ/.

Translates: prevcz_mod.f#PREVCZ_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PREVCZ_MOD import this module and write ``prevcz_mod.<name>``.
Set and read by liburn and liburn2; main resets them when lexcom is set.
"""

rbase0 = 0.0      # radius of the convection zone's base at the step's start
bli6old = 0.0     # CZ-averaged 6Li, 7Li, 9Be burning rates at the step's end
bli7old = 0.0
bbe9old = 0.0
jenvos = 0        # shell of the envelope's base at the step's end
