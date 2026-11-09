"""Python version of ccout_mod.f (YREC 7be4133): the switch for the stitched
output's header, formerly COMMON/CCOUT/.

Translates: ccout_mod.f#CCOUT_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE CCOUT_MOD import this module and write ``ccout_mod.<name>``.
"""

# The switch (Fortran line 24).
lphhd = False          # set by parmin (to lstch), turned off by putstore
