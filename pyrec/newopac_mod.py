"""Python version of newopac_mod.f (YREC 7be4133): two switches of the
opacity setup, formerly parts of COMMON/NEWOPAC/ and /NWLAOL/.

Translates: newopac_mod.f#NEWOPAC_MOD (7be4133)

Routines that USE NEWOPAC_MOD import this module and write
``newopac_mod.<name>``.
"""

# The two switches (Fortran lines 34–35).
l2z = False     # two opacity tables in Z needed: lzramp or ldifz (setupopac)
llaol = False   # obsolete LAOL opacities; set nowhere, unused? (initializer)
