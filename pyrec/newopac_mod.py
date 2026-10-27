"""Python version of newopac_mod.f (YREC 09c3b9b): two switches of the
opacity setup, formerly parts of COMMON/NEWOPAC/ and /NWLAOL/.

Translates: newopac_mod.f#NEWOPAC_MOD (09c3b9b)

Routines that USE NEWOPAC_MOD import this module and write
``newopac_mod.<name>``.
"""

l2z = False     # two opacity tables in Z needed: lzramp or ldifz (setupopac)
llaol = False   # obsolete LAOL opacities; set nowhere (initializer)
