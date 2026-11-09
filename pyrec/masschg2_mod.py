"""Python version of masschg2_mod.f (YREC 7be4133): the thermal response
of the convective envelope to accretion, formerly COMMON/MASSCHG2/.

Translates: masschg2_mod.f#MASSCHG2_MOD (7be4133)

Routines that USE MASSCHG2_MOD import this module and write
``masschg2_mod.<name>``.
"""

# Fortran lines 24–28.
sacc = 0.0      # entropy-like quantity at the envelope's base (massloss); read nowhere
scen = 0.0      # envelope's mean of it (massloss, which then sets it to 0)
smass0 = 0.0    # total mass after the step (mdot); read nowhere
dlogp = 0.0     # change of log P from accretion (mdot; read by coefft)
dlogt = 0.0     # change of log T from accretion (mdot; read by coefft)
