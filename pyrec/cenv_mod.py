"""Python version of cenv_mod.f (YREC 7be4133): the option to move the envelope's fitting point, formerly COMMON/CENV/.

Translates: cenv_mod.f#CENV_MOD (7be4133)

Module-level names, as in the Fortran module, in lower case. Routines that
USE CENV_MOD import this module and write ``cenv_mod.<name>``.
"""

# Fortran lines 24–25.
senv0 = 0.0            # log10(m_fit/M); set by main, adjusted by starin
lsenv0 = False         # set by main, read by starin
