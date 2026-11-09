"""Python version of label_mod.f (YREC 7be4133): the envelope X and Z of
the current run, for the output, formerly COMMON/LABEL/.

Translates: label_mod.f#LABEL_MOD (7be4133)

Routines that USE LABEL_MOD import this module and write
``label_mod.<name>``.
"""

# The envelope X and Z (Fortran lines 34–35).
xenv0 = 0.7     # envelope X of the run (main, rscale); initializer 0.7D0
zenv0 = 0.02    # envelope Z of the run (main, rscale); initializer 0.02D0
