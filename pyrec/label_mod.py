"""Python version of label_mod.f (YREC 09c3b9b): the envelope X and Z of
the current run, for the output, formerly COMMON/LABEL/.

Translates: label_mod.f#LABEL_MOD (09c3b9b)

Routines that USE LABEL_MOD import this module and write
``label_mod.<name>``.
"""

xenv0 = 0.7     # envelope X of the run (main, rscale); initializer 0.7D0
zenv0 = 0.02    # envelope Z of the run (main, rscale); initializer 0.02D0
