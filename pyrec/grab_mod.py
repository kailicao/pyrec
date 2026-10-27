"""Python version of grab_mod.f (YREC 09c3b9b): the He3-burning energy
generation and luminosities, formerly COMMON/GRAB/.

Translates: grab_mod.f#GRAB_MOD (09c3b9b)

Routines that USE GRAB_MOD import this module and write ``grab_mod.<name>``.
"""

import numpy as np

from pyrec import dimensions_mod

_JSON = dimensions_mod.json

enghe3 = 0.0                    # He3+He3 energy generation (engeb)
he3all = 0.0                    # He3+He3 plus He3+He4 (engeb)
umhe33 = np.zeros(_JSON + 1)    # enghe3 * shell mass / Lsun (coefft); read nowhere
umheal = np.zeros(_JSON + 1)    # he3all * shell mass / Lsun (coefft); read nowhere
