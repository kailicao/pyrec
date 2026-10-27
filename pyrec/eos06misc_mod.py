"""Python version of eos06misc_mod.f (YREC 09c3b9b): results and workspace
of the OPAL 2006 EOS interpolation, formerly COMMON/EEOS06/, /AAEOS06/,
/BBEOS06/ and /BEOS06/.

Translates: eos06misc_mod.f#EOS06MISC_MOD (09c3b9b)

Routines that USE EOS06MISC_MOD import this module and write
``eos06misc_mod.<name>``. Arrays have index 0 unused.
"""

import numpy as np

# Table dimensions (PRIVATE in the Fortran).
_MX = 5
_MV = 10
_NR = 169
_NT = 197

esact = 0.0                         # one quantity at (T6, rho) (t6rinteos06)
eos = np.zeros(_MV + 1)             # the ten interpolated quantities (esac06)
q = np.zeros(5)                     # quadratics in T (t6rinteos06)
h = np.zeros(5)                     # quadratics in rho (t6rinteos06)
xxh = 0.0                           # set and read nowhere
l1 = 0                              # density indices around the point (esac06)
l2 = 0
l3 = 0
l4 = 0
k1 = 0                              # temperature indices (esac06)
k2 = 0
k3 = 0
k4 = 0
ip = 0                              # points used in each direction (esac06)
iq = 0
zz = np.zeros(_MX + 1)              # Z of each table (readcoeos06)
iri = np.zeros(11, dtype=int)       # inverse of index (seteos06)
index = np.arange(11)               # storage order of the quantities, 1..10
# Number of temperatures for each density, and of densities for each
# temperature (initializers, expanded from the DATA repeat counts).
nta = np.array([0, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197,
     197, 197, 197, 197, 197, 191, 191, 191, 191, 191, 191, 191, 190, 189,
     189, 185, 179, 170, 149, 149, 133, 125, 123, 122, 120, 115, 113, 107,
     102, 80, 80, 72, 68, 66, 64, 62, 56, 54, 52, 51, 50, 50, 49, 47, 45, 45,
     43, 42, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
     40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 39, 37, 36, 35, 34, 32,
     31, 30, 29, 27, 26])
nra = np.array([0, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169,
     169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169,
     168, 167, 167, 166, 165, 164, 163, 163, 162, 161, 160, 159, 159, 158,
     130, 130, 129, 128, 128, 126, 126, 125, 125, 124, 122, 121, 120, 120,
     119, 119, 118, 118, 118, 118, 118, 118, 117, 117, 116, 116, 115, 115,
     114, 114, 114, 114, 113, 113, 113, 113, 113, 113, 113, 113, 111, 111,
     111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111,
     111, 111, 111, 111, 111, 111, 110, 110, 110, 110, 110, 109, 109, 109,
     109, 109, 109, 108, 108, 107, 107, 107, 107, 107, 106, 106, 105, 104,
     104, 103, 103, 103, 103, 103, 103, 103, 103, 102, 102, 102, 102, 102,
     102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 100, 100, 100,
     100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
     100, 100, 100, 100, 99, 99, 99, 99, 99, 99, 99, 99, 99, 98, 98, 98, 98,
     98, 98, 97, 97, 97, 97, 95, 94, 87, 87, 87, 87, 87, 87])

assert len(nta) == _NR + 1 and len(nra) == _NT + 1
