"""Python version of constants_mod.f (YREC 09c3b9b): the physical and
mathematical constants, formerly COMMON/CONST/, /CONST1/, /CONST2/,
/CONST3/ and C_2 of /CWIND/.

Translates: constants_mod.f#CONSTANTS_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case (all REAL*8).
Routines that USE CONSTANTS_MOD import this module and write
``constants_mod.<name>``. Apart from clsun and crsun, the values are set at
run time by setups (and cln, cc13, cc23, c_2 by parmin).
"""

clsun = 3.8515e33       # solar luminosity (initializer; PHYSICS namelist)
clsunl = 0.0            # log10(clsun) (setups)
clnsun = 0.0            # ln(10)/clsun (setups)
cmsun = 0.0             # solar mass, 1.9891e33 g (setups)
cmsunl = 0.0            # log10(cmsun) (setups; read nowhere)
crsun = 6.9598e10       # solar radius (initializer; PHYSICS namelist)
crsunl = 0.0            # log10(crsun) (setups)
cmbol = 0.0             # solar bolometric magnitude, 4.79 (setups)
cln = 0.0               # ln(10) (parmin, setups)
clni = 0.0              # 1/ln(10) (setups)
c4pi = 0.0              # 4 pi (setups)
c4pil = 0.0             # log10(4 pi)
c4pi3l = 0.0            # log10(4 pi / 3)
cc13 = 0.0              # 1/3 (parmin, setups)
cc23 = 0.0              # 2/3 (parmin, setups)
cpi = 0.0               # pi (setups)
cgas = 0.0              # gas constant (setups)
ca3 = 0.0               # 4 sigma / (3 c) (setups)
ca3l = 0.0              # log10(ca3) (setups; read nowhere)
csig = 0.0              # Stefan-Boltzmann constant (setups)
csigl = 0.0             # log10(csig)
cgl = 0.0               # log10(G) (setups)
cmkh = 0.0              # Saha constant (setups)
cmkhn = 0.0             # set and read nowhere
cdelrl = 0.0            # constant of the radiative gradient (setups)
cmixl2 = 0.0            # 1/3 (setups)
cmixl3 = 0.0            # 16 sqrt(2) sigma (setups)
clndp = 0.0             # ln(10) (setups)
csecyr = 0.0            # seconds per year, 3.1558e7 (setups)
c_2 = 0.0               # centrifugal factor of the wind law, 0.0506 (parmin)
