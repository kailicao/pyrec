"""Python version of pulse2_mod.f (YREC 09c3b9b): the physics at the current point of the envelope integration, for the pulsation output, formerly COMMON/PULSE2/.

Translates: pulse2_mod.f#PULSE2_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case. Routines that
USE PULSE2_MOD import this module and write ``pulse2_mod.<name>``.
Set by qenv and qatm, read by envint. coefft also writes qp, qt and qtl,
which are its own working variables in the Fortran (see the page).
"""

qqdp = 0.0      # dln rho / dln P
qqed = 0.0      # energy generation's derivatives (set to 0)
qqet = 0.0
qqod = 0.0      # opacity's derivatives in rho and T
qqot = 0.0
qdel = 0.0      # nabla
qdela = 0.0     # nabla_ad
qqcp = 0.0      # c_P
qrmu = 0.0
qtl = 0.0       # log T
qpl = 0.0       # log P
qdl = 0.0       # log rho
qo = 0.0        # opacity
qol = 0.0       # log opacity (read nowhere)
qt = 0.0        # T
qp = 0.0        # P
qqdt = 0.0      # dln rho / dln T
qemu = 0.0
qd = 0.0        # rho
qfs = 0.0       # set to 1 by envint and qenv
