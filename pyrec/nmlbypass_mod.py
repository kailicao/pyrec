"""Python version of nmlbypass_mod.f (YREC 09c3b9b): the namelist variables
that only set other quantities, formerly COMMON/NMLBYPASS/.

Translates: nmlbypass_mod.f#NMLBYPASS_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case, with the
module's initializers (the defaults of the PHYSICS and CONTROL namelists,
read by parmin). Routines that USE NMLBYPASS_MOD import this module and
write ``nmlbypass_mod.<name>``.
"""

flag_dx = 0.05
flag_dw = 0.1
flag_dz = 0.05
time_core_min = 0.001
time_dl = 0.02
time_dp = 0.04
time_dr = 0.02
time_dt = 0.02
time_dw_global = 0.08
time_dw_mix = 0.08
time_dx_core_frac = 0.5
time_dx_core_tot = 0.02
time_dx_shell = 0.1
time_dx_total = 0.0015
time_dy_core_frac = 0.5
time_dy_core_tot = 0.02
time_dy_shell = 0.1
time_dy_total = 0.0015
tol_czbase_fine_width = 0.0
tol_dl_max = 0.02
tol_dm_max = 0.08
tol_dm_min = 1e-08
tol_dp_core_max = 0.05
tol_dp_czbase_max = 0.05
tol_dp_env_max = 0.05
tol_dx_max = 1.0
tol_dz_max = 1.0
time_max_dt_frac = 1.5
lstruct_time = False
lnewvars = False
s0_1_1 = 4.01e-22
s0_3_3 = 5210.0
s0_3_4 = 0.56
s0_1_12 = 1.34
s0_1_13 = 7.6
s0_1_14 = 1.66
s0_1_16 = 10.6
s0_pep = 3.5734e-06
s0_1_be7e = 1.7709e-10
s0_1_be7p = 0.0208
s0_hep = 8.6e-20
s0_1_15_c12alp = 73000.0
s0_1_15_o16 = 36.0
s0p_1_1 = 4.49e-24
s0p_3_3 = -4.9
s0p_3_4 = -0.00036
s0p_1_12 = 0.0026
s0p_1_13 = -0.00783
s0p_1_14 = -0.0033
s0p_1_16 = -0.054
s0pp_1_12 = 8.3e-05
s0pp_1_13 = 0.000729
s0pp_1_16 = 0.0
s0p_1_be7p = -3.12e-05
s0pp_1_be7p = -2.288e-07
xsli6 = 5500.0
xsli7 = 52.0
xsbe91 = 1.1
xsbe92 = 15000.0
xsbe93 = 15000.0
lxli6 = False
lxli7 = False
lxbe91 = False
lxbe92 = False
lxbe93 = False
npenv = 0  # no initializer
lnewnuc = False
lenvg = False
isetmix = 0
isetiso = 0
r14_15 = 411.9
r16_17 = 471.4
zxmix = 0.02292
xb10_ini = 1.002e-09
xb11_ini = 4.405e-09
pmma = 2.0
pmmb = 1.0
pmmc = 2.0
pmmd = 0.0
pmmm = 0.5
pmmjd = 1.32e+30
pmmmd = 1270000000000.0
fo = 1.0
