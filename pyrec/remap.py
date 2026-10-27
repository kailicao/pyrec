"""Python version of remap.f (YREC 09c3b9b): translates the newer,
"more intuitively named" namelist parameters into the arrays the rest of
the code uses."""

from pyrec import cross_mod
from pyrec import nmlbypass_mod
from pyrec import nmlpar_mod

# DATA: the reference S-factors that RATES and ENGEB assume (keV b).
S0_1_1_REF = 4.07e-22
S0_3_3_REF = 5.15e3
S0_3_4_REF = 5.4e-1
S0_1_12_REF = 1.45
S0_1_13_REF = 5.5
S0_1_14_REF = 3.32
S0_1_16_REF = 9.4
S0_PEP_REF = 3.4848e-6       # relative to pp
S0_1_BE7E_REF = 1.752e-10    # Be7+e-, in T9 units
S0_1_BE7P_REF = 0.0243
S0_HEP_REF = 8.0e-20
S0_1_15_C12ALP_REF = 6.75e4
S0_1_15_O16_REF = 6.4e1
# DATA: Adelberger et al. 1998 (Solar Fusion I) values; the others are
# commented out in the Fortran (KC 2025-05-31).
S0_1_1_A98 = 4.00e-22
S0_1_BE7P_A98 = 0.019
# First derivatives (b).
QS0E_1_1_A98 = 4.48e-24
QS0E_3_3_A98 = -4.1
QS0E_3_4_A98 = -3.0e-4
QS0E_1_12_A98 = 2.6e-3
QS0E_1_13_A98 = -7.8e-3
QS0E_1_14_A98 = -1.28e-2
QS0E_1_16_A98 = -2.4e-2
QS0E_1_BE7P_A98 = -1.33e-5
# Second derivatives.
QQS0EE_1_12_A98 = 8.3e-5
QQS0EE_1_13_A98 = 7.3e-4
QQS0EE_1_16_A98 = 5.7e-5
QQS0EE_1_BE7P_A98 = 7.22e-8


def _ratio(s0, s0p, ref_deriv, ref_s0):
    """(S'/S0) / (S'/S0)_A98, or 0 when S0 <= 0 (MHP 4/25)."""
    if s0 > 0.0:
        return (s0p / s0) / (ref_deriv / ref_s0)
    return 0.0


def remap():
    """With LNEWNUC, set SSTANDARD, the N15+p branching factors FC12, FO16
    and the derivative factors QS0E, QQS0EE from the individual
    cross sections; with LNEWVARS, copy the new tolerance names into
    HPTTOL, ATIME, LPTIME, DTWIND and DTDIF.

    Translates: remap.f#REMAP (09c3b9b)

    No arguments; nothing is returned. Reads nmlbypass_mod, sets
    cross_mod and nmlpar_mod.
    """
    nb = nmlbypass_mod
    nml = nmlpar_mod
    cr = cross_mod

    for i in range(1, 9):
        cr.qs0e[i] = 1.0
        cr.qqs0ee[i] = 1.0
    cr.fc12 = 1.0
    cr.fo16 = 1.0
    if nb.lnewnuc:
        for i in range(1, 18):
            nml.sstandard[i] = 1.0
        nml.sstandard[1] = nb.s0_1_1 / S0_1_1_REF
        nml.sstandard[2] = nb.s0_3_3 / S0_3_3_REF
        nml.sstandard[3] = nb.s0_3_4 / S0_3_4_REF
        nml.sstandard[4] = nb.s0_1_12 / S0_1_12_REF
        nml.sstandard[5] = nb.s0_1_13 / S0_1_13_REF
        nml.sstandard[6] = nb.s0_1_14 / S0_1_14_REF
        nml.sstandard[7] = nb.s0_1_16 / S0_1_16_REF
        # 8, 10, 11, 12 are resonant (not fit by a simple S0); 9, 13 unused.
        nml.sstandard[14] = nb.s0_pep / S0_PEP_REF
        nml.sstandard[15] = nb.s0_1_be7e / S0_1_BE7E_REF
        nml.sstandard[16] = nb.s0_1_be7p / S0_1_BE7P_REF
        nml.sstandard[17] = nb.s0_hep / S0_HEP_REF
        # The branching of N15+p.
        cr.fc12 = nb.s0_1_15_c12alp / S0_1_15_C12ALP_REF
        cr.fo16 = nb.s0_1_15_o16 / S0_1_15_O16_REF
        # First-derivative factors (multiply Q2 and Q3 in ENGEB and RATES).
        # Question (see the page): 2-7 divide by S0_1_1_A98, not their own.
        cr.qs0e[1] = _ratio(nb.s0_1_1, nb.s0p_1_1, QS0E_1_1_A98, S0_1_1_A98)
        cr.qs0e[2] = _ratio(nb.s0_3_3, nb.s0p_3_3, QS0E_3_3_A98, S0_1_1_A98)
        cr.qs0e[3] = _ratio(nb.s0_3_4, nb.s0p_3_4, QS0E_3_4_A98, S0_1_1_A98)
        cr.qs0e[4] = _ratio(nb.s0_1_12, nb.s0p_1_12, QS0E_1_12_A98,
                            S0_1_1_A98)
        cr.qs0e[5] = _ratio(nb.s0_1_13, nb.s0p_1_13, QS0E_1_13_A98,
                            S0_1_1_A98)
        cr.qs0e[6] = _ratio(nb.s0_1_14, nb.s0p_1_14, QS0E_1_14_A98,
                            S0_1_1_A98)
        cr.qs0e[7] = _ratio(nb.s0_1_16, nb.s0p_1_16, QS0E_1_16_A98,
                            S0_1_1_A98)
        cr.qs0e[8] = _ratio(nb.s0_1_be7p, nb.s0p_1_be7p, QS0E_1_BE7P_A98,
                            S0_1_BE7P_A98)
        # Second-derivative factors (multiply Q4 and Q5); 1-3 and 6 were
        # zeroed out in the 2003 version.
        cr.qqs0ee[4] = _ratio(nb.s0_1_12, nb.s0pp_1_12, QQS0EE_1_12_A98,
                              S0_1_1_A98)
        cr.qqs0ee[5] = _ratio(nb.s0_1_13, nb.s0pp_1_13, QQS0EE_1_13_A98,
                              S0_1_1_A98)
        cr.qqs0ee[7] = _ratio(nb.s0_1_16, nb.s0pp_1_16, QQS0EE_1_16_A98,
                              S0_1_1_A98)
        cr.qqs0ee[8] = _ratio(nb.s0_1_be7p, nb.s0pp_1_be7p,
                              QQS0EE_1_BE7P_A98, S0_1_BE7P_A98)

    # Overwrite the HPTTOL and ATIME vectors with the newer names.
    if nb.lnewvars:
        nml.hpttol[1] = nb.tol_dm_min           # min log mass spacing
        nml.hpttol[2] = nb.tol_dm_max           # max log mass spacing
        nml.hpttol[9] = nb.tol_dl_max           # max L/Lsurf spacing
        nml.hpttol[11] = nb.tol_dp_core_max     # max log P, radiative interior
        nml.hpttol[8] = nb.tol_dp_env_max       # max log P, convective envelope
        nml.hpttol[10] = nb.tol_dp_czbase_max   # max log P, near the CZ base
        nml.hpttol[7] = nb.tol_czbase_fine_width  # width of the fine zone
        nml.hpttol[5] = nb.tol_dx_max           # max X spacing (usually off)
        nml.hpttol[6] = nb.tol_dz_max           # max Z spacing (usually off)
        nml.hpttol[3] = nb.flag_dx              # X difference to flag
        nml.hpttol[4] = nb.flag_dz              # Z difference to flag
        nml.hpttol[12] = nb.flag_dw             # max log omega spacing
        nml.atime[1] = nb.time_core_min      # min central abundance for core dt
        nml.atime[2] = nb.time_dx_core_tot   # max delta X burned in the core
        nml.atime[3] = nb.time_dx_core_frac  # max fraction of X, core
        nml.atime[6] = nb.time_dx_total      # max Msun of X burned
        nml.atime[7] = nb.time_dx_shell      # max fraction of X, shell
        nml.atime[4] = nb.time_dy_core_tot   # max delta Y burned in the core
        nml.atime[5] = nb.time_dy_core_frac  # max fraction of Y, core
        nml.atime[12] = nb.time_dy_shell     # max fraction of Y, shell
        nml.atime[14] = nb.time_dy_total     # max Msun of Y burned in a shell
        # Max changes in T, P, R, L at any shell in a step.
        nml.atime[8] = nb.time_dt
        nml.atime[9] = nb.time_dp
        nml.atime[10] = nb.time_dr
        nml.atime[11] = nb.time_dl
        # Master flag for ATIME(8)-(11) (pre-MS only).
        nml.lptime = nb.lstruct_time
        nml.dtwind = nb.time_dw_global       # max change in omega per step
        nml.dtdif = nb.time_dw_mix           # the same, per mixing substep
        nml.atime[13] = nb.time_max_dt_frac  # max factor between steps
