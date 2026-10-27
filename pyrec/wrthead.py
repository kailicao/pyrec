"""Python version of wrthead.f (YREC 09c3b9b): the header lines of the
output files at the start of each run."""

from pyrec import iounits_mod
from pyrec import label_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod

_RUN_KIND = {1: 'EVOLVE  ', 2: 'RESCALE ', 3: 'RESCALE&EVOLVE '}

# FORMAT 1504 (ITRVER = 0, G. Somers 11/14): one line of column names,
# 16 records of the Fortran joined into one; it does not start with '#'.
_COLUMNS_1504 = (
    '     Step    Shls         Age_gyr       LogL_lsun       LogR_rsun           Log_g'
    '        log_Teff        Mco_core         Mco_env     Rco_env     Tco_env     Dco_env'
    '     Pco_env     Oco_env        LogT_cen        LogD_cen        logP_cen        Beta_cen'
    '         Eta_cen           X_cen           Y_cen           Z_cen        ppI_lsun'
    '       ppII_lsun      ppIII_lsun        CNO_lsun         3a_lsun        HeC_lsun'
    '      Egrav_lsun       Neut_lsun    Cl_snu    Ga_snu   pp_neut  pep_neut  hep_neut'
    '  Be7_neut   B8_neut  N13_neut  O15_neut  F17_neut     diag1     diag2         He3_cen'
    '         C12_cen         C13_cen         N14_cen         N15_cen         O16_cen'
    '         O17_cen         O18_cen         He3_sur         C12_sur         C13_sur'
    '         N14_sur         N15_sur         O16_sur         O17_sur         O18_sur'
    '          H2_sur         Li6_sur         Li7_sur         Be9_sur           X_sur'
    '           Y_sur           Z_sur         Z_X_sur            Jtot      KE_rot_tot'
    '           I_tot            I_cz       Omega_sur       Omega_cen      Prot_sur_d'
    '        Vrot_kms         TauCZ_s    MHshell_base     MHshell_mid     MHshell_top'
    '    RHshell_base     RHShell_mid     RHshell_top       logP_phot       Mass_msun'
)
_FMT_1504 = '# \n' + _COLUMNS_1504 + '\n# '

# FORMAT 1505 (ITRVER = 1 and 2). Two of its lines do not start with '#'.
_FMT_1505 = '\n'.join([
    '# Model #, shells, AGE, log(L/Lsun), log(R/Rsun), log(g),'
    ' log(Teff), Mconv. core, Mconv. env.)',
    '# Central: log(T), log(RHO), log(P), BETA, ETA, X, Y, Z',
    '# Luminosity: ppI, ppII, ppIII, CNO, triple-alpha,'
    ' He-C, gravity, neutrinos (old)',
    ' Neutrinos (1E10 erg/CM^2 at earth): pp, pep, hep, Be7,'
    ' B8, N13, O15, F17',
    '# Central Abundances: He3, C12, C13, N14, N15, O16,'
    ' O17, O18',
    '# Surface Abundances: He3, C12, C13, N14, N15, O16,'
    ' O17,O18',
    '    "        " cont: H2, Li6, Li7, Be9',
    '# H shell loc: mass frac-base, midpoint, top; radius frac-'
    'base, midpoint, top',
])

# FORMAT 1510 (ITRVER = 2): the rotation columns, without '#'.
_FMT_1510 = ' Jtot, K.E. Rotation, OMEGAsurf, OMEGAcenter'

# FORMAT 1515 (ITRVER = 3).
_FMT_1515 = ('#Model #, shells, AGE, log(L/Lsun), log(R/Rsun), log(g),'
             ' log(Teff), Mconv. core, Mconv. env., % Grav Energy, X env')


def wrthead(smass, cmixl, nk):
    """Write the run's line on IOWR and, for a run that reads its starting
    model (LFIRST(NK)) with LTRACK on, the track file's version line and
    column header on ITRACK, in the layout chosen by ITRVER.

    Translates: wrthead.f#WRTHEAD (09c3b9b)

    Arguments: the stellar mass smass (solar masses), the mixing length
    cmixl and the run index nk. Nothing is returned.
    """
    n = nmlpar_mod
    r = runloop_mod
    xenv0 = label_mod.xenv0
    zenv0 = label_mod.zenv0

    # FORMATs 47, 48, 49: the same line but for the run's kind.
    kind = r.iresca[nk]
    if kind in _RUN_KIND:
        iounits_mod.write(
            iounits_mod.iowr,
            f"\n RUN={nk:2d} {_RUN_KIND[kind]} X={xenv0:8.6f}"
            f" Z={zenv0:8.6f} CMIXL={cmixl:8.6f}"
            f" NO.MODS={r.nmodls[nk]:5d}")

    if n.ltrack and r.lfirst[nk]:
        # ITRVER identifies the version of the track file: change it when
        # the track file changes (and WRTOUT with it).
        iounits_mod.write(
            iounits_mod.itrack,
            f"#Version={n.itrver:3d}  Mtot/Msun ={smass:16.8E}"
            f"  Initial: X ={xenv0:16.8E} Z ={zenv0:16.8E}"
            f"  Mix. length ={cmixl:16.8E}")
        if n.itrver == 0:
            # Question (see the page): the comment says to uncomment a
            # block, but FORMAT 1504's second version is active; its line
            # of column names does not start with '#'.
            iounits_mod.write(iounits_mod.itrack, _FMT_1504)
        elif n.itrver == 1:
            iounits_mod.write(iounits_mod.itrack, _FMT_1505)
        elif n.itrver == 2:
            iounits_mod.write(iounits_mod.itrack, _FMT_1505)
            iounits_mod.write(iounits_mod.itrack, _FMT_1510)
        elif n.itrver == 3:
            iounits_mod.write(iounits_mod.itrack, _FMT_1515)
