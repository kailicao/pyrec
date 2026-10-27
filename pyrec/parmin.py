"""Python version of parmin.f (YREC development, 09c3b9b): a run's
parameters, its input and output files, and its run cards."""

import math
import os
import re
import sys
from types import SimpleNamespace

import numpy as np

from pyrec import acdpth_mod
from pyrec import ccout_mod
from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import intpar_mod
from pyrec import iounits_mod
from pyrec import label_mod      # USEd by the Fortran; no name read here
from pyrec import newopac_mod    # USEd by the Fortran; no name read here
from pyrec import nmlbypass_mod
from pyrec import nmlpar_mod
from pyrec import runloop_mod
from pyrec.remap import remap

_this = sys.modules[__name__]

# PARMIN's own locals with DATA defaults, kept between calls (SAVE). Arrays
# have index 0 unused, as everywhere in pyrec. Character values are kept
# without Fortran's trailing blanks.
kindrn = np.ones(51, dtype=int)            # DATA KINDRN/50*1/
rsclm = np.zeros(51)
rsclx = np.zeros(51)
rsclz = np.zeros(51)
rsclcm = np.zeros(51)
aid = ['', 'HE3', 'C12', 'C13', 'N14', 'N15', 'O16', 'O17', 'O18', 'H2',
       'LI6', 'LI7', 'BE9']
anewcp = ''
atmp = ''
# MHP 10/24: the CNO mixture and isotope sets. AMIXT lists the preset
# mixtures, each with its Z/X and C, N, O fractions of the metals: GS98,
# Asplund et al. 2021, and Magg et al. (met, phot).
# Question (see the page): the year of the Magg et al. mixtures.
amix = 'GS98'
aiso = 'L21'
amixt = ['', 'GS98', 'AAG21', 'M22P', 'M22M']
zxmixt = np.array([0.0, 0.02292, 0.0187, 0.0225, 0.0226])
frac_ct = np.array([0.0, 0.172148, 0.184156, 0.19239, 0.191425])
frac_nt = np.array([0.0, 0.050426, 0.050344, 0.059012, 0.058716])
frac_ot = np.array([0.0, 0.468195, 0.416592, 0.415631, 0.413545])
alfa = 1.5          # the old wind law's index
fk = 1.0            # its fudge factor
wmax_sun = 1000.0   # MHP 8/17: saturation threshold in solar units
awind = 'K97'       # the Matt et al. presets: 'K97', 'V13' or 'CUS'
fsound = 'Csound.dat'
descrip = ['', '', '']
fdebug = ffirst = flast = ''
fmilne = fmodpt = fopale06 = fscomp = ''
fshort = fsnu = fstch = fstor = ftrack = fcondopacp = ''
fcalcad = ''        # FCLCD
yrec1 = yrec2 = ''
cmd = ''
sum_frac = 0.0
last_slash_idx = 0  # Question (see the page): unset when FTRACK has no '/'
mrk = 0
ncmxck = 0
lcmxbad = False

# The three namelists, member by member in the Fortran's order, each with
# where the variable lives: 'nml' NMLPAR_MOD, 'byp' NMLBYPASS_MOD, 'con'
# CONSTANTS_MOD, 'ipar' INTPAR_MOD, 'loc' a local of PARMIN (this file),
# 'arg' one of PARMIN's arguments.
_CONTROL = (
    ('nml', 'cmixla'), ('nml', 'calsolage'), ('nml', 'calsolzx'),
    ('loc', 'descrip'), ('nml', 'endage'), ('arg', 'fpurez'),
    ('loc', 'flast'), ('loc', 'ffirst'), ('arg', 'ffermi'), ('loc', 'fdebug'),
    ('loc', 'ftrack'), ('loc', 'fshort'), ('loc', 'fstch'), ('loc', 'fmilne'),
    ('loc', 'fmodpt'), ('loc', 'fstor'), ('arg', 'fpmod'), ('arg', 'fpatm'),
    ('arg', 'fpenv'), ('loc', 'fsnu'), ('loc', 'fscomp'), ('arg', 'fatm'),
    ('loc', 'fsound'), ('arg', 'fallard'), ('arg', 'fscvh'),
    ('arg', 'fscvhe'), ('arg', 'fscvz'), ('arg', 'fliv95'), ('nml', 'ipver'),
    ('nml', 'itrver'), ('loc', 'kindrn'), ('nml', 'ldebug'), ('nml', 'lcorr'),
    ('nml', 'lmilne'), ('nml', 'ltrack'), ('nml', 'lstore'),
    ('nml', 'lfirst'), ('nml', 'lstpch'), ('nml', 'lstch'), ('nml', 'lstatm'),
    ('nml', 'lstenv'), ('nml', 'lstmod'), ('nml', 'lstphys'),
    ('nml', 'lstrot'), ('nml', 'lpulse'), ('nml', 'lzramp'), ('nml', 'lteff'),
    ('nml', 'lcalst'), ('nml', 'lpurez'), ('nml', 'lrwsh'),
    ('nml', 'lsenv0a'), ('nml', 'lpout'), ('nml', 'lcals'),
    ('nml', 'lcalsolzx'), ('nml', 'lopal95'), ('nml', 'npoint'),
    ('byp', 'npenv'), ('nml', 'nprtmod'), ('nml', 'nprtpt'),
    ('nml', 'numrun'), ('nml', 'nmodls'), ('nml', 'poa'), ('nml', 'pob'),
    ('nml', 'poc'), ('nml', 'pomax'), ('loc', 'rsclm'), ('loc', 'rsclx'),
    ('loc', 'rsclz'), ('loc', 'rsclcm'), ('nml', 'rsclzc'),
    ('nml', 'rsclzm1'), ('nml', 'rsclzm2'), ('nml', 'setdt'),
    ('nml', 'senv0a'), ('nml', 'steff'), ('nml', 'sr'), ('nml', 'tolr'),
    ('nml', 'toll'), ('nml', 'tolz'), ('nml', 'xenv0a'), ('nml', 'xls'),
    ('nml', 'xlstol'), ('nml', 'zenv0a'), ('nml', 'zopal951'),
    ('loc', 'fcondopacp'), ('loc', 'fopale06'), ('arg', 'falex06'),
    ('nml', 'lalex06'), ('nml', 'end_dcen'), ('nml', 'end_xcen'),
    ('nml', 'end_ycen'), ('byp', 'isetmix'), ('byp', 'isetiso'),
    ('loc', 'amix'), ('loc', 'aiso'), ('nml', 'frac_c'), ('nml', 'frac_n'),
    ('nml', 'frac_o'), ('nml', 'r12_13'), ('byp', 'r14_15'),
    ('byp', 'r16_17'), ('nml', 'r16_18'), ('byp', 'zxmix'),
    ('nml', 'xh2_ini'), ('nml', 'xhe3_ini'), ('nml', 'xli6_ini'),
    ('nml', 'xli7_ini'), ('nml', 'xbe9_ini'), ('byp', 'xb10_ini'),
    ('byp', 'xb11_ini'),
)

_PHYSICS = (
    ('nml', 'atmmin'), ('nml', 'atmbeg'), ('nml', 'atmerr'),
    ('nml', 'atmmax'), ('nml', 'atmd0'), ('loc', 'anewcp'), ('loc', 'atmp'),
    ('nml', 'acfpft'), ('nml', 'atime'), ('nml', 'alphac'), ('nml', 'alphae'),
    ('loc', 'alfa'), ('nml', 'alpham'), ('nml', 'atmstp'), ('nml', 'abstol'),
    ('nml', 'betac'), ('nml', 'cmin'), ('con', 'clsun'), ('con', 'crsun'),
    ('nml', 'dpenv'), ('nml', 'dtdif'), ('nml', 'dtwind'), ('nml', 'djok'),
    ('nml', 'dt_gs'), ('nml', 'enverr'), ('nml', 'envmax'), ('nml', 'envmin'),
    ('nml', 'envbeg'), ('nml', 'envstp'), ('nml', 'etadh0'),
    ('nml', 'etadh1'), ('nml', 'fcorr0'), ('nml', 'fcorri'), ('loc', 'fk'),
    ('nml', 'fw'), ('nml', 'fc'), ('byp', 'fo'), ('nml', 'fmu'),
    ('nml', 'fes'), ('nml', 'fcore'), ('nml', 'fgsf'), ('nml', 'fss'),
    ('nml', 'fesc'), ('nml', 'fssc'), ('nml', 'fgsfc'), ('nml', 'fgry'),
    ('nml', 'fgrz'), ('nml', 'grtol'), ('nml', 'htoler'), ('nml', 'hpttol'),
    ('nml', 'itfp1'), ('nml', 'itfp2'), ('ipar', 'imax'), ('nml', 'itdif1'),
    ('nml', 'itdif2'), ('nml', 'ies'), ('nml', 'igsf'), ('nml', 'imu'),
    ('nml', 'ilambda'), ('nml', 'kttau'), ('nml', 'kemmax'), ('nml', 'lvfc'),
    ('nml', 'ldifad'), ('nml', 'lnoj'), ('nml', 'lnewdif'), ('nml', 'ldify'),
    ('nml', 'ldifz'), ('nml', 'ldifli'), ('nml', 'lsnu'), ('nml', 'ldh'),
    ('nml', 'lnewcp'), ('nml', 'lkuthe'), ('nml', 'lovstc'),
    ('nml', 'lovste'), ('nml', 'lovstm'), ('nml', 'lovmax'),
    ('nml', 'lexcom'), ('nml', 'lrot'), ('nml', 'lnew0'), ('nml', 'linstb'),
    ('nml', 'lwnew'), ('nml', 'ljdot0'), ('nml', 'lptime'), ('nml', 'ladov'),
    ('nml', 'ltrist'), ('byp', 'lenvg'), ('nml', 'lnulos1'),
    ('nml', 'lthoul'), ('nml', 'lthoulfit'), ('nml', 'lcore'),
    ('nml', 'lsemic'), ('nml', 'lnews'), ('nml', 'mcore'), ('nml', 'niter1'),
    ('nml', 'niter2'), ('nml', 'niter3'), ('nml', 'niter4'), ('ipar', 'nuse'),
    ('nml', 'niter_gs'), ('nml', 'optol'), ('nml', 'rcrit'),
    ('nml', 'reltol'), ('ipar', 'stolr0'), ('nml', 'tcut'), ('nml', 'tscut'),
    ('nml', 'tenv0'), ('nml', 'tenv1'), ('nml', 'tgcut'), ('nml', 'tridt'),
    ('nml', 'tridl'), ('nml', 'tollaol'), ('nml', 'vnew'), ('nml', 'walpcz'),
    ('nml', 'wnew'), ('nml', 'weakscreening'), ('nml', 'xnewcp'),
    ('nml', 'xmin'), ('nml', 'ymin'), ('nml', 'tmolmin'), ('nml', 'tmolmax'),
    ('nml', 'sstandard'), ('nml', 'lscv'), ('nml', 'ldisk'), ('nml', 'tdisk'),
    ('nml', 'pdisk'), ('nml', 'wmax'), ('nml', 'lsolid'), ('nml', 'impjmod'),
    ('nml', 'dmdt0'), ('nml', 'fczdmdt'), ('nml', 'ftotdmdt'),
    ('nml', 'compacc'), ('nml', 'creim'), ('nml', 'lreimer'),
    ('nml', 'lmdot'), ('nml', 'lcondopacp'), ('nml', 'lopale06'),
    ('nml', 'lnumderiv'), ('nml', 'alatm_feh'), ('nml', 'alatm_alpha'),
    ('nml', 'laltptau100'), ('nml', 'cstmixing'), ('nml', 'cstdiffmix'),
    ('nml', 'lsolwind'), ('nml', 'lmwind'), ('nml', 'lrossby'),
    ('nml', 'lbscale'), ('loc', 'awind'), ('byp', 'pmma'), ('byp', 'pmmb'),
    ('byp', 'pmmc'), ('byp', 'pmmd'), ('byp', 'pmmm'), ('byp', 'pmmjd'),
    ('byp', 'pmmmd'), ('nml', 'pmmsolp'), ('nml', 'pmmsolw'),
    ('nml', 'pmmsoltau'), ('nml', 'lcodm'), ('nml', 'codm'),
    ('loc', 'wmax_sun'), ('byp', 'xsli6'), ('byp', 'xsli7'),
    ('byp', 'xsbe91'), ('byp', 'xsbe92'), ('byp', 'xsbe93'), ('byp', 'lxli6'),
    ('byp', 'lxli7'), ('byp', 'lxbe91'), ('byp', 'lxbe92'), ('byp', 'lxbe93'),
    ('nml', 'sli6'), ('nml', 'sli7'), ('nml', 'sbe91'), ('nml', 'sbe92'),
    ('nml', 'sbe93'), ('byp', 'lnewnuc'), ('nml', 'spotf'), ('nml', 'spotx'),
    ('nml', 'lsdepth'), ('byp', 's0_1_1'), ('byp', 's0_3_3'),
    ('byp', 's0_3_4'), ('byp', 's0_1_12'), ('byp', 's0_1_13'),
    ('byp', 's0_1_14'), ('byp', 's0_1_16'), ('byp', 's0_pep'),
    ('byp', 's0_1_be7e'), ('byp', 's0_1_be7p'), ('byp', 's0_hep'),
    ('byp', 's0_1_15_c12alp'), ('byp', 's0_1_15_o16'), ('byp', 's0p_1_1'),
    ('byp', 's0p_3_3'), ('byp', 's0p_3_4'), ('byp', 's0p_1_12'),
    ('byp', 's0p_1_13'), ('byp', 's0p_1_14'), ('byp', 's0p_1_16'),
    ('byp', 's0pp_1_12'), ('byp', 's0pp_1_13'), ('byp', 's0pp_1_16'),
    ('byp', 's0p_1_be7p'), ('byp', 's0pp_1_be7p'), ('byp', 'flag_dx'),
    ('byp', 'flag_dw'), ('byp', 'flag_dz'), ('byp', 'lstruct_time'),
    ('byp', 'time_core_min'), ('byp', 'time_dl'), ('byp', 'time_dp'),
    ('byp', 'time_dr'), ('byp', 'time_dt'), ('byp', 'time_dw_global'),
    ('byp', 'time_dw_mix'), ('byp', 'time_dx_core_frac'),
    ('byp', 'time_dx_core_tot'), ('byp', 'time_dx_shell'),
    ('byp', 'time_dx_total'), ('byp', 'time_dy_core_frac'),
    ('byp', 'time_dy_core_tot'), ('byp', 'time_dy_shell'),
    ('byp', 'time_dy_total'), ('byp', 'tol_czbase_fine_width'),
    ('byp', 'tol_dl_max'), ('byp', 'tol_dm_max'), ('byp', 'tol_dm_min'),
    ('byp', 'tol_dp_core_max'), ('byp', 'tol_dp_czbase_max'),
    ('byp', 'tol_dp_env_max'), ('byp', 'tol_dx_max'), ('byp', 'tol_dz_max'),
    ('byp', 'time_max_dt_frac'), ('byp', 'lnewvars'), ('nml', 'lnewtcz'),
    ('nml', 'lcalcenv'),
)

_DERIVED = (
    ('con', 'cln'), ('con', 'cc13'), ('con', 'cc23'), ('nml', 'constfactor'),
    ('nml', 'excen'), ('nml', 'exl'), ('nml', 'exm'), ('nml', 'exmd'),
    ('nml', 'expr'), ('nml', 'exr'), ('nml', 'extau'), ('nml', 'exw'),
    ('nml', 'wmax'), ('nml', 'ldify'), ('nml', 'ldifz'), ('nml', 'lthoul'),
    ('nml', 'lthoulfit'), ('nml', 'lnew0'), ('nml', 'lnewdif'),
    ('nml', 'ilambda'), ('nml', 'sli6'), ('nml', 'sli7'), ('nml', 'sbe91'),
    ('nml', 'sbe92'), ('nml', 'sbe93'), ('nml', 'frac_c'), ('nml', 'frac_n'),
    ('nml', 'frac_o'), ('byp', 'zxmix'), ('loc', 'sum_frac'),
    ('nml', 'lisotope'), ('nml', 'lmixture'), ('nml', 'inewcp'),
    ('nml', 'lrel'), ('nml', 'kttau0'), ('nml', 'lttau'),
)

# The .store file's header key (format 1012).
_STORE_KEY = (
    "# Header Key\n"
    "# ModType    ModNum    #Shells    M/Msun    log(Teff)    log(L/Lsun)"
    "    log(M/gram)    Age/Gyr    Timestep/yr    log(M_inner/gram)"
    "    log(M_outer/gram)\n"
    "# JCORE  JENV  CMIXL  EOS  ATM  ALOK HIK  LPUREZ  COMPMIX  LEXCOM"
    "  LDIFY  LDIFZ  LSEMIC  LOVSTC  LOVSTE  LOVSTM  LROT  LINSTB  LJDOT0"
    "  LDISK  TDISK  PDISK  WMAX  LSTORE  LSTATM  LSTENV  LSTMOD  LSTPHYS"
    "  LSTROT\n"
    "# PPI_lum    PPII_lum    PPIII_lum    CNO_lum    3He_lum    "
    "OldNeu_lum    Grav_lum\n"
    "# ENV1-3   log(Teff)    log(L/Lsun)    P_base    T_base    R_base"
    "    MatrixElements\n")

# The atmosphere choice by KTTAU (formats 197, 198, 1999, 1888, 1889,
# 1887).
_KTTAU_TEXT = {
    0: " USING EDDINGTON T-TAU RELATION.",
    1: " USING KRISHNA-SWAMY T-TAU RELATION.",
    2: " USING HARVARD-SMITHSONIAN REFERENCE ATMOSPHERE",
    3: " USING KURUCZ ATMOSPHERE TABLE",
    4: " USING ALLARD ATMOSPHERE TABLE",
    5: " USING KURUCZ/CASTELLI ATMOSPHERE TABLE",
}


def _l(flag):
    """A logical in an L format: T or F."""
    return 'T' if flag else 'F'


def _owners(args):
    return {'nml': nmlpar_mod, 'byp': nmlbypass_mod, 'con': constants_mod,
            'ipar': intpar_mod, 'loc': _this, 'arg': args}


def _value(token, like):
    """Convert one namelist value to the type of the variable it sets."""
    if isinstance(like, (bool, np.bool_)):
        t = token.upper().lstrip('.')
        return t.startswith('T')
    if isinstance(like, str):
        if token[:1] in "'\"":
            return token[1:-1].replace(token[0] * 2, token[0])
        return token
    number = token.replace('d', 'e').replace('D', 'E')
    if isinstance(like, (int, np.integer)):
        return int(float(number))
    return float(number)


def _fortran_order(shape, start):
    """The indices of an array (index 0 unused) from `start` on, in
    Fortran's order: the first index varies fastest."""
    idx = list(start)
    while True:
        yield tuple(idx)
        for d in range(len(idx)):
            idx[d] += 1
            if idx[d] < shape[d]:
                break
            idx[d] = 1
        else:
            return


def _read_namelist(path, group, members, owners):
    """READ(UNIT, NML=group): a minimal reader for the namelist files YREC
    uses (name = values, name(i) = values, repeat counts r*v, quoted
    strings, logicals, D exponents, '!' comments). Each value is stored
    where its variable lives."""
    table = {name: (owners[where], name) for where, name in members}
    with open(path) as f:
        # Drop '!' comments, outside quoted strings.
        text = ''.join(re.sub(r"^((?:[^'\"!]|'[^']*'|\"[^\"]*\")*)!.*",
                              r"\1", line) for line in f)
    start = re.search(r'[&$]' + group + r'\b', text, re.IGNORECASE)
    if start is None:
        raise SystemExit(f"namelist /{group}/ not found in {path}")
    body = text[start.end():]
    tokens = re.findall(r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"|"
                        r"[A-Za-z_]\w*\s*(?:\([^)]*\))?\s*=|/|"
                        r"[&$]END\b|[^\s,'\"/=]+", body, re.IGNORECASE)
    target = None
    for tok in tokens:
        if tok == '/' or tok.upper().lstrip('&$') == 'END':
            break
        if tok.endswith('='):
            m = re.match(r'(\w+)\s*(?:\(([^)]*)\))?', tok)
            owner, attr = table[m.group(1).lower()]
            current = getattr(owner, attr)
            if m.group(2):
                first = tuple(int(i) for i in m.group(2).split(','))
            else:
                first = None
            if isinstance(current, np.ndarray):
                place = _fortran_order(current.shape,
                                       first or (1,) * current.ndim)
            elif isinstance(current, list):
                place = iter(range(first[0] if first else 1, len(current)))
            else:
                place = None
            target = (owner, attr, current, place)
            continue
        owner, attr, current, place = target
        count, _, item = tok.rpartition('*') if re.match(r'\d+\*', tok) \
            else ('1', '', tok)
        for _ in range(int(count)):
            if place is None:
                setattr(owner, attr, _value(item, current))
            else:
                where = next(place)
                like = current[where]
                current[where] = _value(item, like)


def _format_value(v):
    if isinstance(v, (bool, np.bool_)):
        return _l(v)
    if isinstance(v, str):
        return "'" + v + "'"
    return repr(v.item() if isinstance(v, np.generic) else v)


def _write_namelist(unit, group, members, owners):
    """WRITE(UNIT, NML=group): every member with its value (the layout
    approximates the compiler's)."""
    lines = [f"&{group}"]
    for where, name in members:
        v = getattr(owners[where], name)
        if isinstance(v, np.ndarray):
            sub = v[(slice(1, None),) * v.ndim].ravel(order='F')
            text = ', '.join(_format_value(x) for x in sub)
        elif isinstance(v, list):
            text = ', '.join(_format_value(x) for x in v[1:])
        else:
            text = _format_value(v)
        lines.append(f" {name.upper()}={text},")
    lines.append("/")
    iounits_mod.write(unit, "\n".join(lines))


def _both(text):
    """WRITE(*, fmt) and WRITE(ISHORT, fmt) of the same record."""
    print(text)
    iounits_mod.write(iounits_mod.ishort, text)


def _screen_and_short(text):
    """WRITE(IOWR, fmt) and WRITE(ISHORT, fmt) of the same record."""
    iounits_mod.write(iounits_mod.iowr, text)
    iounits_mod.write(iounits_mod.ishort, text)


def parmin(falex06, fallard, fatm, ffermi, fliv95, fpatm, fpenv, fpmod,
           fpurez, fscvh, fscvhe, fscvz, yrecver, githash):
    """Read the run's parameters from the CONTROL and PHYSICS namelists into
    the modules, set the unit numbers, expand the file names' placeholders,
    open the output files and some tables, derive settings (wind law,
    mixtures, settling switches), write the settings to the .short file,
    and interpret the run cards.

    Translates: parmin.f#PARMIN (09c3b9b)

    Arguments: the 12 table file names (set from the CONTROL namelist; any
    value), and yrecver, githash, the code's version and git hash.

    Returns (falex06, fallard, fatm, ffermi, fliv95, fpatm, fpenv, fpmod,
    fpurez, fscvh, fscvhe, fscvz): the file names as read and expanded.
    """
    global anewcp, atmp, sum_frac, last_slash_idx, mrk, ncmxck, lcmxbad
    global yrec1, yrec2, cmd, fcalcad
    me = _this
    io = iounits_mod
    nml = nmlpar_mod
    byp = nmlbypass_mod
    c = constants_mod
    args = SimpleNamespace(falex06=falex06, fallard=fallard, fatm=fatm,
                           ffermi=ffermi, fliv95=fliv95, fpatm=fpatm,
                           fpenv=fpenv, fpmod=fpmod, fpurez=fpurez,
                           fscvh=fscvh, fscvhe=fscvhe, fscvz=fscvz)
    owners = _owners(args)

    # All logical unit numbers.
    io.iowr = 6          # output: status file (the screen)
    io.ilast = 11        # output: last model
    io.ifirst = 12       # input: first model
    io.irun = 13         # input: PHYSICS namelist
    io.istand = 14       # input: CONTROL namelist
    io.ifermi = 15       # input: Fermi tables
    io.idebug = 18       # output: debugging
    io.itrack = 19       # output: track
    io.ishort = 20       # output: all diagnostic information
    io.imilne = 21       # output: Milne invariant variables
    io.imodpt = 22       # output: shell by shell information
    io.istor = 23        # output: saved models
    io.iopmod = 24       # output: pulsation, interior
    io.iopenv = 25       # output: pulsation, envelope
    io.iopatm = 26       # output: pulsation, atmosphere
    io.isound = 500      # output: sound-speed profile (WRTMOD, LSOUND)
    io.iliv95 = 48       # input: OPAL95 opacity table
    io.isnu = 33         # output: neutrino fluxes
    io.ioatm = 38        # input: Kurucz atmosphere table
    io.iopale = 49       # input: OPAL equation of state
    io.iopurez = 62      # input: pure C and N opacities
    io.ioatma = 66       # input: Allard model atmospheres
    io.iscvh = 72        # input: SCV EOS tables
    io.iscvhe = 73
    io.iscvz = 74
    io.icondopacp = 75   # input: Potekhin conductive opacities
    # Default [Fe/H] and alpha for the new Allard atmospheres.
    nml.alatm_feh = 0.0
    nml.alatm_alpha = 0.0
    nml.laltptau100 = False
    io.ialex06 = 90      # input: Alexander 2006 opacities

    print(' ')
    print(f" Yale Rotating Evolution Code - YREC, v{yrecver.rstrip()} "
          f"({githash.rstrip()})")
    # VERFMT, the version line (format built at run time, 315).
    verline = f"# YREC v{yrecver.rstrip()} ({githash.rstrip()})"

    # JVS 02/11: the namelist files from the command line, or the defaults.
    # Question (see the page): GETARG is a compiler extension.
    yrec1 = sys.argv[1] if len(sys.argv) > 1 else ''
    if yrec1[:2].strip() == '':
        yrec1 = 'yrec8.nml1'
    print(' ')
    print(f" CONTROL namelist :  {yrec1.rstrip()}")
    yrec2 = sys.argv[2] if len(sys.argv) > 2 else ''
    if yrec2[:2].strip() == '':
        yrec2 = 'yrec8.nml2'
    print(f" PHYSICS namelist :  {yrec2.rstrip()}")
    # OPEN, READ(NML=...) and CLOSE of ISTAND and IRUN.
    _read_namelist(yrec1, 'CONTROL', _CONTROL, owners)
    _read_namelist(yrec2, 'PHYSICS', _PHYSICS, owners)
    # MHP 8/14: map the newer inputs onto the vectors the code uses.
    remap()
    # MHP 06/13: remember the original atmosphere choice.
    nml.kttau0 = nml.kttau
    nml.lttau = False

    # Expand a leading placeholder in every file name.
    for where, name in (('arg', 'falex06'), ('arg', 'fallard'),
                        ('arg', 'fatm'), ('loc', 'fcondopacp'),
                        ('arg', 'ffermi'), ('loc', 'ffirst'),
                        ('loc', 'flast'), ('arg', 'fliv95'),
                        ('loc', 'fmodpt'), ('loc', 'fopale06'),
                        ('arg', 'fpatm'), ('arg', 'fpenv'), ('arg', 'fpmod'),
                        ('arg', 'fpurez'), ('loc', 'fscomp'),
                        ('arg', 'fscvh'), ('arg', 'fscvhe'),
                        ('arg', 'fscvz'), ('loc', 'fshort'),
                        ('loc', 'fsnu'), ('loc', 'fsound'), ('loc', 'fstor'),
                        ('loc', 'ftrack')):
        owner = owners[where]
        setattr(owner, name, expand_value(getattr(owner, name)))

    # Create the output folder, FTRACK up to its last '/'.
    for i in range(len(me.ftrack.rstrip()), 0, -1):
        if me.ftrack[i - 1] == '/':
            last_slash_idx = i
            break
    # label 1250. Question (see the page): no '/' leaves LAST_SLASH_IDX unset.
    cmd = 'mkdir -p "' + me.ftrack[:last_slash_idx] + '"'
    print(f" OUTPUT placed in :  {me.ftrack[:last_slash_idx]}")
    print(' ')
    os.system(cmd)

    # JVS 02/11: the acoustic-depth output, named after the .short file.
    if acdpth_mod.lacout:
        io.iclcd = 91
        mrk = me.fshort.find('short') + 1      # INDEX: 0 if absent
        # Possible bug (update-findings.md): no 'short' in FSHORT gives
        # plain 'calcad' in the current folder.
        fcalcad = me.fshort[:max(mrk - 1, 0)] + 'calcad'
        io.open(io.iclcd, fcalcad, status='UNKNOWN')
        io.write(io.iclcd, ' Acoustic depth calculation output file')
        # The second string continues onto a continuation line in the
        # Fortran; its blanks are kept as in the source.
        io.write(io.iclcd, ' age (Gyr),radius(cm),1/sound speed(s/cm),'
                 'radius (CZ), 1/cs (CZ)       delad,gamma1,P, T, X')

    # G Somers 6/14: Li/Be cross sections given by the user, as factors on
    # the defaults (Li6 5.5 MeV b; Li7 52 keV b; Be9 1.1, 15000, 15000 keV b).
    if byp.lxli6:
        nml.sli6 = byp.xsli6 / 5.5e3
    if byp.lxli7:
        nml.sli7 = byp.xsli7 / 5.2e1
    if byp.lxbe91:
        nml.sbe91 = byp.xsbe91 / 1.1
    if byp.lxbe92:
        nml.sbe92 = byp.xsbe92 / 1.5e4
    if byp.lxbe93:
        nml.sbe93 = byp.xsbe93 / 1.5e4
    # The OPAL 2006 EOS table.
    if nml.lopale06:
        io.open(io.iopale, me.fopale06, status='OLD')

    # The stored models, with their header if LSTORE.
    io.open(io.istor, me.fstor, status='UNKNOWN')
    io.rewind(io.istor)
    if nml.lstore:
        # JvS 08/25: stitched interior and envelope option.
        if nml.lstch:
            ccout_mod.lphhd = True
        else:
            io.write(io.istor, verline)
            io.write(io.istor, _STORE_KEY)       # format 1012

    io.open(io.ishort, me.fshort, status='UNKNOWN')
    io.rewind(io.ishort)
    # The sound-speed profile, opened always (LSOUND is known only later).
    io.open(io.isound, me.fsound, status='UNKNOWN')
    io.rewind(io.isound)
    if nml.ltrack:
        io.open(io.itrack, me.ftrack, status='UNKNOWN')
        io.rewind(io.itrack)
    if nml.lsnu:
        io.open(io.isnu, me.fsnu, status='UNKNOWN')
        io.rewind(io.isnu)
    if nml.ldebug:
        io.open(io.idebug, me.fdebug, status='UNKNOWN')
    if nml.lmilne:
        io.open(io.imilne, me.fmilne, status='UNKNOWN', access='APPEND')
    io.open(io.ifirst, me.ffirst, status='OLD')
    io.open(io.ilast, me.flast, status='UNKNOWN')
    io.open(io.imodpt, me.fmodpt, status='UNKNOWN')
    if nml.lpulse:
        io.open(io.iopmod, args.fpmod, status='UNKNOWN')
        io.open(io.iopenv, args.fpenv, status='UNKNOWN')
        io.open(io.iopatm, args.fpatm, status='UNKNOWN')
    if nml.lcondopacp:
        io.open(io.icondopacp, me.fcondopacp, status='OLD')

    # Semiconvection with any overshoot stops the run.
    if nml.lsemic:
        if nml.lovstc or nml.lovste or nml.lovstm:
            io.write(io.ishort,
                     f" ERROR IN SUBROUTINE PARMIN\nSEMI-CONVECTION AND "
                     f"OVERSHOOT FLAGS BOTH TURNED ON\nFLAGS LSEMIC "
                     f"{_l(nml.lsemic)} OVERSHOOT - CORE,ENVELOPE,"
                     f"INTERMEDIATE- {_l(nml.lovste)} {_l(nml.lovstc)} "
                     f"{_l(nml.lovstm)}\nRUN STOPPED")          # format 2
            raise SystemExit
    io.write(io.ishort,
             " PT TOL" + ''.join(f"{nml.hpttol[i]:6.3f}" for i in range(1, 13))
             + f"\n O.S.ENV{nml.alphae:6.3f} O.S.CORE{nml.alphac:6.3f}"
             f" LINSTB {_l(nml.linstb)} LJDOT {_l(nml.ljdot0)} WIND IND."
             f"{me.alfa:6.3f} FK{me.fk:8.2E}\n  FV{nml.fw:5.2f} FC"
             f"{nml.fc:5.2f} COUPLING{byp.fo:6.3f} F MU{nml.fmu:5.2f}"
             f" RCRIT{nml.rcrit:9.1f}")                        # format 1

    if nml.lrot:
        nml.lnew0 = True
        if nml.walpcz < -2.0:
            nml.walpcz = -2.0
        if nml.walpcz > 0.0:
            nml.walpcz = 0.0
        # JNT 09/2025: LSOLID overrides IMPJMOD.
        if nml.lsolid:
            nml.impjmod = 1
        if not nml.lmwind:
            # The old wind law:
            # JDOT = FK*2.036D33*1.452D9**ALFA*(MDOT/1.0D-14)**(1-2ALFA/3)
            #        *OMEGA**(1+4ALFA/3)*(R/RSUN)**(2-ALFA)*(M/MSUN)**-ALFA/3
            c.cc13 = 1.0 / 3.0
            c.cc23 = 2.0 / 3.0
            nml.constfactor = me.fk * 2.036e33 * 1.452e9**me.alfa
            nml.exmd = 1.0 - c.cc23 * me.alfa
            nml.exw = 1.0 + 2.0 * c.cc23 * me.alfa
            nml.exr = 2.0 - me.alfa
            nml.exm = -c.cc13 * me.alfa
        else:
            # The Matt et al. law: the indices a, b, c, d, m, preset by
            # AWIND ('CUS' keeps the namelist's), give the exponents.
            if me.awind == 'K97':
                byp.pmmm = me.alfa / 3.0
                # MHP 8/17: the default for A is 2, not 0.
                byp.pmma = 2.0
                byp.pmmb = 1.0
                byp.pmmc = 2.0
                byp.pmmd = 0.0
                nml.lbscale = False
            elif me.awind == 'V13':
                byp.pmmm = 0.22
                byp.pmma = 2.0
                byp.pmmb = 1.0
                byp.pmmc = 0.0
                byp.pmmd = 1.0
                nml.lbscale = True
            nml.exw = (1.0 + byp.pmma - 2.0 * byp.pmma * byp.pmmm
                       + 4.0 * byp.pmmm * byp.pmmb)
            # G Somers 8/17: the tau_cz term is now in MWIND/MCOWIND.
            nml.extau = 0.0
            # JvS 09/25 "removed typo". Question (see the page): the
            # comment's derivation gives another exponent.
            nml.exr = 2.0 + 5.0 * byp.pmmm - 4.0 * byp.pmmm * byp.pmmc
            nml.exm = -byp.pmmm
            nml.exl = 1.0 - 2.0 * byp.pmmm
            nml.expr = 2.0 * byp.pmmm * byp.pmmd
            nml.excen = byp.pmmm * byp.pmmd
            # The constant factor of the centrifugal term.
            c.c_2 = 0.0506
            nml.constfactor = me.fk * byp.pmmjd / nml.pmmsolw**nml.exw
            # WMAX_SUN < 1000: saturation, at WMAX_SUN times the solar rate.
            if me.wmax_sun < 1.0e3:
                nml.wmax = me.wmax_sun * nml.pmmsolw
    # Wind law end.
    c.cln = math.log(10.0)

    # A single species to rescale in the starting model.
    if nml.lnewcp:
        nml.lrel = True
        if atmp.rstrip() == 'ABS':
            nml.lrel = False
        for i in range(1, 13):
            if anewcp.rstrip() == aid[i]:
                # INEWCP, the index of the element in HCOMP.
                nml.inewcp = i + 3
                break
        else:
            # ANEWCP is not a recognized element (format 20).
            nml.lnewcp = False
            io.write(io.ishort, f" VARIABLE{anewcp:>4} NOT A RECOGNIZED "
                     f"ELEMENT\n RESCALING NOT PERFORMED")
        # label 30

    # MHP 10/24: the CNO mixture (ISETMIX) and the isotope ratios and
    # light-element abundances (ISETISO) of the starting model.
    nml.lmixture = False
    nml.lisotope = False
    to_606 = True        # False where the Fortran jumps to label 602
    if byp.isetmix == 1:
        if me.amix.rstrip() == 'CUS':
            # A custom mixture: no negative fraction, sum below one.
            if nml.frac_c < 0.0 or nml.frac_n < 0.0 or nml.frac_o < 0.0:
                _both(f"NEGATIVE INPUT CNO FRACTION {nml.frac_c:12.4E}"
                      f"{nml.frac_n:12.4E}{nml.frac_o:12.4E}"
                      " MIX NOT MODIFIED")                       # 591
                to_606 = False
            else:
                sum_frac = nml.frac_c + nml.frac_n + nml.frac_o
                if sum_frac >= 1.0:
                    _both(f"INPUT CNO FRACTION {nml.frac_c:12.4E}"
                          f"{nml.frac_n:12.4E}{nml.frac_o:12.4E}"
                          " EXCEEDS 1. MIX NOT MODIFIED")        # 598
                    to_606 = False
                else:
                    nml.lmixture = True
        else:
            # A preset mixture, if its name is known.
            for i in range(1, 5):
                if me.amix.rstrip() == amixt[i]:
                    byp.zxmix = zxmixt[i]
                    nml.frac_c = frac_ct[i]
                    nml.frac_n = frac_nt[i]
                    nml.frac_o = frac_ot[i]
                    nml.lmixture = True
                    break
            else:
                _both(f"DESIRED CNO MIXTURE {me.amix:8} NOT FOUND."
                      " MIX NOT ALTERED.")                       # 589
    # label 606
    if to_606 and nml.lmixture:
        _both(f"CNO MIXTURE {me.amix:8} C {nml.frac_c:12.4E} N "
              f"{nml.frac_n:12.4E} O {nml.frac_o:12.4E} APPLIED TO "
              "STARTING MODEL.")                                 # 604
    # label 602: the isotope ratios.
    to_605 = True        # False where the Fortran jumps to label 603
    if byp.isetiso == 1:
        if me.aiso.rstrip() == 'CUS':
            ratios = (nml.r12_13, byp.r14_15, byp.r16_17, nml.r16_18)
            light = (nml.xh2_ini, nml.xhe3_ini, nml.xli6_ini, nml.xli7_ini,
                     nml.xbe9_ini, byp.xb10_ini, byp.xb11_ini)
            if any(v < 0.0 for v in ratios + light):
                _both("NEGATIVE INPUT ISOTOPE RATIO OR LIGHT ELEMENT MASS "
                      "FRACTION " + ''.join(f"{v:12.4E}"
                                            for v in ratios + light)
                      + " MIX NOT MODIFIED")                     # 596
                to_605 = False
            else:
                sum_frac = (nml.xh2_ini + nml.xhe3_ini + nml.xli6_ini
                            + nml.xli7_ini + nml.xbe9_ini + byp.xb10_ini
                            + byp.xb11_ini)
                if sum_frac >= 1.0:
                    _both("SUM OF LIGHT ELEMENT MASS FRACTIONS EXCEEDS 1"
                          + ''.join(f"{v:12.4E}" for v in light)
                          + " MIX NOT MODIFIED")                 # 595
                    to_605 = False
        elif me.aiso.rstrip() != 'L21':
            # Only the default (L21) or a custom set are valid.
            to_605 = False
        if to_605:
            nml.lisotope = True
    if to_605 and nml.lisotope:
        _both(f"ISOTOPE AND LIGHT ELEMENT MIXTURE {me.aiso:8} C12/C13 "
              f"{nml.r12_13:12.4E} O16/O18 {nml.r16_18:12.4E} H2 "
              f"{nml.xh2_ini:12.4E} HE3 {nml.xhe3_ini:12.4E} LI6 "
              f"{nml.xli6_ini:12.4E} LI7 {nml.xli7_ini:12.4E} BE9 "
              f"{nml.xbe9_ini:12.4E} APPLIED TO STARTING MODEL.")  # 605
    # label 603: the settling switches. DBG 12/95: Z diffusion needs Y
    # diffusion and Thoul's coefficients.
    if nml.ldifz:
        nml.ldify = True
        nml.lthoul = True
    # G Somers 04/15: light-element diffusion.
    if nml.ldifli:
        nml.lnewdif = True
        nml.ldify = True
        nml.ldifz = True
        nml.lthoul = True
        nml.lthoulfit = False
        nml.ilambda = 4

    # The run's parameters: the namelists as read, and the derived values.
    if nml.npoint <= 0:
        nml.npoint = 9999
    if byp.npenv <= 0:
        byp.npenv = 9999
    if nml.nprtmod <= 0:
        nml.nprtmod = 9999
    if nml.nprtpt <= 0:
        nml.nprtpt = 9999
    _write_namelist(io.ishort, 'PHYSICS', _PHYSICS, owners)
    _write_namelist(io.ishort, 'CONTROL', _CONTROL, owners)
    _write_namelist(io.ishort, 'DERIVED', _DERIVED, owners)

    # 10/15/26: every mixing length the run uses must be positive. Under
    # solar calibration only CMIXLA(1) comes from the namelist.
    if nml.lcals:
        ncmxck = 1
    else:
        ncmxck = min(nml.numrun, 50)
    lcmxbad = False
    for i in range(1, ncmxck + 1):
        if nml.cmixla[i] <= 0.0:
            lcmxbad = True
            text = (f" ERROR IN SUBROUTINE PARMIN\nILLEGAL MIXING LENGTH: "
                    f"CMIXLA({i:2d}) ={nml.cmixla[i]:12.4E} - MUST BE "
                    "POSITIVE")                                  # 8801
            io.write(io.ishort, text)
            print(text)
    if lcmxbad:
        if nml.lcals:
            text = ("SOLAR CALIBRATION IS ON, SO ONLY CMIXLA(1) IS READ FROM"
                    " THE NAMELIST\nRUN STOPPED")                # 8802
        else:
            text = (f"CMIXLA MUST BE SET FOR EACH OF THE{ncmxck:3d} KIND "
                    "CARDS\nRUN STOPPED")                        # 8803
        io.write(io.ishort, text)
        print(text)
        raise SystemExit

    if nml.kttau in _KTTAU_TEXT:
        io.write(io.ishort, _KTTAU_TEXT[nml.kttau])
    if nml.lpulse:
        io.write(io.ishort, "\n CALCULATE PULSATION OUTPUT ON LAST MODEL")
    if nml.lpurez:
        io.write(io.ishort, '  USING PURE C AND N OPACITY TABLES')

    # The run's description (formats 310, 314, 320).
    rule = '#' + '=' * 100
    description = (f"# DESCRIPTION OF RUN:{descrip[1]:<80.80}\n#"
                   + ' ' * 9 + '  ' + ' ' * 8
                   + f": {descrip[2]:<80.80}\n{rule}")
    io.write(io.imodpt, description)
    io.write(io.ishort, "#\n" + rule)
    io.write(io.ishort, verline)
    io.write(io.ishort, description)
    if nml.ltrack:
        io.write(io.itrack, "#\n" + rule)
        io.write(io.itrack, verline)
        io.write(io.itrack, description)
    io.write(io.ishort, " USING OSCILATORY SPLINE INTERPOLATION IN HPOINT")

    # Interpret the run from the sequence of "kind" cards.
    io.write(io.ishort, "\n" + ' ' * 35 + "RUN CARDS\n")
    nml.lfirst[1] = True
    rl = runloop_mod
    for nkind in range(1, nml.numrun + 1):
        rl.iresca[nkind] = kindrn[nkind]
        if kindrn[nkind] == 1:
            # Evolve card. MHP 10/24: stop at an age, or at a central D, X
            # or Y.
            if (nml.endage[nkind] > 0.0 or nml.end_dcen[nkind] > 0.0
                    or nml.end_xcen[nkind] > 0.0
                    or nml.end_ycen[nkind] > 0.0):
                rl.lendag[nkind] = True
            else:
                rl.lendag[nkind] = False
            rl.lsetdt[nkind] = nml.setdt[nkind] > 0.0
            if nml.nmodls[nkind] > 0:
                if nml.lfirst[nkind]:
                    _screen_and_short(
                        f"\n RUN #{nkind:3d}   EVOLVE {nml.nmodls[nkind]:5d}"
                        ' MODELS, STARTING WITH THE INPUT "FIRST MODEL".')
                else:
                    _screen_and_short(
                        f"\n RUN #{nkind:3d}   EVOLVE {nml.nmodls[nkind]:5d}"
                        " MODELS, STARTING WITH THE PREVIOUS RUN'S LAST "
                        "MODEL.")
                if rl.lendag[nkind] or rl.lsetdt[nkind]:
                    _screen_and_short(
                        f" EVOLVE TO AGE {_l(rl.lendag[nkind])} SET DELT "
                        f"{_l(rl.lsetdt[nkind])} FINAL AGE "
                        f"{nml.endage[nkind]:9.2E} FIXED TSTEP "
                        f"{nml.setdt[nkind]:9.2E} CENTRAL D "
                        f"{nml.end_dcen[nkind]:10.4E} CENTRAL X "
                        f"{nml.end_xcen[nkind]:12.4E} CENTRAL Y "
                        f"{nml.end_ycen[nkind]:12.4E}")          # 370
        elif kindrn[nkind] == 2:
            # Rescale card: mass, X, Z and core mass, in RESCAL(1..4, card).
            rl.rescal[1, nkind] = rsclm[nkind]
            rl.rescal[2, nkind] = rsclx[nkind]
            rl.rescal[3, nkind] = rsclz[nkind]
            rl.rescal[4, nkind] = rsclcm[nkind]
            if nml.nmodls[nkind] > 0:
                if nml.lfirst[nkind]:
                    _screen_and_short(f"\n RUN #{nkind:3d}   RESCALE THE "
                                      'INPUT MODEL: "FIRST MODEL".')  # 450
                else:
                    _screen_and_short(f"\n RUN #{nkind:3d}   RESCALE THE "
                                      "PREVIOUS RUN'S LAST MODEL.")   # 451
                _screen_and_short(_relax_text(nml.nmodls[nkind], nkind))
        elif kindrn[nkind] == 3:
            # Rescale-and-evolve card.
            rl.rescal[1, nkind] = rsclm[nkind]
            rl.rescal[2, nkind] = rsclx[nkind]
            rl.rescal[3, nkind] = rsclz[nkind]
            rl.rescal[4, nkind] = rsclcm[nkind]
            if nml.lfirst[nkind]:
                _screen_and_short(f"\n RUN #{nkind:3d}   RESCALE & EVOLVE "
                                  'THE INPUT MODEL: "FIRST MODEL".')  # 550
            else:
                # Question (see the page): format 451, its own 551 is
                # commented out.
                _screen_and_short(f"\n RUN #{nkind:3d}   RESCALE THE "
                                  "PREVIOUS RUN'S LAST MODEL.")
            _screen_and_short(_relax_text(nml.nmodls[nkind], nkind))
        # Question (see the page): runs for every card, evolve cards too.
        if rl.rescal[3, nkind] >= 0.0:
            comp_mod.zenv = rl.rescal[3, nkind]
    # label 1000
    a = args
    return (a.falex06, a.fallard, a.fatm, a.ffermi, a.fliv95, a.fpatm,
            a.fpenv, a.fpmod, a.fpurez, a.fscvh, a.fscvhe, a.fscvz)


def _relax_text(nmodls, nkind):
    """Format 452: the rescaling of a rescale card."""
    r = runloop_mod.rescal
    return (f" RELAX RESCALED MODEL{nmodls:3d} TIMES. RESCALE THE FOLLOWING"
            f"(0=USE CURRENT VALUE):\n MASS {r[1, nkind]:9.6f}   X"
            f"{r[2, nkind]:9.6f}   Z{r[3, nkind]:9.6f}   CORE MASS"
            f"{r[4, nkind]:9.6f}")


# The placeholders EXPAND_VALUE knows, with the defaults used when the
# environment variable of that name is not set.
_PHOLDERS = ('{YREC_INPUT}', '{YREC_START}', '{YREC_OUTPUT}')
_DEFAULTS = ('../../input', '../../startmodels', 'output')


def expand_value(nml_var):
    """Replace a placeholder at the start of a file name, {YREC_INPUT},
    {YREC_START} or {YREC_OUTPUT}, by the environment variable of that
    name, or by its default if it is not set.

    Translates: parmin.f#EXPAND_VALUE (09c3b9b)

    Arguments: nml_var, the file name.

    Returns nml_var, expanded.
    """
    for i in range(len(_PHOLDERS)):
        pholder = _PHOLDERS[i]
        ph_len = len(pholder)
        if nml_var[:ph_len] == pholder:
            temp = os.environ.get(pholder[1:ph_len - 1], '').rstrip()
            if len(temp) != 0:
                # The environment variable is set.
                temp = temp + nml_var.rstrip()[ph_len:]
            else:
                # Not set: use the default path.
                temp = _DEFAULTS[i] + nml_var.rstrip()[ph_len:]
            nml_var = temp
    # label 5000
    return nml_var
