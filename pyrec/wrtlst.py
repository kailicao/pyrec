"""Python version of wrtlst.f (YREC 7be4133): write a converged model in
MODEL2 format."""

import numpy as np

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.putmodel2 import putmodel2

# ATM is SAVEd: it keeps its last value when KTTAU is outside 0-5.
_saved = {'atm': None}

_ATM_LABELS = {0: 'EDD ', 1: 'KS  ', 2: 'HRA ', 3: 'KUR ', 4: 'ALL ',
               5: 'K/C '}  # 5: Kurucz/Castelli (JNT 06/14)


def wrtlst(iwrite: int, hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray,
           hs: np.ndarray, ht: np.ndarray, lc: np.ndarray, trit: np.ndarray, tril: np.ndarray, ps: np.ndarray,
           ts: np.ndarray, rs: np.ndarray, cfenv: np.ndarray, ftri: float, tlumx: np.ndarray, jcore: int, jenv: int,
           model: int, m: int, smass: float, teffl: float, bl: float, hstot: float, dage: float, ddage: float,
           omega: np.ndarray, cmixl: float, compmix: str) -> None:
    """Set the header labels for the atmosphere, EOS and opacities, and
    write the model to unit iwrite in MODEL2 format through putmodel2.

    Translates: wrtlst.f#WRTLST (7be4133)

    Arguments: the unit iwrite; the model (arrays hcomp ... omega, the
    envelope triangle trit, tril, ps, ts, rs, cfenv, ftri, tlumx) and its
    scalars; the mixing length cmixl; compmix, the initial composition's
    code (a string). Nothing is returned.

    Units are Fortran unit numbers; writes to the unit numbers iowr and
    ishort are shown as print().

    """

    n = nmlpar_mod

    # Physics flags. Atmosphere flag ATM (Fortran lines 91–106).
    # Suspected bug (see the page): KTTAU outside 0-5 leaves ATM unchanged.
    if n.kttau in _ATM_LABELS:
        _saved['atm'] = _ATM_LABELS[n.kttau]
    atm = _saved['atm']

    # Equation of state flag EOS (Fortran lines 107–124).
    eos = 'SAHA  '
    if n.ldh:
        eos = 'SAH+DH'
    if n.lscv:
        eos = 'SCV   '
        if n.ldh:
            # Suspected bug (see the page): 'SCV+O6' is overwritten at once.
            if n.lopale06:
                eos = 'SCV+O6'
            eos = 'SCV+DH'
            if n.lopale06:
                eos = 'SCDHO6'
    else:
        if n.lopale06:
            eos = 'OPAL06'
            if n.ldh:
                eos = 'OP6+DH'

    # Low-temperature opacity flag ALOK (no options left) (Fortran lines
    # 125–130).
    alok = 'NONE'

    # High-temperature opacity flag HIK (Fortran lines 131–135).
    hik = 'NONE'
    if n.lopal95:
        hik = 'OP95'

    # The first three lines are YREC7 inputs, the last two MODEL2 add-ons
    # (Fortran lines 137–143).
    putmodel2(bl, cfenv, cmixl, dage, ddage, ftri, hcomp, hd, hl,
              hp, hr, hs, hstot, ht, iwrite, iounits_mod.ishort, jcore, jenv,
              lc, n.lexcom, n.lrot, m,
              model, omega, ps, rs, smass, teffl, tlumx, tril, trit, ts,
              atm, eos, hik, n.ldify, n.ldifz, n.ldisk, n.linstb, n.ljdot0,
              alok, n.lovstc, n.lovste, n.lovstm, n.lpurez, n.lsemic,
              compmix, n.pdisk, n.tdisk, n.wmax)

    # To unit IOWR (format 360) (Fortran lines 147–149).
    print(f"{model:6d}  #SHELLS={m:4d}  LogTeff={teffl:8.5f}"
          f"  Log(L/Lsun)={bl:8.5f}  Age={dage:12.5f}")

    # To unit ISHORT; unit 11 is ILAST (Fortran lines 150–156).
    if iwrite == 11:
        print(f" DUMPED MODEL{model:5d}  FILE{iwrite:3d}")          # 330
    else:
        print(f" DUMPED MODEL{model:5d} AGE{dage:13.9f}  FILE{iwrite:3d}")

    # Fortran line 157.
    if iwrite == iounits_mod.ilast:
        pass  # REWIND ILAST: the unit holds only the latest model.
