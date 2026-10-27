"""Python version of wrtlst.f (YREC 09c3b9b): write a converged model in
MODEL2 format."""

from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec.putmodel2 import putmodel2

# ATM is SAVEd: it keeps its last value when KTTAU is outside 0-5.
_saved = {'atm': None}

_ATM_LABELS = {0: 'EDD ', 1: 'KS  ', 2: 'HRA ', 3: 'KUR ', 4: 'ALL ',
               5: 'K/C '}  # 5: JNT 06/14


def wrtlst(iwrite, hcomp, hd, hl, hp, hr, hs, ht, lc, trit, tril, ps,
           ts, rs, cfenv, ftri, tlumx, jcore, jenv, model, m, smass, teffl,
           bl, hstot, dage, ddage, omega, cmixl, compmix):
    """Set the header labels for the atmosphere, EOS and opacities, and
    write the model to unit iwrite in MODEL2 format through putmodel2.

    Translates: wrtlst.f#WRTLST (09c3b9b)

    Arguments: the unit iwrite; the model (arrays hcomp ... omega, the
    envelope triangle trit, tril, ps, ts, rs, cfenv, ftri, tlumx) and its
    scalars; the mixing length cmixl; compmix, the initial composition's
    code (a string). Nothing is returned.

    Units are Fortran unit numbers; writes to the unit numbers iowr and
    ishort are shown as print().
    """
    n = nmlpar_mod
    # Physics flags. Atmosphere flag ATM.
    # Question (see the page): KTTAU outside 0-5 leaves ATM unchanged.
    if n.kttau in _ATM_LABELS:
        _saved['atm'] = _ATM_LABELS[n.kttau]
    atm = _saved['atm']

    # Equation of state flag EOS.
    eos = 'SAHA  '
    if n.ldh:
        eos = 'SAH+DH'
    if n.lscv:
        eos = 'SCV   '
        if n.ldh:
            # Question (see the page): 'SCV+O6' is overwritten at once.
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

    # Low-temperature opacity flag ALOK (no options left since 09c3b9b).
    alok = 'NONE'
    # High-temperature opacity flag HIK.
    hik = 'NONE'
    if n.lopal95:
        hik = 'OP95'

    # The first three lines are YREC7 inputs, the last two MODEL2 add-ons.
    putmodel2(bl, cfenv, cmixl, dage, ddage, ftri, hcomp, hd, hl,
              hp, hr, hs, hstot, ht, iwrite, iounits_mod.ishort, jcore, jenv,
              lc, n.lexcom, n.lrot, m,
              model, omega, ps, rs, smass, teffl, tlumx, tril, trit, ts,
              atm, eos, hik, n.ldify, n.ldifz, n.ldisk, n.linstb, n.ljdot0,
              alok, n.lovstc, n.lovste, n.lovstm, n.lpurez, n.lsemic,
              compmix, n.pdisk, n.tdisk, n.wmax)

    # To unit IOWR (format 360).
    print(f"{model:6d}  #SHELLS={m:4d}  LogTeff={teffl:8.5f}"
          f"  Log(L/Lsun)={bl:8.5f}  Age={dage:12.5f}")
    # To unit ISHORT.
    # Question (see the page): the hard-coded unit 11.
    if iwrite == 11:
        print(f" DUMPED MODEL{model:5d}  FILE{iwrite:3d}")          # 330
    else:
        print(f" DUMPED MODEL{model:5d} AGE{dage:13.9f}  FILE{iwrite:3d}")
    if iwrite == iounits_mod.ilast:
        pass  # REWIND ILAST: the unit holds only the latest model.
