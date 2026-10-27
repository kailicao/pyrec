"""Python version of pdist.f (YREC 09c3b9b): when to write a model for
pulsation calculations."""

from pyrec import iounits_mod
from pyrec import nmlpar_mod


def _fmt276(polen, d1, d2, d3):
    # FORMAT 276: 5X,' LEN SQ=',1PE10.3,'  D L=',1PE10.3,'  D Te=',1PE10.3,'  D AGE=',1PE10.3
    return (" " * 5 + f" LEN SQ={polen:10.3E}  D L={d1:10.3E}"
            f"  D Te={d2:10.3E}  D AGE={d3:10.3E}")


def _pulse_name(base, nk, modeln):
    """The file name <base>NN_MMMM (or NN_MMMMM from model 10000 on).

    base is a CHARACTER*256 in Fortran, padded with blanks: the name ends
    at its first blank, INDEX(base, ' ') - 1. A Python string has no
    padding, so a name without blanks is taken whole."""
    ioccol = base.find(" ") if " " in base else len(base)
    if modeln < 10000:
        ctemp = f"{nk:02d}_{modeln:04d}"
    else:
        ctemp = f"{nk:02d}_{modeln:05d}"
    return base[:ioccol] + ctemp


def pdist(pol1, pot1, poa1, polen, bl, teffl, modeln, nk, fpatm, fpenv,
          fpmod, dage):
    """Accumulate the squared distance the star has moved in the HR diagram
    and in age since the last pulsation model; when it reaches POMAX, open
    three new pulsation output files and set LPULSE, otherwise close them
    and clear LPULSE.

    Translates: pdist.f#PDIST (09c3b9b)

    Arguments: pol1, pot1, poa1 (log L, log Teff and age at the previous
    model), polen (the accumulated squared distance); bl, teffl, the
    model number modeln, the run index nk and the age dage now; the base
    file names fpatm, fpenv, fpmod.

    Returns (pol1, pot1, poa1, polen), updated. Sets nmlpar_mod.lpulse.
    """
    tm1 = bl - pol1
    tm2 = teffl - pot1
    tm3 = dage - poa1
    tmp1 = nmlpar_mod.poa * tm1 * nmlpar_mod.poa * tm1
    tmp2 = nmlpar_mod.pob * tm2 * nmlpar_mod.pob * tm2
    tmp3 = nmlpar_mod.poc * tm3 * nmlpar_mod.poc * tm3
    # Question (see the page): POMAX is compared with a squared distance.
    polen = polen + tmp1 + tmp2 + tmp3
    pol1 = bl
    pot1 = teffl
    poa1 = dage
    if polen >= nmlpar_mod.pomax:
        nmlpar_mod.lpulse = True
        polen = 0.0
        # Question (see the page): STATUS='NEW' stops if the file exists.
        iounits_mod.open(iounits_mod.iopmod, _pulse_name(fpmod, nk, modeln),
                         status="new")
        iounits_mod.open(iounits_mod.iopenv, _pulse_name(fpenv, nk, modeln),
                         status="new")
        iounits_mod.open(iounits_mod.iopatm, _pulse_name(fpatm, nk, modeln),
                         status="new")
    else:
        iounits_mod.close(iounits_mod.iopmod)
        iounits_mod.close(iounits_mod.iopenv)
        iounits_mod.close(iounits_mod.iopatm)
        nmlpar_mod.lpulse = False
    # Question (see the page): the differences go to IOWR, their weighted
    # squares to ISHORT, under the same labels.
    iounits_mod.write(iounits_mod.iowr, _fmt276(polen, tm1, tm2, tm3))
    iounits_mod.write(iounits_mod.ishort, _fmt276(polen, tmp1, tmp2, tmp3))
    return pol1, pot1, poa1, polen
