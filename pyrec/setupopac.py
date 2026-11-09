"""Python version of setupopac.f (YREC 7be4133): read the opacity tables
chosen for the run."""

from pyrec import newopac_mod
from pyrec import nmlpar_mod
from pyrec.ll95tbl import ll95tbl
from pyrec.op95xtab import op95xtab
from pyrec.rdzlaol import rdzlaol
from pyrec.readalex06 import readalex06
from pyrec.zsulaol import zsulaol


def setupopac(xenv: float, falex06: str, fliv95: str, fpurez: str) -> None:
    """Read the interior (OPAL95, and the pure-Z table) and low-temperature
    (Ferguson et al. 2005) opacity tables the namelist switches on, and set
    L2Z, which tells the readers to read a second table at another Z.

    Translates: setupopac.f#SETUPOPAC (7be4133)

    Arguments: xenv, the envelope's X; the file names falex06, fliv95,
    fpurez. None is set.

    """

    n = nmlpar_mod

    # With LZRAMP or LDIFZ, read a second set of tables at another Z (Fortran line 64).
    newopac_mod.l2z = n.lzramp or n.ldifz

    # Interior tables: OPAL95 (Fortran lines 69–72).
    if n.lopal95:
        ll95tbl(fliv95)
        op95xtab(xenv)

    # The pure-Z table (Fortran lines 79–82).
    if n.lpurez:
        rdzlaol(fpurez)
        zsulaol()

    # Low-temperature tables: Alexander/Ferguson 2006 (Fortran lines 87–93).
    if n.lalex06:
        readalex06(falex06)
