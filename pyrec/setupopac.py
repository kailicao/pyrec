"""Python version of setupopac.f (YREC 09c3b9b): read the opacity tables
chosen for the run."""

from pyrec import newopac_mod
from pyrec import nmlpar_mod
from pyrec.ll95tbl import ll95tbl
from pyrec.op95xtab import op95xtab
from pyrec.rdzlaol import rdzlaol
from pyrec.readalex06 import readalex06
from pyrec.zsulaol import zsulaol


def setupopac(xenv, falex06, fliv95, fpurez):
    """Read the interior (OPAL95, and the pure-Z table) and low-temperature
    (Ferguson et al. 2005) opacity tables the namelist switches on, and set
    L2Z, which tells the readers to read a second table at another Z.

    Translates: setupopac.f#SETUPOPAC (09c3b9b)

    Arguments: xenv, the envelope's X; the file names falex06, fliv95,
    fpurez. None is set.
    """
    n = nmlpar_mod
    # With LZRAMP or LDIFZ, read a second set of tables at another Z.
    newopac_mod.l2z = n.lzramp or n.ldifz
    # Interior tables: OPAL95.
    if n.lopal95:
        ll95tbl(fliv95)
        op95xtab(xenv)
    # The pure-Z table.
    if n.lpurez:
        rdzlaol(fpurez)
        zsulaol()
    # Low-temperature tables: Alexander/Ferguson 2006.
    if n.lalex06:
        readalex06(falex06)
