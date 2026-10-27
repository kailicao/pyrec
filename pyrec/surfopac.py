"""Python version of surfopac.f (YREC development, 09c3b9b)."""

from pyrec import nmlpar_mod
from pyrec.op95xtab import op95xtab


def surfopac(xenv):
    """Set up the opacity tables for the surface hydrogen abundance.

    Assumes the tables have already been read. Only the OPAL 1995
    tables are interpolated in X here.

    Translates: surfopac.f#SURFOPAC (09c3b9b)

    Parameters
    ----------
    xenv : float
        Surface (envelope) hydrogen mass fraction X.
    """
    # Interior tables: set up the OPAL95 tables.
    if nmlpar_mod.lopal95:
        op95xtab(xenv)
