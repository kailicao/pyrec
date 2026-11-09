"""Python version of surfopac.f (YREC development, 7be4133)."""

from pyrec import nmlpar_mod
from pyrec.op95xtab import op95xtab


def surfopac(xenv: float) -> None:
    """Set up the opacity tables for the surface hydrogen abundance.

    Assumes the tables have already been read. Only the OPAL 1995
    tables are interpolated in X here.

    Translates: surfopac.f#SURFOPAC (7be4133)

    Parameters
    ----------
    xenv : float
        Surface (envelope) hydrogen mass fraction X.

    """
    # Interior tables: set up the OPAL95 tables (Fortran lines 36–41).
    if nmlpar_mod.lopal95:
        op95xtab(xenv)
