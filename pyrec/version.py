"""Python version of version.f (YREC 09c3b9b): the YREC version and the git commit it was built from.

The Fortran takes both from preprocessor macros: YREC_VERSION is written in
the file, GIT_HASH is passed by the build (-DGIT_HASH, the short hash with
"-dirty" appended, or "unknown"), defaulting to 'unknown'.
"""

YREC_VERSION = '202602'
GIT_HASH = 'unknown'     # the file's default; the build passes the real hash


def setversion(yrecver, githash):
    """Set the version string and the git hash (both only set).

    Translates: version.f#SETVERSION (09c3b9b)

    Arguments: yrecver (CHARACTER*10), githash (CHARACTER*20), any value.
    Returns (yrecver, githash), padded with blanks to 10 and 20 characters
    as Fortran's fixed-length assignment does.
    """
    yrecver = YREC_VERSION.ljust(10)[:10]
    githash = GIT_HASH.ljust(20)[:20]
    return yrecver, githash
