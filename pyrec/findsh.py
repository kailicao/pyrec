"""Python version of findsh.f (YREC 09c3b9b): the edges of the
convective core, the surface convection zone and the hydrogen shell."""

from pyrec import nmlpar_mod

# DATA: the shell ends where L grows by less than dlmin * L_surf, or X is
# within dxmin of its surface value.
dlmin = 1.0e-5
dxmin = 1.0e-5


def findsh(hcomp, hl, lc, m, jcore, jenv, jxbeg, jxend, jxmid, lshell):
    """Find the outermost point of the convective core, the innermost
    point of the surface convection zone and, once the central X is at
    most nmlpar_mod.atime[1], the start, middle and end of the hydrogen
    shell.

    Translates: findsh.f#FINDSH (09c3b9b)

    Arguments: hcomp (abundances, row 1 is X), hl (luminosities), lc
    (convective flags), m (number of points); jcore, jenv, jxbeg, jxend,
    jxmid, lshell are only set.

    Returns (jcore, jenv, jxbeg, jxend, jxmid, lshell).
    """
    jxbeg = 1
    jxmid = 1
    jxend = 1
    lshell = False

    if hcomp[1, 1] <= nmlpar_mod.atime[1]:
        lshell = True
        xmid = 0.50 * hcomp[1, m]
        xlend = dlmin * hl[m]
        # Question (see the page): jxmid counts points rather than
        # recording an index.
        i = m  # if no end is found
        for ii in range(1, m + 1):  # label 10
            if hcomp[1, ii] <= nmlpar_mod.atime[1]:
                jxbeg = jxbeg + 1
                jxmid = jxmid + 1
            elif hcomp[1, ii] <= xmid:
                jxmid = jxmid + 1
            elif hl[ii] - hl[ii - 1] < xlend:
                i = ii
                break  # GOTO 20
            elif hcomp[1, m] - hcomp[1, ii] < dxmin:
                i = ii
                break  # GOTO 20
        jxend = i  # label 20

    # The core: the last point of the convective run from the centre.
    i = m + 1  # Fortran's index after a loop that runs out
    for ii in range(1, m + 1):  # label 30
        if not lc[ii]:
            i = ii
            break  # GOTO 40
    if i > 1:  # label 40
        jcore = i - 1
    else:
        jcore = 1

    # The envelope: the first point of the convective run from the surface.
    i = 0
    for ii in range(m, 0, -1):  # label 50
        if not lc[ii]:
            i = ii
            break  # GOTO 60
    if i < m:  # label 60
        jenv = i + 1
    else:
        jenv = m

    # Fully convective: JCORE is used in the nuclear time step, and
    # putting it at the surface does strange things.
    if jenv == 1:
        jcore = 1
    return jcore, jenv, jxbeg, jxend, jxmid, lshell
