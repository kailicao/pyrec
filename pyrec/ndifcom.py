"""Python version of ndifcom.f (YREC 09c3b9b): composition mixed across
each region made unstable by rotation, through MIXGRID and MIXCOM."""

from pyrec.mixcom import mixcom
from pyrec.mixgrid import mixgrid

# The Fortran's locals (SAVE).
ldum = False
ldum2 = False


def ndifcom(dt, cod2, ecod2, em, hd, hl, hp, hr, hs, hs1, hs2, hstot, hv,
            ibeg, iend, imax, imin, lcz, lok, m, hcomp, jbeg, jend):
    """Before the final iteration (lok False), mix rows jbeg..jend over the
    region found earlier (ibeg..iend, with its grid ecod2, em). On the
    final iteration, find each run of interfaces with hv > 0 between imin
    and imax, set up its grid (MIXGRID) and mix it (MIXCOM); regions of a
    single unstable interface are skipped.

    Translates: ndifcom.f#NDIFCOM (09c3b9b)

    Arguments: the step dt (s); the coefficients cod2 at the model points;
    ecod2 and em on the equal grid (kept between calls by the caller); the
    model hd, hl, hp, hr, hs, hs1, hs2, hstot; the velocities hv; ibeg,
    iend (the last region); the search range imin..imax; the flags lcz; lok;
    m shells; hcomp[k, j]; the rows jbeg..jend. Arrays have index 0 unused.

    Sets in place: hcomp, and ecod2, em for the last region.

    Returns (ibeg, iend), the last region.
    """
    global ldum, ldum2
    if not lok:
        # Question (see the page): only the last region before the final
        # call.
        mixcom(dt, ecod2, em, hs2, ibeg, iend, lcz, lok, m, hcomp, jbeg, jend)
        return ibeg, iend
    # Find the unstable regions in order: ibeg is the zone below the first
    # non-zero v, iend the zone above the last.
    ldum = False
    istart = imin
    while True:                                          # label 60
        ltest = False
        ended = False
        for j in range(istart, imax + 1):
            if hv[j] > 0.0:
                ltest = True
                if not ldum:
                    # Start of an unstable region.
                    ldum = True
                    ibeg = j - 1
            elif ldum:
                # End of an unstable region.
                iend = j - 1
                ldum = False
                istart = j + 1
                ended = True
                break                                    # GOTO 80
        if not ended:
            # The last interface is unstable: iend is IMAX.
            if ldum:
                iend = imax
            istart = imax + 1
        # label 80: no non-zero v, exit.
        if not ltest:
            break                                        # GOTO 90
        ldum2 = mixgrid(cod2, hd, hl, hp, hr, hs, hs1, hs2, hstot, ibeg,
                        iend, lcz, m, ecod2, em, ldum2)
        # Question (see the page): a single unstable interface is skipped.
        if ldum2:
            if istart <= imax:
                continue                                 # GOTO 60
            break                                        # GOTO 90
        mixcom(dt, ecod2, em, hs2, ibeg, iend, lcz, lok, m, hcomp, jbeg, jend)
        # Back for the next region, if any.
        if istart > imax:
            break
    # label 90
    return ibeg, iend
