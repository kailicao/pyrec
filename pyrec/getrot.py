"""Python version of getrot.f (YREC 09c3b9b): the rotation curve from the
specific angular momentum."""

from pyrec.solid import solid
from pyrec.wczimp import wczimp


def getrot(hd, hjm, hr, hs, hs2, lcz, m, eta2, hi, omega, qiw, r0):
    """Turn the specific angular momentum j(m) into Omega(m), outward:
    each radiative shell alone with SOLID, each convection zone as a
    whole with WCZIMP when the next radiative shell (or the surface) is
    reached.

    Translates: getrot.f#GETROT (09c3b9b)

    Arguments: hd, hr, hs, hs2 (the structure), lcz (convective flags),
    m (number of shells).
    Arrays changed in place: hjm (rewritten as I Omega / dm), eta2, hi,
    omega, qiw, r0.

    Returns nothing.
    """
    lconvz = False
    jstart = 0
    for i in range(1, m + 1):  # label 20
        if not lconvz:
            if not lcz[i]:
                # A radiative shell.
                jstart = i
                jend = i
                solid(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega,
                      qiw, r0, m)
            else:
                # The start of a convection zone.
                # Question (see the page): a zone that starts at m is
                # never treated.
                lconvz = True
                jstart = i
        elif not lcz[i]:
            # The convection zone ends below this radiative shell.
            jend = i - 1
            wczimp(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega, qiw,
                   r0, m)
            jstart = i
            jend = i
            solid(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega, qiw,
                  r0, m)
            lconvz = False
        elif i == m:
            # A convection zone reaching the surface.
            jend = m
            wczimp(hd, hjm, hr, hs, hs2, jstart, jend, eta2, hi, omega, qiw,
                   r0, m)
