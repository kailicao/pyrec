"""Python version of ccoeft.f (YREC 09c3b9b): the implicit tridiagonal
system for one step of composition mixing on the even grid."""

from pyrec import tridi_mod


def ccoeft(cod, dr, dt, ecomp, em, ntot):
    """Fill the tridiagonal system a, b, c, d of tridi_mod for one fully
    implicit diffusion step of a mass fraction on the even grid in chi,
    with no flux across either end.

    Translates: ccoeft.f#CCOEFT (09c3b9b)

    Arguments: cod, the diffusion coefficients at the interfaces (cod[i]
    between points i-1 and i); dr, the grid spacing (MIXCOM passes the
    spacing in chi); dt, the time step (s); ecomp, the mass fractions at
    the start of the step; em, the shell masses (g); ntot, the number of
    points. None is set; the results go to tridi_mod.
    """
    a = tridi_mod.a
    b = tridi_mod.b
    c = tridi_mod.c
    d = tridi_mod.d
    # The composition at the start of the step.
    for i in range(1, ntot + 1):
        d[i] = ecomp[i]
    # First shell: no flow below it.
    fact0 = dt / dr
    fact = fact0 / em[1]
    a[1] = 0.0
    b[1] = 1.0 + fact * cod[2]
    c[1] = -fact * cod[2]
    for i in range(2, ntot):
        fact = fact0 / em[i]
        a[i] = -fact * cod[i]
        b[i] = 1.0 + fact * (cod[i] + cod[i + 1])
        c[i] = -fact * cod[i + 1]
    # Last shell: no flow above it.
    fact = fact0 / em[ntot]
    a[ntot] = -fact * cod[ntot]
    b[ntot] = 1.0 + fact * cod[ntot]
    c[ntot] = 0.0
