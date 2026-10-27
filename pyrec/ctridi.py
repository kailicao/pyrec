"""Python version of ctridi.f (YREC 09c3b9b): the tridiagonal solver
(Thomas algorithm) for the composition mixing step."""

from pyrec import tridi_mod


def ctridi(n):
    """Solve the tridiagonal system a, b, c, d of tridi_mod (set by
    ccoeft) by forward elimination and back substitution, without
    pivoting; the solution goes to tridi_mod.u.

    Translates: ctridi.f#CTRIDI (09c3b9b)

    Arguments: n, the number of points. None is returned.
    """
    a = tridi_mod.a
    b = tridi_mod.b
    c = tridi_mod.c
    d = tridi_mod.d
    u = tridi_mod.u
    gama = tridi_mod.gama
    bet = b[1]
    u[1] = d[1] / bet
    for j in range(2, n + 1):
        gama[j] = c[j - 1] / bet
        bet = b[j] - a[j] * gama[j]
        # Question (see the page): only an exact zero is caught.
        if bet == 0:
            raise SystemExit('#TRIDIA:SINGULAR MATRIX')
        u[j] = (d[j] - a[j] * u[j - 1]) / bet
    for j in range(n - 1, 0, -1):
        u[j] = u[j] - gama[j + 1] * u[j + 1]
