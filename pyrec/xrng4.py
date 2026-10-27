"""Python version of xrng4.f (YREC 09c3b9b): four consecutive table
indices around an index, for 4-point interpolation."""


def xrng4(jx, numx, jx1, jx2):
    """Choose four consecutive indices jx1..jx2 around jx (jx-1..jx+2),
    shifted to 1..4 near the start and numx-3..numx near the end.

    Translates: xrng4.f#XRNG4 (09c3b9b)

    Arguments: jx, the index; numx, the table's length; jx1 and jx2 are
    passed in with any value.

    Returns (jx1, jx2).
    """
    # Question (see the page): numx < 4 gives indices outside the table.
    if jx <= 2:
        jx1 = 1
        jx2 = 4
    elif jx >= numx - 2:
        jx1 = numx - 3
        jx2 = numx
    else:
        jx1 = jx - 1
        jx2 = jx + 2
    return jx1, jx2
