"""Python version of simeqc.f (YREC 09c3b9b): a small linear system solved
in place by Gaussian elimination with partial pivoting."""

from pyrec import iounits_mod


def simeqc(a, m, n):
    """Solve n linear equations whose matrix and m - n right-hand sides are
    stored together by columns in a[1..56] (n rows, m columns); the
    solutions overwrite the right-hand-side columns. A zero pivot stops the
    run.

    Translates: simeqc.f#SIMEQC (09c3b9b)

    Arguments: a[1..56], changed in place; m, the number of columns; n, the
    number of equations. None is returned.
    """
    # Forward elimination, one column j at a time.
    jj = -n
    for j in range(1, n + 1):
        jy = j + 1
        jj = jj + n + 1
        biga = 0.0
        it = jj - j
        for i in range(j, n + 1):
            ij = it + i
            if abs(biga) >= abs(a[ij]):
                continue
            biga = a[ij]
            imax = i
        # Question (see the page): an exact-zero test only.
        if abs(biga) == 0.0:
            # Question (see the page): the message goes to unit 5.
            iounits_mod.write(5, " STOPPED AT 1010")  # label 1010
            raise SystemExit(29)
        # Swap the pivot row in and divide it by the pivot (label 1012).
        ia = j + n * (j - 2)
        it = imax - j
        for k in range(j, m + 1):
            ia = ia + n
            ib = ia + it
            save = a[ia]
            a[ia] = a[ib]
            a[ib] = save
            a[ia] = a[ia] / biga
        if j == n:
            break
        # Eliminate the column below the pivot.
        ia = n * (j - 1)
        for ix in range(jy, n + 1):
            ib = ia + ix
            it = j - ix
            for jx in range(jy, m + 1):
                ic = n * (jx - 1) + ix
                id_ = ic + it
                a[ic] = a[ic] - a[ib] * a[id_]

    # Back substitution (label 70).
    ny = n - 1
    it = n * n
    for j in range(1, ny + 1):
        ia = it - j
        ic = n * m
        ib = ic - j
        for k in range(1, j + 1):
            a[ib] = a[ib] - a[ia] * a[ic]
            n1 = m - 1
            ig = ib
            ih = ic
            while n1 > n:  # label 75
                ig = ig - n
                ih = ih - n
                a[ig] = a[ig] - a[ia] * a[ih]
                n1 = n1 - 1
            ia = ia - n  # label 78
            ic = ic - 1
