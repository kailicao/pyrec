"""Python version of hsolve.f (YREC 09c3b9b): the last stage of the
Henyey solution, the surface condition and the back-substitution."""


def hsolve(m, hha, hhb, hsrf):
    """Solve the 2x2 surface system for the two corrections at point m
    (qqqp, qqqt), then back-substitute inward through the eliminated
    rows; on return hhb[1:5, i] holds the four corrections at point i.

    Translates: hsolve.f#HSOLVE (09c3b9b)

    Arguments: m (number of points), hha (the eliminated coefficients,
    hha[4, 2, i]), hsrf (the six surface boundary-condition entries).
    Array changed in place: hhb (right-hand sides, then corrections).

    Returns nothing.
    """
    # The surface 2x2 system.
    # Question (see the page): no check for a zero pivot.
    div = 1.0 / (hsrf[1] - hha[3, 1, m])
    pivot = (hsrf[4] - hha[4, 1, m]) * div
    hhb[1, 1] = hsrf[3] - hhb[3, m]
    hhb[2, 1] = (hsrf[6] - hhb[4, m]) - pivot * hhb[1, 1]
    qtemp = hsrf[2] - hha[3, 2, m]
    qqqt = hhb[2, 1] / (hsrf[5] - hha[4, 2, m] - pivot * qtemp)
    qqqp = div * (hhb[1, 1] - qqqt * qtemp)
    # Inward from the surface.
    for i in range(2, m + 1):  # label 10
        im = m + 2 - i
        pstore = hhb[1, im] - qqqp * hha[1, 1, im] - qqqt * hha[1, 2, im]
        tstore = hhb[2, im] - qqqp * hha[2, 1, im] - qqqt * hha[2, 2, im]
        hhb[3, im] = hhb[3, im] - qqqp * hha[3, 1, im] - qqqt * hha[3, 2, im]
        hhb[4, im] = hhb[4, im] - qqqp * hha[4, 1, im] - qqqt * hha[4, 2, im]
        hhb[1, im] = qqqp
        hhb[2, im] = qqqt
        qqqp = pstore
        qqqt = tstore
    # The centre.
    hhb[1, 1] = qqqp
    hhb[2, 1] = qqqt
    hhb[3, 1] = hhb[3, 1] - qqqp * hha[3, 1, 1] - qqqt * hha[3, 2, 1]
    hhb[4, 1] = hhb[4, 1] - qqqp * hha[4, 1, 1] - qqqt * hha[4, 2, 1]
