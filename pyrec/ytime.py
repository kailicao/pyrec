"""Python version of ytime.f (YREC 09c3b9b): the time step allowed by
helium burning."""

from pyrec import constants_mod
from pyrec import nmlpar_mod
from pyrec.engeb import engeb

# SAVE: the extended-composition abundances keep their values between
# calls, so that without LEXCOM ENGEB gets the last ones set (or 0).
xh2 = 0.0
xli6 = 0.0
xli7 = 0.0
xbe9 = 0.0


def ytime(eg, hcomp, hd, hl, hs1, ht, jcore, m, deltsy,
          hr1, hr2, hr3, hr4, hr5, hr6, hr7, hr8, hr9, hr10, hr11, hr12, hr13,
          hf1, hf2, jxbeg):
    """The helium-burning time step DELTSY: in the helium flash (central
    log rho >= 5), 1e15 over the triple-alpha rate at the temperature
    maximum; otherwise the time to burn a fraction of the core's helium.

    Translates: ytime.f#YTIME (09c3b9b)

    Arguments: the energy rates eg(1..6) (set here, in the flash), the model
    (hcomp, hd = log rho, hl, hs1 = m, ht = log T, jcore, m, jxbeg), the
    step deltsy, and the rate work arrays hr1..hr13, hf1, hf2 passed on to
    ENGEB. Returns deltsy.
    """
    global xh2, xli6, xli7, xbe9
    n = nmlpar_mod
    if hd[1] >= 5.0:
        # Helium flash: the energy rates at the temperature maximum.
        htmax = 0.0
        itmax = 1
        for it in range(1, m + 1):  # label 210
            if ht[it] > htmax:
                itmax = it
                htmax = ht[it]
        dl = hd[itmax]
        tl = ht[itmax]
        x = hcomp[1, itmax]
        y = hcomp[2, itmax]
        z = hcomp[3, itmax]
        xhe3 = hcomp[4, itmax]
        xc12 = hcomp[5, itmax]
        xc13 = hcomp[6, itmax]
        xn14 = hcomp[7, itmax]
        xn15 = hcomp[8, itmax]
        xo16 = hcomp[9, itmax]
        xo17 = hcomp[10, itmax]
        xo18 = hcomp[11, itmax]
        if n.lexcom:
            xh2 = hcomp[12, itmax]
            xli6 = hcomp[13, itmax]
            xli7 = hcomp[14, itmax]
            xbe9 = hcomp[15, itmax]
        iu = itmax
        # ENGEB's set scalars are passed in with any value; BE7MFLOC is not
        # used here.
        (egg11, egg22, egg33, egg44, egg55, qed, qet, sum1, be7mfloc,
         ealpca, enu) = engeb(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dl,
                              tl, x, y, xhe3, xc12, xc13, xn14, xo16,
                              xo18, xh2, iu, hr1, hr2, hr3, hr4, hr5, hr6,
                              hr7, hr8, hr9, hr10, hr11, hr12, hr13, hf1,
                              hf2, 0.0, 0.0, 0.0)
        e = sum1
        eg[1] = egg11
        eg[2] = egg22
        eg[3] = egg33
        eg[4] = egg44
        eg[5] = egg55
        eg[6] = enu
        eca = ealpca
        if eg[5] < 1.0e-22:
            eg[5] = 1.0e-22
        # Question (see the page): the constant 1e15 and the floor 1e-22.
        deltsy = 1.0e15 / eg[5]
    else:
        # Limits on helium burning for stars outside the flash.
        # Question (see the page): 1 - Z is X + Y, not Y.
        ycore = 1.0 - hcomp[3, jcore]
        if ycore >= n.atime[1]:
            deltsy = min(n.atime[4], n.atime[5] * ycore)
            deltsy = ((5.85e17 / constants_mod.clsun) * deltsy
                      * (hs1[jcore] / hl[jcore]))
        else:
            # Question (see the page): deltsy is the value passed in.
            deltsy = ((5.85e17 / constants_mod.clsun) * deltsy
                      * (hs1[jxbeg - 1] / hl[jxbeg - 1]))
    return deltsy
