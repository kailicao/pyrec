"""Python version of trapzd.f (YREC 7be4133): one refinement of the
trapezoidal rule for FPFT's integral between two shells."""

# The Fortran's SAVEd local: the number of intervals of the next refinement.
it = 0


def trapzd(b1: float, b2: float, s: float, n: int, rho: float, rhop: float, sm: float, smp: float, w2: float,
           w2p: float, eta22: float, eta22p: float, q: float, qp: float) -> tuple[float, float]:
    """One refinement step of the trapezoidal rule (Numerical Recipes'
    trapzd) for the integral of rho*Omega^2*r0^7/m*(5+eta2)/(2+eta2)
    between r0 = b1 and b2: the first call (n = 1) takes the two ends,
    each later call adds the midpoints of the current intervals.

    Translates: trapzd.f#TRAPZD (7be4133)

    Arguments: the shells' radii b1 (inner), b2 (outer); s, the running
    estimate; n, the refinement (1, 2, 3, ...); rho, sm, w2, eta22 at the
    outer shell and rhop, smp, w2p, eta22p at the inner one (eta22 is
    eta2 + 2); q, the integrand at the outer shell (set when n = 1, passed
    in with any value otherwise); qp, the integrand at the inner shell.

    Returns (s, q).

    """

    global it
    r0 = b2
    r03 = r0**3

    if n == 1:
        # The first call: the integrand at the outer shell and the
        # trapezoid over the whole step (Fortran lines 22–30).
        q = (rho * w2 * r03 * (3.0 + eta22) / (sm * eta22)) * r03 * r0
        s = 0.5 * (b2 - b1) * (qp + q)
        it = 1
    else:
        # A refinement: add the midpoints of the current intervals
        # (Fortran lines 31–55).
        tnm = float(it)
        dr = b2 - b1
        delta = dr / tnm                 # DEL
        y = b1 + 0.5 * delta
        total = 0.0                      # SUM
        drho = (rho - rhop) / dr
        dm = (sm - smp) / (b2**2 - b1**2)
        deta2 = (eta22 - eta22p) / dr
        dw2 = (w2 - w2p) / dr
        for j in range(1, it + 1):      # loop 10
            r03t = y**3
            # Suspected bug (see the page): rho, w2 and eta2 are interpolated
            # at b1 + delta for every midpoint; only m uses y.
            rhot = rhop + drho * delta
            smt = smp + dm * (y**2 - b1**2)
            w2t = w2p + dw2 * delta
            eta22t = eta22p + deta2 * delta
            q0 = (rhot * w2t * r03t * (3.0 + eta22t) / (smt * eta22t)) * r03t * y
            total = total + q0
            y = y + delta
        s = 0.5 * (s + delta * total)
        it = it * 2
    return s, q
