"""Python version of setupv.f (YREC 09c3b9b): the structural factors at
the shell interfaces that the rotational-mixing criteria and circulation
velocities need, computed once per rotational-mixing step."""

import math

import numpy as np

from pyrec import constants_mod
from pyrec import difad3_mod
from pyrec import dimensions_mod
from pyrec import egridder_mod
from pyrec import intvar_mod
from pyrec import mdphy_mod
from pyrec import nmlpar_mod
from pyrec import oldrot2_mod
from pyrec import rotder_mod
from pyrec import taukh_mod
from pyrec import vfact_mod

# SAVE: the local arrays keep their values between calls (elements the
# current call does not set are read only if IMIN < 2 or IMAX > M).
cpmi = np.zeros(dimensions_mod.json + 1)   # c_P at the interfaces
opm = np.zeros(dimensions_mod.json + 1)    # opacity at the interfaces (unused)


def _lagrange(hru, rm, i, j, faci):
    """Fill faci[1..4, i], the 4-point Lagrange weights at rm[i] over the
    shells j .. j+3 (the Fortran repeats this block three times)."""
    a43 = hru[j + 3] - hru[j + 2]
    a42 = hru[j + 3] - hru[j + 1]
    a41 = hru[j + 3] - hru[j]
    a32 = hru[j + 2] - hru[j + 1]
    a31 = hru[j + 2] - hru[j]
    a21 = hru[j + 1] - hru[j]
    d1 = -a21*a31*a41
    d2 = a21*a32*a42
    d3 = -a31*a32*a43
    d4 = a41*a42*a43
    xa1 = rm[i] - hru[j]
    xa2 = rm[i] - hru[j + 1]
    xa3 = rm[i] - hru[j + 2]
    xa4 = rm[i] - hru[j + 3]
    faci[1, i] = (xa2*xa3*xa4)/d1
    faci[2, i] = (xa1*xa3*xa4)/d2
    faci[3, i] = (xa1*xa2*xa4)/d3
    faci[4, i] = (xa1*xa2*xa3)/d4


def _interp(y, faci, i, j):
    """y interpolated to interface i from shells j .. j+3, summed in the
    Fortran's order."""
    return (y[j]*faci[1, i] + y[j + 1]*faci[2, i]
            + y[j + 2]*faci[3, i] + y[j + 3]*faci[4, i])


def _interface_values(i, j, hd, hg, hl, hp, hs1, ht, faci):
    """The physical variables at interface i, from shells j .. j+3 (the
    Fortran repeats this block three times: first, last and inner
    interfaces). P, T and rho are interpolated in the log."""
    c = constants_mod
    iv = intvar_mod
    md = mdphy_mod
    iv.pm[i] = math.exp(c.cln*_interp(hp, faci, i, j))
    iv.tm[i] = math.exp(c.cln*_interp(ht, faci, i, j))
    iv.dm[i] = math.exp(c.cln*_interp(hd, faci, i, j))
    # del (radiative) is interpolated; del is set to min(dela, delr) later.
    iv.delmi[i] = _interp(md.delrm, faci, i, j)
    iv.delami[i] = _interp(md.delam, faci, i, j)
    iv.qdtmi[i] = _interp(md.qdtm, faci, i, j)       # d ln rho / d ln T
    iv.hs3[i] = _interp(hs1, faci, i, j)             # mass inside (g)
    iv.epsilm[i] = _interp(md.esumm, faci, i, j)     # energy generation
    iv.alm[i] = c.clsun*_interp(hl, faci, i, j)      # luminosity (erg/s)
    iv.hgm[i] = _interp(hg, faci, i, j)              # gravity
    opm[i] = _interp(md.om, faci, i, j)              # opacity
    cpmi[i] = _interp(md.cpm, faci, i, j)            # specific heat


def setupv(hd, hg, hl, hp, hr, hs1, ht, imin, imax, m, hru, qwrmax, faci):
    """Prepare, for the interfaces imin..imax (interface i lies between
    shells i-1 and i), the factors of the stability criteria and
    circulation velocities that do not change during a diffusion time
    step; VCIRC multiplies them by the current Omega gradients.

    Translates: setupv.f#SETUPV (09c3b9b)

    Arguments: the model's log rho (hd), gravity (hg), L in solar units
    (hl), log P (hp), log r (hr), enclosed mass in g (hs1), log T (ht),
    arrays with index 0 unused; the interface range imin, imax; m, the
    number of shells. Set in place: hru (unlogged radii), qwrmax (the
    largest stable d Omega / d ln r against the dynamical shear) and
    faci[1..4, i], the interpolation weights (owned by SECULR, read by
    VCIRC). Also set: the interface values in intvar_mod, the velocity
    factors in vfact_mod, taukh_mod, oldrot2_mod and, with ldifad, in
    difad3_mod and egridder_mod.

    Returns nothing.
    """
    c = constants_mod
    iv = intvar_mod
    md = mdphy_mod
    vf = vfact_mod
    tk = taukh_mod
    o2 = oldrot2_mod
    da = difad3_mod
    rd = rotder_mod
    eg = egridder_mod

    # Unlogged radii of the model points (label 5).
    for i in range(1, m + 1):
        hru[i] = math.exp(c.cln*hr[i])
    # Radii of the midpoints of the interfaces (label 10).
    for i in range(2, m + 1):
        iv.rm[i] = 0.5*(hru[i] + hru[i - 1])
    # 4-point Lagrange interpolation: the first interface uses the first 4
    # points, the last the last 4, and interface i uses i-2, i-1, i, i+1.
    if imin < 3:
        _lagrange(hru, iv.rm, 2, 1, faci)
        i0 = 3
        _interface_values(2, 1, hd, hg, hl, hp, hs1, ht, faci)
    else:
        i0 = imin
    if imax == m:
        _lagrange(hru, iv.rm, m, m - 3, faci)
        i1 = m - 1
        _interface_values(m, m - 3, hd, hg, hl, hp, hs1, ht, faci)
    else:
        i1 = imax
    # Interpolation factors for all other interfaces (label 20).
    for i in range(i0, i1 + 1):
        _lagrange(hru, iv.rm, i, i - 2, faci)
    cg = math.exp(c.cln*c.cgl)
    cg2 = cg**2   # unused, as in the Fortran
    # Physical variables at the interfaces (label 30).
    for i in range(i0, i1 + 1):
        _interface_values(i, i - 2, hd, hg, hl, hp, hs1, ht, faci)
    # The actual del is min(del_rad, del_ad) (label 35).
    for i in range(imin, imax + 1):
        iv.delmi[i] = min(iv.delmi[i], iv.delami[i])
    for i in range(1, m + 1):
        tk.fes1[i] = 0.0
        tk.fes2[i] = 0.0
        tk.fes3[i] = 0.0
        tk.fgsfj[i] = 0.0
        tk.fgsfk[i] = 0.0
        vf.fact1[i] = 0.0
        vf.fact2[i] = 0.0
        vf.fact3[i] = 0.0
        vf.fact4[i] = 0.0
        vf.fact5[i] = 0.0
        o2.festime[i] = 0.0
        o2.thm[i] = 0.0
        vf.fact7[i] = 0.0
        vf.fact8[i] = 0.0
    # The d theta/dt term of Zahn and Maeder 1998: r^4 / (sqrt(m) G).
    for i in range(2, m + 1):
        fac = (2.0*(hr[i] + hr[i - 1])
               - 0.5*(math.log10(hs1[i]) + math.log10(hs1[i - 1])) - c.cgl)
        o2.thm[i] = math.exp(c.cln*fac)

    # Structural quantities for the velocities at all interfaces (label 40).
    cpigi = 4.0/c.c4pi/cg
    for i in range(imin, imax + 1):
        if not nmlpar_mod.ldifad:
            delmu = ((math.log10(md.amum[i]) - math.log10(md.amum[i - 1]))
                     / (hp[i] - hp[i - 1]))
        else:
            delmu = 0.0
        if i == imin:
            dmin = max(1.0e-3, 0.5*(md.delam[imin + 1] - md.delm[imin + 1])
                       + delmu)
        elif i == imax:
            dmin = max(1.0e-3, 0.5*(md.delam[imax - 1] - md.delm[imax - 1])
                       + delmu)
        else:
            dmin = 1.0e-3
        ddel = max(iv.delami[i] - iv.delmi[i], dmin)
        o2.ddelm[i] = max(iv.delami[i] - iv.delmi[i] + delmu, dmin)
        o2.festime[i] = iv.pm[i]/(iv.hgm[i]*ddel*iv.dm[i])
        hpsc = iv.pm[i]*iv.rm[i]**2/iv.dm[i]/iv.hs3[i]/cg
        htsc = hpsc/iv.delmi[i]
        # Question (see the page): the factor 2 in eta.
        eta = 2.0*c.cc23*iv.rm[i]**3/cg/iv.hs3[i]
        tk.fgsfk[i] = 8.0*hpsc*eta/ddel
        ff = iv.pm[i]/(iv.hgm[i]*ddel*cpmi[i]*iv.dm[i]*iv.tm[i])
        fac = 2.0*eta*ff
        em = iv.alm[i]/iv.hs3[i]
        tk.fes1[i] = fac*(em - iv.epsilm[i])
        tk.fes2[i] = -0.5*fac*em*cpigi/iv.dm[i]
        tk.fes3[i] = em*ff*c.cc13*cpigi*iv.rm[i]/iv.dm[i]
        da.fv0[i] = ff
        da.fv1a[i] = tk.fes1[i]/ff
        da.fv1b[i] = tk.fes2[i]/ff
        da.fv2a[i] = em*cpigi*3.0*htsc/iv.dm[i]
        da.fv2b[i] = -em*cpigi*iv.rm[i]/iv.dm[i]
        tk.fgsfj[i] = abs(fac*em*htsc)*iv.rm[i]
        # Eddington circulation: v_ES = fact1*fact2*Omega**2 (Endal and
        # Sofia, paper II).
        # Question (see the page): the floor 1e-3 caps 1/(del_ad - del)
        # at 1e3, not the comment's 1e6.
        vf.fact1[i] = 1.0/iv.qdtmi[i]/max(iv.delami[i] - iv.delmi[i], 1.0e-3)
        fac = iv.delami[i]/iv.hgm[i]
        vf.fact2[i] = fac*eta*(iv.epsilm[i] - em)
        tk.fact6[i] = -fac*em/c.c4pi/cg/iv.dm[i]/iv.rm[i]
        # GSF: fact3 is the temperature scale height times the radius.
        vf.fact3[i] = abs(htsc*iv.rm[i])
        # Fictitious mu currents (Kippenhahn, IAU 66).
        vf.fact4[i] = (iv.rm[i]**2/iv.hs3[i]/cg/iv.dm[i]
                       * (iv.pm[i] - c.ca3*iv.tm[i]**4)/math.sqrt(iv.hgm[i]))
        # Dynamical shear: stable when the Richardson number exceeds 1/4.
        qwrmax[i] = (math.sqrt(iv.dm[i]/iv.pm[i]
                               * max(iv.delami[i] - iv.delmi[i], 1.0e-3))
                     * 2.0*iv.hgm[i])
        # Secular shear inhibited by mu gradients.
        vf.fact5[i] = (iv.dm[i]*iv.hgm[i]**2/iv.pm[i]/c.cln
                       / (hp[i] - hp[i - 1])
                       * (c.ca3*iv.tm[i]**4/iv.pm[i] - 1.0)
                       / (1.0 + 3.0*c.ca3*iv.tm[i]**4/iv.pm[i]))

    # Terms of order d Omega/dr from Zahn and Maeder 1998.
    if nmlpar_mod.ldifad:
        # Question (see the page): reads shell imin-1, index 0 if imin = 1.
        htm = (math.exp(c.cln*(hp[imin - 1] + 2.0*hr[imin - 1]
                               - hd[imin - 1]))
               / hs1[imin - 1]/cg/md.delrm[imin - 1])
        eg.qqqcod[imin - 1] = htm
        for i in range(imin, imax + 1):
            hpsc = iv.pm[i]*iv.rm[i]**2/iv.dm[i]/iv.hs3[i]/cg
            htsc = hpsc/iv.delmi[i]
            f1 = iv.pm[i]/(iv.hgm[i]*o2.ddelm[i]*cpmi[i]*iv.dm[i]*iv.tm[i])
            f2 = iv.alm[i]*iv.rm[i]/iv.hs3[i]/3.0
            f3 = 0.75*iv.hs3[i]/(c.cpi*iv.dm[i]*iv.rm[i]**3)
            # Question (see the page): v0, qc1r, qqc1rr, qqchitr, htp2 and
            # qhtr feed only commented-out formulas.
            v0 = -f1*f2*f3
            ff = iv.rm[i]**3/iv.hs3[i]
            c1 = c.cc23*iv.rm[i]**4/cg/iv.hs3[i]
            eg.qqcod[i] = c1
            qc1r = 4.0*c.cc23*ff/cg*(1.0 - c.cpi*iv.dm[i]*ff)
            dr = 10.0**hr[i] - 10.0**hr[i - 1]
            qdr = (10.0**hd[i] - 10.0**hd[i - 1])/dr
            qqc1rr = (8.0*ff/cg/iv.rm[i]
                      * (1.0 - c.cc13*c.cpi*ff*(1.0e1*iv.dm[i] - iv.rm[i]*qdr)
                         + c.cc13*8.0*(c.cpi*iv.dm[i]*ff)**2))
            # d ln chi / d ln T = 3 - d ln kappa / d ln T
            qchit = 3.0 - 0.5*(rd.sqot[i] + rd.sqot[i - 1])
            qqchitr = (rd.sqot[i - 1] - rd.sqot[i])/dr
            htp = (math.exp(c.cln*(hp[i] + 2.0*hr[i] - hd[i]))
                   / hs1[i]/cg/md.delm[i])
            eg.qqqcod[i] = htp
            htp2 = dr/c.cln/(ht[i - 1] - ht[i])
            qhtr = (abs(htp) - abs(htm))/dr
            htm = htp
            qetm = 0.5*(rd.sqet[i] + rd.sqet[i - 1])
            fnu = 0.5*(rd.fracnuc[i] + rd.fracnuc[i - 1])
            vf.fact7[i] = (-f1*c1*(fnu*qetm
                                   + iv.epsilm[i]*(1.0 - fnu - qchit)))
            da.fv2b[i] = da.fv2b[i] + vf.fact7[i]/f1
            vf.fact8[i] = 0.0
            fl = iv.alm[i]/c.c4pi/iv.dm[i]/iv.rm[i]**2
            da.facd2[i] = f1*(fl*qchit - iv.epsilm[i]*htsc)
            da.facd3[i] = -f1*fl
