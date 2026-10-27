"""Python version of sneut.f (YREC 09c3b9b): the thermal neutrino
energy-loss rate from the fits of Itoh et al. (1996, ApJS 102, 411), with
its derivatives."""

from math import cos, exp, log10, sin, sqrt

from pyrec.ifermi12 import ifermi12
from pyrec.zfermim12 import zfermim12

# The Fortran SAVEs all its locals, but each is set before it is read, so
# no state is kept between calls here.

# Mathematical constants (PARAMETER).
pi = 3.1415926535897932384
fac1 = 5.0 * pi / 3.0
fac2 = 10.0 * pi
fac3 = pi / 5.0
oneth = 1.0/3.0
twoth = 2.0/3.0
con1 = 1.0/5.9302
sixth = 1.0/6.0
iln10 = 4.342944819032518e-1

# Weinberg angle, number of neutrino flavors and the coupling factors
# (PARAMETER): tfac1 = C+, tfac2 = C-, tfac3 = C-/C+, tfac4 = C+/2,
# tfac5 = C-/2.
theta = 0.2319
xnufam = 3.0
cv = 0.5 + 2.0 * theta
cvp = 1.0 - cv
ca = 0.5
cap = 1.0 - ca
tfac1 = cv*cv + ca*ca + (xnufam-1.0) * (cvp*cvp+cap*cap)
# Question (see the page): cap-cap (= 0) where cap*cap is expected.
tfac2 = cv*cv - ca*ca + (xnufam-1.0) * (cvp*cvp - cap-cap)
tfac3 = tfac2/tfac1
tfac4 = 0.5 * tfac1
tfac5 = 0.5 * tfac2
tfac6 = cv*cv + 1.5*ca*ca + (xnufam - 1.0)*(cvp*cvp + 1.5*cap*cap)


def sneut(temp, den, abar, zbar, snu, dsnudt, dsnudd, dsnuda, dsnudz):
    """The thermal neutrino loss of stellar matter, the sum of pair,
    plasma, photo, bremsstrahlung and recombination neutrinos, with its
    derivatives by the chain rule; zero below 1e7 K.

    Translates: sneut.f#SNEUT (09c3b9b)

    Arguments: temp (K), den (g/cm^3), abar and zbar (mean atomic weight
    and charge); snu and its derivatives are only set.

    Returns (snu, dsnudt, dsnudd, dsnuda, dsnudz): the loss in erg/g/s and
    its derivatives with respect to temp, den, abar and zbar.
    """
    # Initialize.
    spair = 0.0
    spairdt = 0.0
    spairdd = 0.0
    spairda = 0.0
    spairdz = 0.0

    splas = 0.0
    splasdt = 0.0
    splasdd = 0.0
    splasda = 0.0
    splasdz = 0.0

    sphot = 0.0
    sphotdt = 0.0
    sphotdd = 0.0
    sphotda = 0.0
    sphotdz = 0.0

    sbrem = 0.0
    sbremdt = 0.0
    sbremdd = 0.0
    sbremda = 0.0
    sbremdz = 0.0

    sreco = 0.0
    srecodt = 0.0
    srecodd = 0.0
    srecoda = 0.0
    srecodz = 0.0

    snu = 0.0
    dsnudt = 0.0
    dsnudd = 0.0
    dsnuda = 0.0
    dsnudz = 0.0

    if temp < 1.0e7:
        return snu, dsnudt, dsnudd, dsnuda, dsnudz

    # Some frequent factors.
    deni = 1.0/den
    tempi = 1.0/temp
    abari = 1.0/abar
    zbari = 1.0/zbar

    # Electron mole number and mean molecular weight per electron.
    ye = zbar*abari
    xmue = abar*zbari  # unused, as in the Fortran

    # Temperature and its powers; xl = kT/(m_e c^2).
    t9 = temp * 1.0e-9
    xl = t9 * con1
    xldt = 1.0e-9 * con1
    xlp5 = sqrt(xl)
    xl2 = xl*xl
    xl3 = xl2*xl
    xl4 = xl3*xl
    xl5 = xl4*xl
    xl6 = xl5*xl
    xl7 = xl6*xl
    xl8 = xl7*xl
    xl9 = xl8*xl
    xlmp5 = 1.0/xlp5
    xlm1 = 1.0/xl
    xlm2 = xlm1*xlm1
    xlm3 = xlm1*xlm2
    xlm4 = xlm1*xlm3

    # rho * ye and its derivatives.
    rm = den*ye
    rmdd = ye
    rmda = -rm*abari
    rmdz = den*abari
    rmi = 1.0/rm

    # zeta = (rho ye / 1e9)^(1/3) / lambda.
    a0 = rm * 1.0e-9
    a1 = a0**oneth
    zeta = a1 * xlm1
    zetadt = -a1 * xlm2 * xldt
    a2 = oneth * a1*rmi * xlm1
    zetadd = a2 * rmdd
    zetada = a2 * rmda
    zetadz = a2 * rmdz

    zeta2 = zeta * zeta
    zeta3 = zeta2 * zeta

    # ---- Pair neutrinos ----
    # g(lambda), eq. 2.8.
    gl = 1.0 - 13.04*xl2 + 133.5*xl4 + 1534.0*xl6 + 918.6*xl8
    gldt = xldt*(-26.08*xl + 534.0*xl3 + 9204.0*xl5 + 7348.8*xl7)

    # f(pair), eq. 2.7: numerator.
    a1 = 6.002e19 + 2.084e20*zeta + 1.872e21*zeta2
    a2 = 2.084e20 + 2.0*1.872e21*zeta

    if t9 < 10.0:
        b1 = exp(-5.5924*zeta)
        b2 = -b1*5.5924
    else:
        b1 = exp(-4.9924*zeta)
        b2 = -b1*4.9924

    xnum = a1 * b1
    c = a2*b1 + a1*b2
    xnumdt = c*zetadt
    xnumdd = c*zetadd
    xnumda = c*zetada
    xnumdz = c*zetadz

    # Denominator.
    if t9 < 10.0:
        a1 = 9.383e-1*xlm1 - 4.141e-1*xlm2 + 5.829e-2*xlm3
        a2 = -9.383e-1*xlm2 + 2.0*4.141e-1*xlm3 - 3.0*5.829e-2*xlm4
    else:
        a1 = 1.2383*xlm1 - 8.141e-1*xlm2
        a2 = -1.2383*xlm2 + 2.0*8.141e-1*xlm3

    b1 = 3.0*zeta2

    xden = zeta3 + a1
    xdendt = b1*zetadt + a2*xldt
    xdendd = b1*zetadd
    xdenda = b1*zetada
    xdendz = b1*zetadz

    a1 = 1.0/xden
    fpair = xnum*a1
    fpairdt = (xnumdt - fpair*xdendt)*a1
    fpairdd = (xnumdd - fpair*xdendd)*a1
    fpairda = (xnumda - fpair*xdenda)*a1
    fpairdz = (xnumdz - fpair*xdendz)*a1

    # q(pair), eq. 2.6.
    a1 = 10.7480*xl2 + 0.3967*xlp5 + 1.005
    a2 = xldt*(2.0*10.7480*xl + 0.5*0.3967*xlmp5)
    xnum = 1.0/a1
    xnumdt = -xnum*xnum*a2

    a1 = 7.692e7*xl3 + 9.715e6*xlp5
    a2 = xldt*(3.0*7.692e7*xl2 + 0.5*9.715e6*xlmp5)

    c = 1.0/a1
    b1 = 1.0 + rm*c

    xden = b1**(-0.3)

    d = -0.3*xden/b1
    xdendt = -d*rm*c*c*a2
    xdendd = d*rmdd*c
    xdenda = d*rmda*c
    xdendz = d*rmdz*c

    qpair = xnum*xden
    qpairdt = xnumdt*xden + xnum*xdendt
    qpairdd = xnum*xdendd
    qpairda = xnum*xdenda
    qpairdz = xnum*xdendz

    # Eq. 2.5: exp(-2/lambda) * f(pair) ...
    a1 = exp(-2.0*xlm1)
    a2 = a1*2.0*xlm2*xldt

    spair = a1*fpair
    spairdt = a2*fpair + a1*fpairdt
    spairdd = a1*fpairdd
    spairda = a1*fpairda
    spairdz = a1*fpairdz

    # ... times g(lambda) ...
    a1 = spair
    spair = gl*a1
    spairdt = gl*spairdt + gldt*a1
    spairdd = gl*spairdd
    spairda = gl*spairda
    spairdz = gl*spairdz

    # ... times (C+ + C- q(pair)) / 2.
    a1 = tfac4*(1.0 + tfac3 * qpair)
    a2 = tfac4*tfac3

    a3 = spair
    spair = a1*a3
    spairdt = a1*spairdt + a2*qpairdt*a3
    spairdd = a1*spairdd + a2*qpairdd*a3
    spairda = a1*spairda + a2*qpairda*a3
    spairdz = a1*spairdz + a2*qpairdz*a3

    # ---- Plasma neutrinos ----
    # gamma, eq. 4.6.
    a1 = 1.019e-6*rm
    a2 = a1**twoth
    a3 = twoth*a2/a1

    b1 = sqrt(1.0 + a2)
    b2 = 1.0/b1

    c00 = 1.0/(temp*temp*b1)

    gl2 = 1.1095e11 * rm * c00

    gl2dt = -2.0*gl2*tempi
    d = rm*c00*b2*0.5*b2*a3*1.019e-6
    gl2dd = 1.1095e11 * (rmdd*c00 - d*rmdd)
    gl2da = 1.1095e11 * (rmda*c00 - d*rmda)
    gl2dz = 1.1095e11 * (rmdz*c00 - d*rmdz)

    gl = sqrt(gl2)
    gl12 = sqrt(gl)
    gl32 = gl * gl12
    gl72 = gl2 * gl32
    gl6 = gl2 * gl2 * gl2

    # Transverse, eq. 4.7.
    ft = 2.4 + 0.6*gl12 + 0.51*gl + 1.25*gl32
    gum = 1.0/gl2
    a1 = (0.25*0.6*gl12 + 0.5*0.51*gl + 0.75*1.25*gl32)*gum
    ftdt = a1*gl2dt
    ftdd = a1*gl2dd
    ftda = a1*gl2da
    ftdz = a1*gl2dz

    # Longitudinal, eq. 4.8.
    a1 = 8.6*gl2 + 1.35*gl72
    a2 = 8.6 + 1.75*1.35*gl72*gum

    b1 = 225.0 - 17.0*gl + gl2
    b2 = -0.5*17.0*gl*gum + 1.0

    c = 1.0/b1
    fl = a1*c

    d = (a2 - fl*b2)*c
    fldt = d*gl2dt
    fldd = d*gl2dd
    flda = d*gl2da
    fldz = d*gl2dz

    # xnum and xden, eqs. 4.9 and 4.10.
    cc = log10(2.0*rm)
    xlnt = log10(temp)

    xnum = sixth * (17.5 + cc - 3.0*xlnt)
    xnumdt = -iln10*0.5*tempi
    a2 = iln10*sixth*rmi
    xnumdd = a2*rmdd
    xnumda = a2*rmda
    xnumdz = a2*rmdz

    xden = sixth * (-24.5 + cc + 3.0*xlnt)
    xdendt = iln10*0.5*tempi
    xdendd = a2*rmdd
    xdenda = a2*rmda
    xdendz = a2*rmdz

    # The correction fxy, eq. 4.11.
    if abs(xnum) > 0.7 or xden < 0.0:
        fxy = 1.0
        fxydt = 0.0
        fxydd = 0.0
        fxydz = 0.0
        fxyda = 0.0
    else:
        a1 = 0.39 - 1.25*xnum - 0.35*sin(4.5*xnum)
        a2 = -1.25 - 4.5*0.35*cos(4.5*xnum)

        b1 = 0.3 * exp(-1.0*(4.5*xnum + 0.9)**2)
        b2 = -b1*2.0*(4.5*xnum + 0.9)*4.5

        c = min(0.0, xden - 1.6 + 1.25*xnum)
        if c == 0.0:
            dumdt = 0.0
            dumdd = 0.0
            dumda = 0.0
            dumdz = 0.0
        else:
            dumdt = xdendt + 1.25*xnumdt
            dumdd = xdendd + 1.25*xnumdd
            dumda = xdenda + 1.25*xnumda
            dumdz = xdendz + 1.25*xnumdz

        d = 0.57 - 0.25*xnum
        a3 = c/d
        c00 = exp(-1.0*a3**2)

        f1 = -c00*2.0*a3/d
        c01 = f1*(dumdt + a3*0.25*xnumdt)
        c02 = f1*(dumdd + a3*0.25*xnumdd)
        c03 = f1*(dumda + a3*0.25*xnumda)
        c04 = f1*(dumdz + a3*0.25*xnumdz)

        fxy = 1.05 + (a1 - b1)*c00
        fxydt = (a2*xnumdt - b2*xnumdt)*c00 + (a1-b1)*c01
        fxydd = (a2*xnumdd - b2*xnumdd)*c00 + (a1-b1)*c02
        fxyda = (a2*xnumda - b2*xnumda)*c00 + (a1-b1)*c03
        fxydz = (a2*xnumdz - b2*xnumdz)*c00 + (a1-b1)*c04

    # Eqs. 4.1 and 4.5: (ft + fl) * fxy ...
    splas = (ft + fl) * fxy
    splasdt = (ftdt + fldt)*fxy + (ft+fl)*fxydt
    splasdd = (ftdd + fldd)*fxy + (ft+fl)*fxydd
    splasda = (ftda + flda)*fxy + (ft+fl)*fxyda
    splasdz = (ftdz + fldz)*fxy + (ft+fl)*fxydz

    # ... times exp(-gamma) ...
    a2 = exp(-gl)
    a3 = -0.5*a2*gl*gum

    a1 = splas
    splas = a2*a1
    splasdt = a2*splasdt + a3*gl2dt*a1
    splasdd = a2*splasdd + a3*gl2dd*a1
    splasda = a2*splasda + a3*gl2da*a1
    splasdz = a2*splasdz + a3*gl2dz*a1

    # ... times gamma^6 ...
    a2 = gl6
    a3 = 3.0*gl6*gum

    a1 = splas
    splas = a2*a1
    splasdt = a2*splasdt + a3*gl2dt*a1
    splasdd = a2*splasdd + a3*gl2dd*a1
    splasda = a2*splasda + a3*gl2da*a1
    splasdz = a2*splasdz + a3*gl2dz*a1

    # ... times the prefactor and lambda^9.
    a2 = 0.93153 * 3.0e21 * xl9
    a3 = 0.93153 * 3.0e21 * 9.0*xl8*xldt

    a1 = splas
    splas = a2*a1
    splasdt = a2*splasdt + a3*a1
    splasdd = a2*splasdd
    splasda = a2*splasda
    splasdz = a2*splasdz

    # ---- Photoneutrinos ----
    # tau (eq. 3.8), cc (eq. 3.6) and table 2, for three temperature ranges.
    if temp >= 1.0e7 and temp < 1.0e8:
        tau = log10(temp * 1.0e-7)
        cc = 0.5654 + tau
        c00 = 1.008e11
        c01 = 0.0
        c02 = 0.0
        c03 = 0.0
        c04 = 0.0
        c05 = 0.0
        c06 = 0.0
        c10 = 8.156e10
        c11 = 9.728e8
        c12 = -3.806e9
        c13 = -4.384e9
        c14 = -5.774e9
        c15 = -5.249e9
        c16 = -5.153e9
        c20 = 1.067e11
        c21 = -9.782e9
        c22 = -7.193e9
        c23 = -6.936e9
        c24 = -6.893e9
        c25 = -7.041e9
        c26 = -7.193e9
        dd01 = 0.0
        dd02 = 0.0
        dd03 = 0.0
        dd04 = 0.0
        dd05 = 0.0
        dd11 = -1.879e10
        dd12 = -9.667e9
        dd13 = -5.602e9
        dd14 = -3.370e9
        dd15 = -1.825e9
        dd21 = -2.919e10
        dd22 = -1.185e10
        dd23 = -7.270e9
        dd24 = -4.222e9
        dd25 = -1.560e9

    elif temp >= 1.0e8 and temp < 1.0e9:
        tau = log10(temp * 1.0e-8)
        cc = 1.5654
        c00 = 9.889e10
        c01 = -4.524e8
        c02 = -6.088e6
        c03 = 4.269e7
        c04 = 5.172e7
        c05 = 4.910e7
        c06 = 4.388e7
        c10 = 1.813e11
        c11 = -7.556e9
        c12 = -3.304e9
        c13 = -1.031e9
        c14 = -1.764e9
        c15 = -1.851e9
        c16 = -1.928e9
        c20 = 9.750e10
        c21 = 3.484e10
        c22 = 5.199e9
        c23 = -1.695e9
        c24 = -2.865e9
        c25 = -3.395e9
        c26 = -3.418e9
        dd01 = -1.135e8
        dd02 = 1.256e8
        dd03 = 5.149e7
        dd04 = 3.436e7
        dd05 = 1.005e7
        dd11 = 1.652e9
        dd12 = -3.119e9
        dd13 = -1.839e9
        dd14 = -1.458e9
        dd15 = -8.956e8
        dd21 = -1.549e10
        dd22 = -9.338e9
        dd23 = -5.899e9
        dd24 = -3.035e9
        dd25 = -1.598e9

    elif temp >= 1.0e9:
        tau = log10(t9)
        cc = 1.5654
        c00 = 9.581e10
        c01 = 4.107e8
        c02 = 2.305e8
        c03 = 2.236e8
        c04 = 1.580e8
        c05 = 2.165e8
        c06 = 1.721e8
        c10 = 1.459e12
        c11 = 1.314e11
        c12 = -1.169e11
        c13 = -1.765e11
        c14 = -1.867e11
        c15 = -1.983e11
        c16 = -1.896e11
        c20 = 2.424e11
        c21 = -3.669e9
        c22 = -8.691e9
        c23 = -7.967e9
        c24 = -7.932e9
        c25 = -7.987e9
        c26 = -8.333e9
        dd01 = 4.724e8
        dd02 = 2.976e8
        dd03 = 2.242e8
        dd04 = 7.937e7
        dd05 = 4.859e7
        dd11 = -7.094e11
        dd12 = -3.697e11
        dd13 = -2.189e11
        dd14 = -1.273e11
        dd15 = -5.705e10
        dd21 = -2.254e10
        dd22 = -1.551e10
        dd23 = -7.793e9
        dd24 = -4.489e9
        dd25 = -2.185e9

    taudt = iln10*tempi

    # The trigonometric terms of eq. 3.7, computed once.
    cos1 = cos(fac1*tau)
    cos2 = cos(fac1*2.0*tau)
    cos3 = cos(fac1*3.0*tau)
    cos4 = cos(fac1*4.0*tau)
    cos5 = cos(fac1*5.0*tau)
    last = cos(fac2*tau)

    sin1 = sin(fac1*tau)
    sin2 = sin(fac1*2.0*tau)
    sin3 = sin(fac1*3.0*tau)
    sin4 = sin(fac1*4.0*tau)
    sin5 = sin(fac1*5.0*tau)
    xast = sin(fac2*tau)

    # The three fit functions a0, a1, a2 and their temperature
    # derivatives f0, f1, f2.
    a0 = (0.5*c00
          + c01*cos1 + dd01*sin1 + c02*cos2 + dd02*sin2
          + c03*cos3 + dd03*sin3 + c04*cos4 + dd04*sin4
          + c05*cos5 + dd05*sin5 + 0.5*c06*last)

    f0 = (taudt*fac1*(-c01*sin1 + dd01*cos1 - c02*sin2*2.0
                      + dd02*cos2*2.0 - c03*sin3*3.0 + dd03*cos3*3.0
                      - c04*sin4*4.0 + dd04*cos4*4.0
                      - c05*sin5*5.0 + dd05*cos5*5.0)
          - 0.5*c06*xast*fac2*taudt)

    a1 = (0.5*c10
          + c11*cos1 + dd11*sin1 + c12*cos2 + dd12*sin2
          + c13*cos3 + dd13*sin3 + c14*cos4 + dd14*sin4
          + c15*cos5 + dd15*sin5 + 0.5*c16*last)

    f1 = (taudt*fac1*(-c11*sin1 + dd11*cos1 - c12*sin2*2.0
                      + dd12*cos2*2.0 - c13*sin3*3.0 + dd13*cos3*3.0
                      - c14*sin4*4.0 + dd14*cos4*4.0 - c15*sin5*5.0
                      + dd15*cos5*5.0) - 0.5*c16*xast*fac2*taudt)

    a2 = (0.5*c20
          + c21*cos1 + dd21*sin1 + c22*cos2 + dd22*sin2
          + c23*cos3 + dd23*sin3 + c24*cos4 + dd24*sin4
          + c25*cos5 + dd25*sin5 + 0.5*c26*last)

    f2 = (taudt*fac1*(-c21*sin1 + dd21*cos1 - c22*sin2*2.0
                      + dd22*cos2*2.0 - c23*sin3*3.0 + dd23*cos3*3.0
                      - c24*sin4*4.0 + dd24*cos4*4.0 - c25*sin5*5.0
                      + dd25*cos5*5.0) - 0.5*c26*xast*fac2*taudt)

    # f(photo), eq. 3.4: numerator ...
    dum = a0 + a1*zeta + a2*zeta2
    dumdt = f0 + f1*zeta + a1*zetadt + f2*zeta2 + 2.0*a2*zeta*zetadt
    dumdd = a1*zetadd + 2.0*a2*zeta*zetadd
    dumda = a1*zetada + 2.0*a2*zeta*zetada
    dumdz = a1*zetadz + 2.0*a2*zeta*zetadz

    z = exp(-cc*zeta)

    xnum = dum*z
    xnumdt = dumdt*z - dum*z*cc*zetadt
    xnumdd = dumdd*z - dum*z*cc*zetadd
    xnumda = dumda*z - dum*z*cc*zetada
    xnumdz = dumdz*z - dum*z*cc*zetadz

    # ... denominator ...
    xden = zeta3 + 6.290e-3*xlm1 + 7.483e-3*xlm2 + 3.061e-4*xlm3

    dum = 3.0*zeta2
    xdendt = dum*zetadt - xldt*(6.290e-3*xlm2
                                + 2.0*7.483e-3*xlm3 + 3.0*3.061e-4*xlm4)
    xdendd = dum*zetadd
    xdenda = dum*zetada
    xdendz = dum*zetadz

    dum = 1.0/xden
    fphot = xnum*dum
    fphotdt = (xnumdt - fphot*xdendt)*dum
    fphotdd = (xnumdd - fphot*xdendd)*dum
    fphotda = (xnumda - fphot*xdenda)*dum
    fphotdz = (xnumdz - fphot*xdendz)*dum

    # q(photo), eq. 3.3: numerator ...
    a0 = 1.0 + 2.045 * xl
    xnum = 0.666*a0**(-2.066)
    xnumdt = -2.066*xnum/a0 * 2.045*xldt

    # ... denominator ...
    dum = 1.875e8*xl + 1.653e8*xl2 + 8.449e8*xl3 - 1.604e8*xl4
    dumdt = xldt*(1.875e8 + 2.0*1.653e8*xl + 3.0*8.449e8*xl2
                  - 4.0*1.604e8*xl3)

    z = 1.0/dum
    xden = 1.0 + rm*z
    xdendt = -rm*z*z*dumdt
    xdendd = rmdd*z
    xdenda = rmda*z
    xdendz = rmdz*z

    z = 1.0/xden
    qphot = xnum*z
    qphotdt = (xnumdt - qphot*xdendt)*z
    dum = -qphot*z
    qphotdd = dum*xdendd
    qphotda = dum*xdenda
    qphotdz = dum*xdendz

    # Eq. 3.2: lambda^5 * f(photo) ...
    sphot = xl5 * fphot
    sphotdt = 5.0*xl4*xldt*fphot + xl5*fphotdt
    sphotdd = xl5*fphotdd
    sphotda = xl5*fphotda
    sphotdz = xl5*fphotdz

    # ... times rho*ye ...
    a1 = sphot
    sphot = rm*a1
    sphotdt = rm*sphotdt
    sphotdd = rm*sphotdd + rmdd*a1
    sphotda = rm*sphotda + rmda*a1
    sphotdz = rm*sphotdz + rmdz*a1

    # ... times (C+ - C- q(photo)) / 2.
    a1 = tfac4*(1.0 - tfac3 * qphot)
    a2 = -tfac4*tfac3

    a3 = sphot
    sphot = a1*a3
    sphotdt = a1*sphotdt + a2*qphotdt*a3
    sphotdd = a1*sphotdd + a2*qphotdd*a3
    sphotda = a1*sphotda + a2*qphotda*a3
    sphotdz = a1*sphotdz + a2*qphotdz*a3

    if sphot <= 0.0:
        sphot = 0.0
        sphotdt = 0.0
        sphotdd = 0.0
        sphotda = 0.0
        sphotdz = 0.0

    # ---- Bremsstrahlung neutrinos ----
    den6 = den * 1.0e-6
    t8 = temp * 1.0e-8
    t812 = sqrt(t8)
    t832 = t8 * t812
    t82 = t8*t8
    t83 = t82*t8
    t85 = t82*t83
    t86 = t85*t8
    t8m1 = 1.0/t8
    t8m2 = t8m1*t8m1
    t8m3 = t8m2*t8m1
    t8m5 = t8m3*t8m2
    t8m6 = t8m5*t8m1

    # The Fermi temperature, eq. 4.3.
    tfermi = 5.9302e9*(sqrt(1.0+1.018*(den6*ye)**twoth)-1.0)

    if temp > 0.01 * tfermi:
        # "Weak" degenerate electrons only; eta, eq. 5.3.
        dum = 7.05e6 * t832 + 5.12e4 * t83
        dumdt = (1.5*7.05e6*t812 + 3.0*5.12e4*t82)*1.0e-8

        z = 1.0/dum
        eta = rm*z
        etadt = -rm*z*z*dumdt
        etadd = rmdd*z
        etada = rmda*z
        etadz = rmdz*z

        etam1 = 1.0/eta
        etam2 = etam1 * etam1
        etam3 = etam2 * etam1

        # F, eq. 5.2.
        a0 = 23.5 + 6.83e4*t8m2 + 7.81e8*t8m5
        f0 = (-2.0*6.83e4*t8m3 - 5.0*7.81e8*t8m6)*1.0e-8
        xnum = 1.0/a0

        dum = 1.0 + 1.47*etam1 + 3.29e-2*etam2
        z = -1.47*etam2 - 2.0*3.29e-2*etam3
        dumdt = z*etadt
        dumdd = z*etadd
        dumda = z*etada
        dumdz = z*etadz

        c00 = 1.26 * (1.0+etam1)
        z = -1.26*etam2
        c01 = z*etadt
        c02 = z*etadd
        c03 = z*etada
        c04 = z*etadz

        z = 1.0/dum
        xden = c00*z
        xdendt = (c01 - xden*dumdt)*z
        xdendd = (c02 - xden*dumdd)*z
        xdenda = (c03 - xden*dumda)*z
        xdendz = (c04 - xden*dumdz)*z

        fbrem = xnum + xden
        fbremdt = -xnum*xnum*f0 + xdendt
        fbremdd = xdendd
        fbremda = xdenda
        fbremdz = xdendz

        # G, eq. 5.9.
        a0 = 230.0 + 6.7e5*t8m2 + 7.66e9*t8m5
        f0 = (-2.0*6.7e5*t8m3 - 5.0*7.66e9*t8m6)*1.0e-8

        z = 1.0 + rm*1.0e-9
        dum = a0*z
        dumdt = f0*z
        z = a0*1.0e-9
        dumdd = z*rmdd
        dumda = z*rmda
        dumdz = z*rmdz

        xnum = 1.0/dum
        z = -xnum*xnum
        xnumdt = z*dumdt
        xnumdd = z*dumdd
        xnumda = z*dumda
        xnumdz = z*dumdz

        c00 = 7.75e5*t832 + 247.0*t8**(3.85)
        dd00 = (1.5*7.75e5*t812 + 3.85*247.0*t8**(2.85))*1.0e-8

        c01 = 4.07 + 0.0240 * t8**(1.4)
        dd01 = 1.4*0.0240*t8**(0.4)*1.0e-8

        c02 = 4.59e-5 * t8**(-0.110)
        dd02 = -0.11*4.59e-5 * t8**(-1.11)*1.0e-8

        z = den**(0.656)
        dum = c00*rmi + c01 + c02*z
        dumdt = dd00*rmi + dd01 + dd02*z
        z = -c00*rmi*rmi
        dumdd = z*rmdd + 0.656*c02*den**(-0.454)
        dumda = z*rmda
        dumdz = z*rmdz

        xden = 1.0/dum
        z = -xden*xden
        xdendt = z*dumdt
        xdendd = z*dumdd
        xdenda = z*dumda
        xdendz = z*dumdz

        gbrem = xnum + xden
        gbremdt = xnumdt + xdendt
        gbremdd = xnumdd + xdendd
        gbremda = xnumda + xdenda
        gbremdz = xnumdz + xdendz

        # Eq. 5.1.
        dum = 0.5738*zbar*ye*t86*den
        dumdt = 0.5738*zbar*ye*6.0*t85*den*1.0e-8
        dumdd = 0.5738*zbar*ye*t86
        dumda = -dum*abari
        dumdz = 0.5738*2.0*ye*t86*den

        z = tfac4*fbrem - tfac5*gbrem
        sbrem = dum * z
        sbremdt = dumdt*z + dum*(tfac4*fbremdt - tfac5*gbremdt)
        sbremdd = dumdd*z + dum*(tfac4*fbremdd - tfac5*gbremdd)
        sbremda = dumda*z + dum*(tfac4*fbremda - tfac5*gbremda)
        sbremdz = dumdz*z + dum*(tfac4*fbremdz - tfac5*gbremdz)

    else:
        # Liquid metal, with C12 parameters (not too different for other
        # elements); u, eqs. 5.18 and 5.16.
        u = fac3 * (log10(den) - 3.0)
        a0 = iln10*fac3*deni

        cos1 = cos(u)
        cos2 = cos(2.0*u)
        cos3 = cos(3.0*u)
        cos4 = cos(4.0*u)
        cos5 = cos(5.0*u)

        sin1 = sin(u)
        sin2 = sin(2.0*u)
        sin3 = sin(3.0*u)
        sin4 = sin(4.0*u)
        sin5 = sin(5.0*u)

        # fb (eq. 5.21), with its density derivative c00.
        fb = (0.5 * 0.17946 + 0.00945*u + 0.34529
              - 0.05821*cos1 - 0.04969*sin1
              - 0.01089*cos2 - 0.01584*sin2
              - 0.01147*cos3 - 0.00504*sin3
              - 0.00656*cos4 - 0.00281*sin4
              - 0.00519*cos5)

        c00 = a0*(0.00945
                  + 0.05821*sin1 - 0.04969*cos1
                  + 0.01089*sin2*2.0 - 0.01584*cos2*2.0
                  + 0.01147*sin3*3.0 - 0.00504*cos3*3.0
                  + 0.00656*sin4*4.0 - 0.00281*cos4*4.0
                  + 0.00519*sin5*5.0)

        # ft (eq. 5.22), with its derivative c01.
        ft = (0.5 * 0.06781 - 0.02342*u + 0.24819
              - 0.00944*cos1 - 0.02213*sin1
              - 0.01289*cos2 - 0.01136*sin2
              - 0.00589*cos3 - 0.00467*sin3
              - 0.00404*cos4 - 0.00131*sin4
              - 0.00330*cos5)

        c01 = a0*(-0.02342
                  + 0.00944*sin1 - 0.02213*cos1
                  + 0.01289*sin2*2.0 - 0.01136*cos2*2.0
                  + 0.00589*sin3*3.0 - 0.00467*cos3*3.0
                  + 0.00404*sin4*4.0 - 0.00131*cos4*4.0
                  + 0.00330*sin5*5.0)

        # gb (eq. 5.23), with its derivative c02.
        gb = (0.5 * 0.00766 - 0.01259*u + 0.07917
              - 0.00710*cos1 + 0.02300*sin1
              - 0.00028*cos2 - 0.01078*sin2
              + 0.00232*cos3 + 0.00118*sin3
              + 0.00044*cos4 - 0.00089*sin4
              + 0.00158*cos5)

        c02 = a0*(-0.01259
                  + 0.00710*sin1 + 0.02300*cos1
                  + 0.00028*sin2*2.0 - 0.01078*cos2*2.0
                  - 0.00232*sin3*3.0 + 0.00118*cos3*3.0
                  - 0.00044*sin4*4.0 - 0.00089*cos4*4.0
                  - 0.00158*sin5*5.0)

        # gt (eq. 5.24), with its derivative c03.
        gt = (-0.5 * 0.00769 - 0.00829*u + 0.05211
              + 0.00356*cos1 + 0.01052*sin1
              - 0.00184*cos2 - 0.00354*sin2
              + 0.00146*cos3 - 0.00014*sin3
              + 0.00031*cos4 - 0.00018*sin4
              + 0.00069*cos5)

        c03 = a0*(-0.00829
                  - 0.00356*sin1 + 0.01052*cos1
                  + 0.00184*sin2*2.0 - 0.00354*cos2*2.0
                  - 0.00146*sin3*3.0 - 0.00014*cos3*3.0
                  - 0.00031*sin4*4.0 - 0.00018*cos4*4.0
                  - 0.00069*sin5*5.0)

        # dum and the powers of its inverse.
        dum = 2.275e-1 * zbar * zbar*t8m1 * (den6*abari)**oneth
        dumdt = -dum*tempi
        dumdd = oneth*dum*deni
        dumda = -oneth*dum*abari
        dumdz = 2.0*dum*zbari

        gm1 = 1.0/dum
        gm2 = gm1*gm1
        gm13 = gm1**oneth
        gm23 = gm13 * gm13
        gm43 = gm13*gm1
        gm53 = gm23*gm1

        # v and w, eqs. 5.25 and 5.26, with their derivatives a0, a1.
        v = -0.05483 - 0.01946*gm13 + 1.86310*gm23 - 0.78873*gm1
        a0 = oneth*0.01946*gm43 - twoth*1.86310*gm53 + 0.78873*gm2

        w = -0.06711 + 0.06859*gm13 + 1.74360*gm23 - 0.74498*gm1
        a1 = -oneth*0.06859*gm43 - twoth*1.74360*gm53 + 0.74498*gm2

        # Eqs. 5.19 and 5.20.
        fliq = v*fb + (1.0 - v)*ft
        fliqdt = a0*dumdt*(fb - ft)
        fliqdd = a0*dumdd*(fb - ft) + v*c00 + (1.0 - v)*c01
        fliqda = a0*dumda*(fb - ft)
        fliqdz = a0*dumdz*(fb - ft)

        gliq = w*gb + (1.0 - w)*gt
        gliqdt = a1*dumdt*(gb - gt)
        gliqdd = a1*dumdd*(gb - gt) + w*c02 + (1.0 - w)*c03
        gliqda = a1*dumda*(gb - gt)
        gliqdz = a1*dumdz*(gb - gt)

        # Eq. 5.17.
        dum = 0.5738*zbar*ye*t86*den
        dumdt = 0.5738*zbar*ye*6.0*t85*den*1.0e-8
        dumdd = 0.5738*zbar*ye*t86
        dumda = -dum*abari
        dumdz = 0.5738*2.0*ye*t86*den

        z = tfac4*fliq - tfac5*gliq
        sbrem = dum * z
        sbremdt = dumdt*z + dum*(tfac4*fliqdt - tfac5*gliqdt)
        sbremdd = dumdd*z + dum*(tfac4*fliqdd - tfac5*gliqdd)
        sbremda = dumda*z + dum*(tfac4*fliqda - tfac5*gliqda)
        sbremdz = dumdz*z + dum*(tfac4*fliqdz - tfac5*gliqdz)

    # ---- Recombination neutrinos ----
    # Eq. 6.11, solved for nu.
    xnum = 1.10520e8 * den * ye / (temp*sqrt(temp))
    xnumdt = -1.50*xnum*tempi
    xnumdd = xnum*deni
    xnumda = -xnum*abari
    xnumdz = xnum*zbari

    # The chemical potential nu ...
    nu = ifermi12(xnum)

    # ... and its derivatives, through the integral of order -1/2.
    a0 = 1.0/(0.5*zfermim12(nu))
    nudt = a0*xnumdt
    nudd = a0*xnumdd
    nuda = a0*xnumda
    nudz = a0*xnumdz

    nu2 = nu * nu
    nu3 = nu2 * nu

    # Table 12.
    if nu >= -20.0 and nu < 0.0:
        a1 = 1.51e-2
        a2 = 2.42e-1
        a3 = 1.21
        b = 3.71e-2
        c = 9.06e-1
        d = 9.28e-1
        f1 = 0.0
        f2 = 0.0
        f3 = 0.0
    elif nu >= 0.0 and nu <= 10.0:
        a1 = 1.23e-2
        a2 = 2.66e-1
        a3 = 1.30
        b = 1.17e-1
        c = 8.97e-1
        d = 1.77e-1
        f1 = -1.20e-2
        f2 = 2.29e-2
        f3 = -1.04e-3

    # Question (see the page): zero outside -20 <= nu <= 10.
    if nu >= -20.0 and nu <= 10.0:
        # Eqs. 6.7, 6.13 and 6.14.
        zeta = 1.579e5*zbar*zbar*tempi
        zetadt = -zeta*tempi
        zetadd = 0.0
        zetada = 0.0
        zetadz = 2.0*zeta*zbari

        c00 = 1.0/(1.0 + f1*nu + f2*nu2 + f3*nu3)
        c01 = f1 + f2*2.0*nu + f3*3.0*nu2
        dum = zeta*c00
        dumdt = zetadt*c00 + zeta*c01*nudt
        dumdd = zeta*c01*nudd
        dumda = zeta*c01*nuda
        dumdz = zetadz*c00 + zeta*c01*nudz

        z = 1.0/dum
        dd00 = dum**(-2.25)
        dd01 = dum**(-4.55)
        c00 = a1*z + a2*dd00 + a3*dd01
        c01 = -(a1*z + 2.25*a2*dd00 + 4.55*a3*dd01)*z

        z = exp(c*nu)
        dd00 = b*z*(1.0 + d*dum)
        gum = 1.0 + dd00
        gumdt = dd00*c*nudt + b*z*d*dumdt
        gumdd = dd00*c*nudd + b*z*d*dumdd
        gumda = dd00*c*nuda + b*z*d*dumda
        gumdz = dd00*c*nudz + b*z*d*dumdz

        # Eq. 6.5.
        z = exp(nu)
        a1 = 1.0/gum

        bigj = c00 * z * a1
        bigjdt = c01*dumdt*z*a1 + c00*z*nudt*a1 - c00*z*a1*a1 * gumdt
        bigjdd = c01*dumdd*z*a1 + c00*z*nudd*a1 - c00*z*a1*a1 * gumdd
        bigjda = c01*dumda*z*a1 + c00*z*nuda*a1 - c00*z*a1*a1 * gumda
        bigjdz = c01*dumdz*z*a1 + c00*z*nudz*a1 - c00*z*a1*a1 * gumdz

        # Eq. 6.5.
        z = exp(zeta + nu)
        dum = 1.0 + z
        a1 = 1.0/dum
        a2 = 1.0/bigj

        sreco = tfac6 * 2.649e-18 * ye * zbar**13 * den * bigj*a1
        srecodt = sreco*(bigjdt*a2 - z*(zetadt + nudt)*a1)
        srecodd = sreco*(1.0*deni + bigjdd*a2 - z*(zetadd + nudd)*a1)
        srecoda = sreco*(-1.0*abari + bigjda*a2 - z*(zetada+nuda)*a1)
        srecodz = sreco*(14.0*zbari + bigjdz*a2 - z*(zetadz+nudz)*a1)

    # Convert from erg/cm^3/s to erg/g/s (comment these out to duplicate
    # the Itoh et al. plots).
    spair = spair*deni
    spairdt = spairdt*deni
    spairdd = spairdd*deni - spair*deni
    spairda = spairda*deni
    spairdz = spairdz*deni

    splas = splas*deni
    splasdt = splasdt*deni
    splasdd = splasdd*deni - splas*deni
    splasda = splasda*deni
    splasdz = splasdz*deni

    sphot = sphot*deni
    sphotdt = sphotdt*deni
    sphotdd = sphotdd*deni - sphot*deni
    sphotda = sphotda*deni
    sphotdz = sphotdz*deni

    sbrem = sbrem*deni
    sbremdt = sbremdt*deni
    sbremdd = sbremdd*deni - sbrem*deni
    sbremda = sbremda*deni
    sbremdz = sbremdz*deni

    sreco = sreco*deni
    srecodt = srecodt*deni
    srecodd = srecodd*deni - sreco*deni
    srecoda = srecoda*deni
    srecodz = srecodz*deni

    # The total and its derivatives.
    snu = splas + spair + sphot + sbrem + sreco
    dsnudt = splasdt + spairdt + sphotdt + sbremdt + srecodt
    dsnudd = splasdd + spairdd + sphotdd + sbremdd + srecodd
    dsnuda = splasda + spairda + sphotda + sbremda + srecoda
    dsnudz = splasdz + spairdz + sphotdz + sbremdz + srecodz

    return snu, dsnudt, dsnudd, dsnuda, dsnudz
