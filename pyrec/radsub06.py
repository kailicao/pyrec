"""Python version of radsub06.f (YREC 09c3b9b): radiation added to the
OPAL 2006 gas equation of state."""

from pyrec import eos06misc_mod

# DATA: sigma/c = a/4 in Mb per unit T6^4, and the gas constant in
# Mb cc per mole and unit T6.
SIGMACC = 1.8914785e-3
APROP = 83.14510


def radsub06(irad, t6, density, moles, tmass):
    """Rewrite EOS(IRI(1..7)) in EOS06MISC_MOD with radiation added (when
    irad != 0): pressure, energy, entropy, dE/drho, c_V, chi_rho, chi_T,
    and add the radiation's differences to Gamma_1 and Gamma_2/(Gamma_2-1).
    With irad = 0 only the gas values are rewritten (c_V times Nk/M).

    Translates: radsub06.f#RADSUB06 (09c3b9b)

    Arguments: irad; t6, T in 10^6 K; density (g/cc); moles; tmass.
    Returns nothing.
    """
    eos = eos06misc_mod.eos
    iri = eos06misc_mod.iri
    rat = SIGMACC
    molenak = moles * APROP     # Mb-cc/unit T6

    pt = eos[iri[1]]
    et = eos[iri[2]]
    st = eos[iri[3]]
    dedrho = eos[iri[4]]
    chir = eos[iri[6]] * eos[iri[1]] / pt
    chitt = (eos[iri[1]] * eos[iri[7]]) / pt
    cvtt = (eos[iri[5]] * molenak / tmass)
    gam3pt_norad = pt * chitt / (cvtt * density * t6)
    gam1t_norad = chir + chitt * gam3pt_norad
    gam2pt_norad = gam1t_norad / gam3pt_norad

    if irad != 0:
        pr = 4.0 / 3.0 * rat * t6**4    # Mb
        er = 3.0 * pr / density         # Mb-cc/gm
        sr = 4.0 / 3.0 * er / t6        # Mb-cc/(gm-unit T6)
        pt = eos[iri[1]] + pr
        et = eos[iri[2]] + er
        st = eos[iri[3]] + sr
        dedrho = eos[iri[4]] - er / density
        chir = eos[iri[6]] * eos[iri[1]] / pt
        chitt = (eos[iri[1]] * eos[iri[7]] + 4.0 * pr) / pt
        cvtt = (eos[iri[5]] * molenak / tmass + 4.0 * er / t6)
        gam3pt = pt * chitt / (cvtt * density * t6)
        gam1t = chir + chitt * gam3pt   # eq 16.16 Landau & Lifshitz
        gam2pt = gam1t / gam3pt

        eos[iri[8]] = eos[iri[8]] + gam1t - gam1t_norad
        eos[iri[9]] = eos[iri[9]] + gam2pt - gam2pt_norad
        # Question (see the page): Gamma_3 (IRI(10)) is not corrected.
    eos[iri[1]] = pt
    eos[iri[2]] = et
    eos[iri[3]] = st
    eos[iri[4]] = dedrho
    eos[iri[5]] = cvtt
    eos[iri[6]] = chir
    eos[iri[7]] = chitt
