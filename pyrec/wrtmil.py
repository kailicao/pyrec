"""Python version of wrtmil.f (YREC 7be4133): the Milne (homology)
invariants of one model, to the Milne output file."""

import math

import numpy as np

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod


def _row(i: int, hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray,
         hs1: np.ndarray) -> str:
    """One shell's line (FORMAT 10): U, V, W and n + 1 = 1/nabla (the
    Fortran repeats these lines at each shell: lines 71–79, 86–94, and
    99–107).

    """

    c = constants_mod
    s = scrtch_mod
    d = math.exp(c.cln * hd[i])
    p = math.exp(c.cln * hp[i])
    r = math.exp(c.cln * hr[i])
    u = c.c4pi * d * r**3 / hs1[i]
    v = math.exp(c.cln * c.cgl) * hs1[i] * d / (p * r)
    # Question (see the page): W includes the gravitational term seg[7].
    w = u * hs1[i] * (s.sesum[i] + s.seg[7, i]) / (hl[i] * c.clsun)
    np1 = 1.0 / s.sdel[2, i]
    vals = (hs1[i], r, p, d, hcomp[1, i], hl[i], u, v, w, np1)
    return f' {i:4d}' + ''.join(f'{x:11.3E}' for x in vals)


def _d(x: float) -> str:
    """Fortran's 1PD10.3: like 1PE10.3, with a D exponent.

    """

    return f'{x:10.3E}'.replace('E', 'D')


def wrtmil(hcomp: np.ndarray, hd: np.ndarray, hl: np.ndarray, hp: np.ndarray, hr: np.ndarray, hs1: np.ndarray, m: int,
           model: int) -> None:
    """Write the model's header and U, V, W and n + 1 at the first shell,
    every nprtpt-th shell and the last one, then close the unit.

    Translates: wrtmil.f#WRTMIL (7be4133)

    Arguments: the model hcomp[species, shell], log rho hd, L/Lsun hl,
    log P hp, log r hr, m in grams hs1 (index 0 unused); m, the number of
    shells; model, the model's number. None is returned.

    """

    c = constants_mod
    io = iounits_mod
    nprtpt = nmlpar_mod.nprtpt

    # The header; Fortran's 1PD10.3 writes a D exponent (Fortran lines
    # 66–69).
    smtot = math.exp(c.cln * comp_mod.stotal) / c.cmsun
    io.write(io.imilne,
             f'{"":10}MODEL{model:5d}  XENV ={_d(comp_mod.xenv)}'
             f'  ZENV ={_d(comp_mod.zenv)}  MASS(SOLAR UNITS) ={_d(smtot)}')

    # The first shell, always written (Fortran lines 70–79).
    io.write(io.imilne, _row(1, hcomp, hd, hl, hp, hr, hs1))

    # Every nprtpt-th shell (Fortran lines 81–96).
    iend = 1
    if nprtpt <= m:
        ibeg = max(2, nprtpt)
        iend = m - m % nprtpt
        for i in range(ibeg, iend + 1, nprtpt):     # loop 20
            io.write(io.imilne, _row(i, hcomp, hd, hl, hp, hr, hs1))

    # The last shell, if nprtpt does not divide m (Fortran lines 97–108).
    if iend < m:
        io.write(io.imilne, _row(m, hcomp, hd, hl, hp, hr, hs1))

    # Suspected bug (see the page): closed after each model, opened only
    # once (Fortran line 109).
    io.close(io.imilne)
