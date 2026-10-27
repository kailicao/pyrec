"""Python version of wrtmil.f (YREC 09c3b9b): the Milne (homology)
invariants of one model, to the Milne output file."""

import math

from pyrec import comp_mod
from pyrec import constants_mod
from pyrec import iounits_mod
from pyrec import nmlpar_mod
from pyrec import scrtch_mod


def _row(i, hcomp, hd, hl, hp, hr, hs1):
    """One shell's line (FORMAT 10): U, V, W and n + 1 = 1/nabla."""
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


def wrtmil(hcomp, hd, hl, hp, hr, hs1, m, model):
    """Write the model's header and U, V, W and n + 1 at the first shell,
    every nprtpt-th shell and the last one, then close the unit.

    Translates: wrtmil.f#WRTMIL (09c3b9b)

    Arguments: the model hcomp[species, shell], log rho hd, L/Lsun hl,
    log P hp, log r hr, m in grams hs1 (index 0 unused); m, the number of
    shells; model, the model's number. None is returned.
    """
    c = constants_mod
    io = iounits_mod
    nprtpt = nmlpar_mod.nprtpt
    smtot = math.exp(c.cln * comp_mod.stotal) / c.cmsun
    io.write(io.imilne,
             f'{"":10}MODEL{model:5d}  XENV ={comp_mod.xenv:10.3E}'
             f'  ZENV ={comp_mod.zenv:10.3E}  MASS(SOLAR UNITS) ={smtot:10.3E}')
    io.write(io.imilne, _row(1, hcomp, hd, hl, hp, hr, hs1))
    iend = 1
    if nprtpt <= m:
        ibeg = max(2, nprtpt)
        iend = m - m % nprtpt
        for i in range(ibeg, iend + 1, nprtpt):     # loop 20
            io.write(io.imilne, _row(i, hcomp, hd, hl, hp, hr, hs1))
    if iend < m:
        io.write(io.imilne, _row(m, hcomp, hd, hl, hp, hr, hs1))
    # Question (see the page): closed after each model, opened only once.
    io.close(io.imilne)
