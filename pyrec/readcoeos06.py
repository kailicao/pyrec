"""Python version of readcoeos06.f (YREC 09c3b9b): reads the OPAL 2006
equation-of-state tables."""

import numpy as np

from pyrec import aeos06_mod
from pyrec import eeeos06_mod
from pyrec import eos06misc_mod
from pyrec import iounits_mod

MX = 5
MV = 10
NR = 169
NT = 197

# Locals kept between calls (SAVE), formerly in part COMMON/EEEEOS06/.
moles = np.zeros(MX + 1)
xin = np.zeros(MX + 1)
tmass = np.zeros(MX + 1)
amu_m = np.zeros((NR + 1, NT + 1))
alogne = np.zeros((NR + 1, NT + 1))
rhogr = np.zeros((MX + 1, NR + 1))
frac = np.zeros((MX + 1, 7))
alogr = np.zeros((NR + 1, NT + 1))
icycuse = np.zeros((MX + 1, NR + 1), dtype=int)
itimeco = 0
j = 0     # saved: the leftover index used in alogne[jcs, j] below


def _fields(line, widths, start=0, ints=()):
    """Split a fixed-format record into numbers (a blank field reads as 0,
    as in Fortran); the fields numbered in ints are integers."""
    out, pos = [], start
    for n, w in enumerate(widths):
        text = line[pos:pos + w].strip()
        if n in ints:
            out.append(int(text) if text else 0)
        else:
            out.append(float(text.replace('D', 'E')) if text else 0.0)
        pos += w
    return out


def readcoeos06():
    """Read the OPAL 2006 EOS tables (five hydrogen abundances) from unit
    IOPALE into aeos06_mod (xz, t6list, and the grids t6a, rho with their
    inverse spacings dfs, dfsr, dfsx) and eos06misc_mod (zz); on the
    first call, fill xz with 1e35 and set eos06misc_mod.eos to 1.

    Translates: readcoeos06.f#READCOEOS06 (09c3b9b)

    No arguments; nothing is returned.
    """
    global itimeco, j
    a = aeos06_mod
    e = eos06misc_mod
    iopale = iounits_mod.iopale
    ishort = iounits_mod.ishort

    if itimeco != 12345678:
        for i in range(1, MX + 1):
            for j in range(1, MV + 1):
                for k in range(1, NT + 1):
                    for l in range(1, NR + 1):
                        a.xz[i, j, k, l] = 1.0e35
        for j in range(1, MV + 1):
            e.eos[j] = 1.0
        j = MV + 1          # the value a Fortran DO leaves
        itimeco = 12345678

    # Question (see the page): unit 2 is closed though never opened here.
    iounits_mod.close(2)
    # Read the tables (the file is opened in parmin).
    for midx in range(1, MX + 1):
        a.midx = midx
        xin[midx], e.zz[midx], moles[midx], tmass[midx] = _fields(
            iounits_mod.read(iopale), [6, 3, 12, 11, 10, 17, 10],
            start=3)[0::2]
        frac[midx, 1:7] = _fields(
            iounits_mod.read(iopale),
            [14, 4, 14, 3, 11, 3, 11, 3, 11, 4, 11], start=21)[0::2]
        iounits_mod.read(iopale)
        for jcs in range(1, NR + 1):
            numtot, icycuse[midx, jcs], dum, dum, _, rhogr[midx, jcs] = \
                _fields(iounits_mod.read(iopale), [5, 5, 12, 12, 17, 15],
                        ints=(0, 1))
            if numtot != jcs:
                iounits_mod.write(
                    ishort, f" OEOS06 Data file incorrect: numtot,jcs= "
                    f"{numtot:5d}{jcs:5d}")
                raise SystemExit
            iounits_mod.read(iopale)
            iounits_mod.read(iopale)
            if icycuse[midx, jcs] < e.nta[jcs]:
                # Question (see the page): "OEOS96" for the 2006 tables.
                iounits_mod.write(
                    ishort, f"Problem with OEOS96 data files: X="
                    f"{xin[midx]:6.4f} density={rhogr[midx, jcs]:14.4E}")
                raise SystemExit
            for i in range(1, icycuse[midx, jcs] + 1):
                if i > e.nta[jcs]:
                    iounits_mod.read(iopale)
                    continue                      # go to 4
                vals = _fields(iounits_mod.read(iopale),
                               [11, 1, 6, 11, 13, 13, 11, 11,
                                10, 10, 10, 10, 10])
                a.t6list[jcs, i] = vals[0]
                amu_m[jcs, i] = vals[2]
                # Question (see the page): j, not i (j is left over).
                alogne[jcs, j] = vals[3]
                for iv in range(1, 10):
                    a.xz[midx, e.index[iv], i, jcs] = vals[3 + iv]
            iounits_mod.read(iopale)
            iounits_mod.read(iopale)
            iounits_mod.read(iopale)
        iounits_mod.read(iopale)

    for i in range(1, NT + 1):
        if a.t6list[1, i] == 0.0:
            iounits_mod.write(
                ishort, f"READCOEOS06: Error:{i:4d}-th T6 value is zero")
            raise SystemExit
        a.t6a[i] = a.t6list[1, i]
    for i in range(2, NT + 1):
        a.dfs[i] = 1.0 / (a.t6a[i] - a.t6a[i - 1])
    a.rho[1] = rhogr[1, 1]
    for i in range(2, NR + 1):
        a.rho[i] = rhogr[1, i]
        a.dfsr[i] = 1.0 / (a.rho[i] - a.rho[i - 1])
    for i in range(2, MX + 1):
        a.dfsx[i] = 1.0 / (eeeos06_mod.xx[i] - eeeos06_mod.xx[i - 1])
    iounits_mod.close(iopale)
