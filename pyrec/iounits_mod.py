"""Python version of iounits_mod.f (YREC 09c3b9b): the logical unit numbers
of every file, formerly COMMON/IOUNITS/ (and before it /LUOUT/, /LUNUM/
and eight other blocks).

Translates: iounits_mod.f#IOUNITS_MOD (09c3b9b)

Module-level names, as in the Fortran module, in lower case (all INTEGER*4;
parmin assigns them). Routines that USE IOUNITS_MOD import this module and
write ``iounits_mod.<name>``.

The functions below are not in the Fortran module: they stand in for
Fortran's own I/O statements (OPEN, CLOSE, REWIND, READ, WRITE on a unit),
so that every translated routine does its I/O the same way. Records are
lines of text; FORMATs are applied by the callers.
"""

import builtins
import os
import sys

ilast = 0           # last model (11)
idebug = 0          # debugging output (18)
itrack = 0          # track (19)
ishort = 0          # diagnostic output, the .short file (20)
imilne = 0          # Milne invariant variables (21)
imodpt = 0          # shell-by-shell output (22)
istor = 0           # saved models (23)
iowr = 0            # status output, the screen (6; main sets 9 first)
ifirst = 0          # starting model (12)
irun = 0            # PHYSICS namelist (13)
istand = 0          # CONTROL namelist (14)
ifermi = 0          # Fermi-Dirac tables (15)
iopmod = 0          # pulsation output, interior (24)
iopenv = 0          # pulsation output, envelope (25)
iopatm = 0          # pulsation output, atmosphere (26)
isnu = 0            # neutrino fluxes (33)
ialex06 = 0         # Alexander 2006 opacities (90)
iclcd = 0           # acoustic-depth output (91, only with lacout)
ijlast = 0          # used by wrtout (ljlast), assigned nowhere
iopurez = 0         # pure-CN opacities (62)
icondopacp = 0      # Potekhin conductive opacities (75); IcondOpacP
iliv95 = 0          # OPAL95 opacities (48)
ioatm = 0           # Kurucz atmospheres (38)
ioatma = 0          # Allard atmospheres (66)
iopale = 0          # OPAL 2006 EOS (49)
iscvh = 0           # SCV EOS, hydrogen (72)
iscvhe = 0          # SCV EOS, helium (73)
iscvz = 0           # SCV EOS, metals (74)
isound = 0          # sound-speed profile (500)


# --- Stand-ins for Fortran I/O (Python only) -------------------------------

_units = {}         # unit number -> open text file


def open(unit, name, status='UNKNOWN', form='FORMATTED', access='SEQUENTIAL'):
    """OPEN(unit, FILE=name, STATUS=status): connect a unit to a file.

    status is Fortran's: 'OLD' (must exist), 'NEW' (must not exist),
    'REPLACE', 'SCRATCH' or 'UNKNOWN' (kept if it exists, positioned at its
    start). form and access are accepted and ignored (formatted, sequential).
    """
    close(unit)
    name = str(name).strip()
    status = status.upper()
    if status == 'OLD':
        mode = 'r+' if os.access(name, os.W_OK) else 'r'
    elif status == 'NEW':
        mode = 'x+'
    elif status in ('REPLACE', 'SCRATCH'):
        mode = 'w+'
    else:   # UNKNOWN
        mode = 'r+' if os.path.exists(name) else 'w+'
    _units[unit] = builtins.open(name, mode)


def _file(unit):
    """The file of a unit, connecting an unopened unit to fort.<unit> as
    Fortran does."""
    if unit not in _units:
        open(unit, f'fort.{unit}')
    return _units[unit]


def close(unit):
    """CLOSE(unit)."""
    f = _units.pop(unit, None)
    if f is not None:
        f.close()


def rewind(unit):
    """REWIND(unit)."""
    if unit in _units:
        _units[unit].seek(0)


def read(unit):
    """READ(unit, ...): the next record as a string, without its newline,
    or None at the end of the file (Fortran's END= branch). Unit 5 is
    standard input."""
    if unit == 5:
        line = sys.stdin.readline()
    else:
        line = _file(unit).readline()
    if line == '':
        return None
    return line.rstrip('\n')


def write(unit, text):
    """WRITE(unit, ...): write text as one or more records. Units 6 and '*'
    are the screen. As in a sequential file, what follows the record just
    written is cut off."""
    if unit in (6, '*') and unit not in _units:
        print(text)
        return
    f = _file(unit)
    f.write(f'{text}\n')
    f.truncate()
    f.flush()
