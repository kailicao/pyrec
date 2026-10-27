# PyREC: the Yale Rotating Evolution Code in Python

**This repository is under development and contains unvetted AI-generated code. It is not part of YREC; please refer to the [YREC documentation](https://yrec.readthedocs.io/en/latest/).**

PyREC (the package `pyrec`) is a Python version of [YREC](https://github.com/yreclab/yrec), the Yale Rotating Evolution Code, written by Kaili Cao with AI (Claude) as scaffolding for overhauling the code. YREC's developers have not vetted it.

## What it is

A translation that keeps YREC's current structure: one Python file per Fortran file (`pyrec/amcalc.py` for `src/amcalc.f`) and one function per routine, made from YREC's development branch at 09c3b9b (Oct 20, 2026). It follows YREC as YREC changes. Its structure will change only once the whole translation has been checked; the new structure planned for it is drafted, with Python stubs, on its [documentation site](https://github.com/kailicao/pyrec-map).

## Status

Every one of YREC's 284 source files at 09c3b9b is translated (Oct 27, 2026), and the calls agree with what each routine returns (checked by script). Nothing has been run: the translation is checked only by compiling it, and is to be tested routine by routine against the Fortran.

## Conventions (provisional)

- Each function has its routine's name and arguments, in lower case and in the same order. Arrays are changed in place; scalars a routine sets are returned, in argument order.
- A Fortran module becomes a Python module of the same name (`comp_mod.f` becomes `comp_mod.py`), and routines refer to its variables as `comp_mod.xenv`.
- Indices keep Fortran's values: arrays have one extra element, so that index 0 is unused, and keep Fortran's index order.
- Plain Python and NumPy; the Fortran's apparent bugs are reproduced and marked with comments.
