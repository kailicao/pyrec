# PyREC: the Yale Rotating Evolution Code in Python

**This repository is under development and contains unvetted AI-generated code. It is not part of YREC; please refer to the [YREC documentation](https://yrec.readthedocs.io/en/latest/).**

PyREC (the package `pyrec`) is a Python version of [YREC](https://github.com/yreclab/yrec), the Yale Rotating Evolution Code, written by Kaili Cao with AI (Claude) as scaffolding for overhauling the code. YREC's developers have not vetted it.

## What it is

A translation that keeps YREC's current structure: one Python file per Fortran file (`pyrec/amcalc.py` for `src/amcalc.f`) and one function per routine, made from YREC's development branch at 09c3b9b (Oct 20, 2026), brought to a8a5950 (Oct 30, 2026) in the eight source files changed between the two, and to 7be4133 (Nov 2, 2026) in the 67 changed since (two unused modules removed). It follows YREC at the commit Kaili Cao chooses (7be4133 since Nov 9, 2026), moved only when the work needs a later one. Its structure will change only once the whole translation has been checked; the new structure planned for it is drafted, with Python stubs, on its [documentation site](https://github.com/kailicao/pyrec-map).

## Status

Every one of YREC's 282 source files is translated (at 09c3b9b, Oct 27, 2026; the files changed since updated to a8a5950 and 7be4133), and the calls agree with what each routine returns and with its arguments (checked by script). Nothing has been run: The translation is checked only by compiling and linting it, and is to be tested routine by routine against the Fortran.

## Conventions (provisional)

The choices made for the translation, so that each Python routine can be read, and later tested, beside its Fortran. Examples: `hpoint.py`, `mdot.py`, `esac06.py`.

**Structure**

- **One function per Fortran unit,** with the unit's name and arguments in lower case, in the same order. Arrays are changed in place, as in Fortran. Scalars a routine sets are returned, as a tuple in argument order, and are also passed in (with any value), so that argument lists match the Fortran's.
- **A FUNCTION that also sets arguments** returns `(value, <set arguments in order>)`. **An alternate return** (`CALL X(…, *999)`) becomes a boolean returned last, `alt1` (true when the Fortran would take it; `alt2` for a second).
- **Modules:** a Fortran module becomes a Python module of the same name (`comp_mod.f` becomes `comp_mod.py`), its variables module-level names in lower case. A routine imports the module (`from pyrec import comp_mod`) and writes `comp_mod.xenv`, so that each name's origin stays visible.
- **Calls** import the routine from its file: `from pyrec.physic import physic`.
- **Locals kept between calls** (`SAVE`, or locals initialized by `DATA`) are module-level variables of the routine's file. An array declared with one element and used as a scalar is a Python scalar.
- **Names:** YREC's, in lower case; what each means is in the site's glossary.

**Values**

- **Indices keep Fortran's values:** arrays have one extra element, so that index 0 is unused, and keep Fortran's index order (`hcomp[k, j]` for `HCOMP(K, J)`).
- **Bug for bug:** the Fortran is reproduced exactly, its apparent bugs included. A constant the Fortran writes in single precision keeps its single-precision value (`float(np.float32(1e-2))`), and sums keep the Fortran's order (an explicit loop, not `np.sum`), so that results can match to the last bit in differential tests.

**Control flow and I/O**

- **GOTOs** become `while` loops, `break`, `continue`, early returns, or flags, with the statement label in a comment where it helps (`# label 20`). A few private helpers (a leading underscore) take blocks the Fortran repeats.
- **`STOP`** becomes `raise SystemExit(<exit status>)`, after the message.
- **Fortran I/O** goes through `iounits_mod`: `write(unit, text)`, `open(unit, name, status=...)`, `close(unit)`, `rewind(unit)`, and `read(unit)` (one record, a string, parsed as the FORMAT says; `None` at end of file, for `END=`). FORMATs become f-strings.

**Reading it**

- **Docstrings:** the routine's purpose in a sentence or two, then `Translates: <file>#<UNIT> (<commit>)`, the YREC commit it follows, then the arguments and what is returned. A blank line before a multi-line docstring's closing quotes, and one after the docstring.
- **Type annotations** on every routine's arguments and returns (`float`, `int`, `bool`, `str`, `np.ndarray`), from YREC's declarations.
- **Blocks:** a blank line between logical blocks, and each block's opening comment ends with the Fortran lines it translates: `# Mass change in g (Fortran lines 210–212).`
- **The site's points:** `# Suspected bug (see the page): …` where the page lists a suspected bug, `# Question (see the page): …` where it asks a question.
- **Checked** by compiling, by ruff (`pyproject.toml`: one-line statements and YREC's names allowed; a local the Fortran also never reads is kept, marked `# noqa: F841`), and by pyrec-map's `scripts/check_returns.py` (each call's unpacking and arguments against the routine it calls). Nothing is run until the differential tests.

## Where the Python differs from the Fortran

On six paths where the Fortran has a suspected bug, its behavior is undefined (a read or write outside an array) or IEEE-specific (NaN, overflow), and the Python does not yet reproduce it; the site's page of each file gives the bug.

- `wtime.py` raises `ZeroDivisionError` where the Fortran gives NaN (wtime.f, bug 1).
- `hpoint.py` sizes `nflag` with 102 elements, so it stays in bounds where the Fortran writes past `NFLAG(101)` (hpoint.f, bug 2).
- `liburn.py` and `liburn2.py` raise `OverflowError` where the Fortran's `EXP` overflows (liburn.f, bug 3; liburn2.f, bug 5).
- `liburn2.py` reads `hro[0]` as 0 where the Fortran reads outside its array (liburn2.f, bug 4).
- `microdiff_run.py` prints 0 where the Fortran prints values left from an earlier call.
- `meval.py` reads `xval[0]`, and would raise `IndexError` at `xval[nume+1]`, where the Fortran reads outside its array (meval.f, bug 1).

Proposed: NaN and infinities reproduced with NumPy floats, where IEEE arithmetic defines them; the undefined reads marked as such, since no translation can match them.
