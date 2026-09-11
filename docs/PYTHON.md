# Python

**Yes, and it costs a loader.** A REBOUND user in Python reaches a
custom C integrator by name, with no binding layer, no fork of
REBOUND's Python package and no change to it. This file is what was
measured, where the one obstacle is, and what a user installs.

    import rebound
    import cft_rebound                     # python/cft_rebound.py

    cft_rebound.load("build/libcft_ias15.so")
    sim = rebound.Simulation()
    sim.integrator = "ias15_cft"           # works: the name goes straight to C

Investigated 2026-09-10 against REBOUND 5.1.1 - the pinned clone at
bdfda4bd and, separately, the official PyPI wheels, which are at
githash 33549d1d. Extended to Windows 2026-09-11, against the
`cp312-win_amd64` wheel of that same 5.1.1; the loading rule there is
the opposite of the POSIX one and has its own section. The numbers and
the failures are in docs/VALIDATION.md.

## Why it works

REBOUND's Python layer keeps no list of integrators. `Simulation.__setattr__`
(`rebound/simulation.py`) intercepts `integrator`, lowercases the
string, takes five WHFast/SABA shortcuts, and otherwise passes the name
through untouched:

    clibrebound.reb_simulation_set_integrator.argtypes = [POINTER(Simulation), c_char_p]
    clibrebound.reb_simulation_set_integrator(byref(self), c_char_p(value.encode("ascii")))
    self.process_messages()

`reb_simulation_set_integrator` (REBOUND's own src/simulation.c, in the
pinned clone) tries the built-ins and then walks
`reb_integrator_configurations_custom`, the global that
`reb_integrator_register` writes. So the question was never about
Python. It was whether a separately loaded library's registration lands
in **that** global - the one belonging to the librebound the Python
package already has open.

It does, and then everything else in REBOUND's Python layer works with
no special-casing at all, because that layer drives the integrator
entirely through the `struct reb_integrator` you registered. Measured
with a probe integrator:

- `sim.integrator` reads back the registered name;
- `sim.integrate()` and `sim.steps()` call your `step`;
- `sim.integrator.<field>` reads and writes your state through your
  `field_descriptor_list`, by the name the descriptor gives - so for
  REBOUND's own ias15, `sim.integrator.epsilon = 0.0` is how a Python
  user turns off adaptive stepping. For this integrator the same
  mechanism works and the **name is different**: every entry in
  `cft_ias15_field_descriptor_list` is `cft_`-prefixed, so
  `sim.integrator.epsilon` raises `AttributeError: Field 'epsilon' not
  found` and `sim.integrator.cft_epsilon` resolves. Measured against
  the 5.1.1 wheel: setting `cft_format`, `cft_epsilon` and
  `cft_max_iter` that way gives a run identical to
  `cft_rebound.configure()`'s (docs/VALIDATION.md entry 29). Use
  whichever you prefer; `configure()` validates the values and applies
  the per-format `max_iter` default, and the attribute path is there
  when you want REBOUND's own idiom;
- `repr(sim.integrator)` prints every field, and `sim.integrator.__doc__`
  is generated from your `documentation` string and your field
  documentation;
- `sim.status()` names it;
- `sim.save_to_file()` writes your fields into the Simulationarchive as
  `integrator.<name>.<field>`, and a process that has the library loaded
  restores them.

## The one obstacle, and it is POSIX-only

`rebound/__init__.py` loads the library with `cdll.LoadLibrary`, which
is `RTLD_LOCAL`. Its symbols are therefore not in the global scope, and
the obvious snippet fails:

    >>> import rebound, ctypes
    >>> ctypes.CDLL("libcft_ias15.so")
    OSError: libcft_ias15.so: undefined symbol: reb_integrator_register

The fix is one line. Re-opening the same file with `RTLD_GLOBAL`
promotes the mapping already in the process - it is not a second copy;
the address of `reb_integrator_configurations_custom` seen through both
handles is the same - and our library then resolves against it:

    ctypes.CDLL(rebound.__libpath__, mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL("libcft_ias15.so")        # registers "ias15_cft"

That is `cft_rebound.load()`, and it is the entire Python-side cost.

**On Windows this particular obstacle does not arise** - and a
different one does. See "The Windows rule" below; it works now, and it
is a different rule rather than the same one with a flag.

## What a user installs, in what order

1. `pip install rebound`. The wheel puts **both** halves you need side
   by side in site-packages: `librebound.<abi>.so` (or `.pyd`), which is
   `rebound.__libpath__`, and a `src/` directory with all of REBOUND's
   headers, which is `cft_rebound.include_dir()`. No REBOUND source
   checkout is needed.
2. Build this repository's integrator as a **shared library**, compiled
   against those headers - not against a different REBOUND's, because
   `struct reb_simulation`'s layout must match the library that will be
   calling you. On POSIX link against nothing; on Windows link against
   `rebound.__libpath__`.
3. `python python/cft_rebound.py path/to/the/library` to check: it
   prints the REBOUND it found, the header directory, on Windows the
   librebound the library was linked against, and the integrator names
   now registered.
4. `python python/example_equivalence.py --library path/to/the/library`
   for the worked example: REBOUND's own `ias15` and `ias15_cft` on the
   same problem, the same fixed step, printed as exact hex floats. At
   binary64 they must be identical, which is parcel A's gate seen from
   Python.

Step 2 is done on all three platforms. `make python-lib` builds
`build/libcft_ias15.so` (`.dylib` on macOS, `.dll` on Windows) and
`make check-python` runs the equivalence example against it. Measured
on Linux, 2026-09-10, against the rebound 5.1.1 wheel:

```
IDENTICAL: 13 values, bit for bit.
```

and on Windows 11 with MSYS2 mingw64 gcc 16.1.0 and CPython 3.12.9,
2026-09-11, against the `cp312-win_amd64` wheel of the same 5.1.1:

```
IDENTICAL: 13 values, bit for bit.
```

The 13 are the same 13. `build/ias15_ref --problem
data/problems/kepler.txt --dt 0.05 --epsilon 0 --steps 252 --sample
252` on the same host prints the same hex floats, and that program is
REBOUND's own ias15 compiled from the
pinned clone by mingw and linked statically - so the wheel's
MSVC-compiled ias15, the mingw-compiled ias15, and `ias15_cft` agree to
the last bit, and the Windows numbers are the project's numbers rather
than a locally consistent set.

**And the control, which is what stops that being a pass against
itself.** `make check-python-control` is the same run at binary128 with
`--expect-differ`, and it must NOT be identical - if a format change
leaves 13 values unmoved then the format is being accepted and ignored,
and the binary64 comparison above is binary64 against binary64. On
Windows, 2026-09-11:

```
DIFFERS: 9 of 13 values, largest absolute difference 4.263e-14

Expected, and this is the control passing: fp128 changed the answer, so
the format reaches the arithmetic and the binary64 comparison is a real
one.
```

(Nine and not thirteen because `kepler.txt` is planar: the four `z` and
`vz` values are exactly zero at every format and cannot move. The other
nine - `t` and the eight in-plane coordinates - all do, `t` included,
because the clock is accumulated in the wide format too.)

Two different C runtimes on the two sides of that call and the same
bits out of it, which is less lucky than it sounds: the only libm
functions IAS15 reaches are `sqrt`, `fabs`, `copysign` and `isnormal`,
and the first is correctly rounded by IEEE-754 while the rest are
exact. `pow(x, 1./7.)` would have been the exception, and REBOUND does
not call it - `integrator_ias15.c:195` replaces it with a Newton
iteration for that reason, in its own words "machine independent". The
run is at a fixed step because that is what makes the comparison a
comparison; the adaptive step controller uses only those same four, but
it has not been run from Python on Windows and is not claimed here.

That is REBOUND's own `ias15` against `ias15_cft` selected by name from
Python, same problem, same fixed step, at binary64. The wide formats
run too, and the format reaches the arithmetic rather than being
accepted and ignored - 252 steps take 5.2 s at binary64, 15.3 s at
binary128 and 49.4 s at binary256, which is the software backend's
known ratio.

**Four things it needed, and one of them was a correction.**

1. **PIC objects of the four drop-in sources**, which is `PIC_OBJ` in
   the Makefile - the same list as `DROPIN_OBJ` with `-fPIC`, the
   archive added, and `src/cft_ias15_shared.c`. (`-fPIC` is `$(SHPIC)`
   now, and empty on Windows, where all code is position independent
   and gcc only warns that the flag means nothing.)

2. **libcft as a shared object, not the archive.** *This is a
   correction.* An earlier version of this section said `libcft.a`
   links into a shared object as it is, on the evidence that `device.o`
   carries no absolute relocations. One member is not the archive:
   `backend_xrt.o` is C++ and the link stops on `R_X86_64_PC32 against
   symbol _ZSt7nothrow@@GLIBCXX_3.4`. cft-fp256 builds a proper
   `libcft.so` from its own PIC objects, with the card backend when
   `XRT=1`, so the Makefile builds and links that. **On Windows the
   correction does not apply and the original claim is true** - see
   "The Windows rule".

3. **REBOUND's symbols left undefined**, resolved because
   `cft_rebound.load()` promotes the wheel's librebound to
   `RTLD_GLOBAL` first. The library registers itself in a constructor
   (`src/cft_ias15_shared.c`); there is nowhere else it could happen,
   because a Python caller's first contact is assigning a string to
   `sim.integrator`.

4. **A way to configure it, which did not exist.** Selecting the
   integrator always worked. The next line did not, though not for the
   reason first recorded here: `sim.integrator` resolves a field by the
   name its descriptor gives, and ours are all `cft_`-prefixed, so
   `sim.integrator.epsilon` - the spelling every REBOUND example uses -
   does not resolve, and nothing in REBOUND's Python layer knows to try
   `cft_epsilon`. A reader who knew to try it could always have
   configured the integrator; a reader following any REBOUND example
   could not. Python could reach
   binary64 with the default step control and nothing else - the one
   thing this repository is for was unreachable from Python.
   `cft_ias15_configure()` takes format, epsilon and max_iter by value,
   which is what a ctypes caller can call, and
   `cft_rebound.configure(lib, sim, format="fp256", epsilon=0.0)` is the
   Python side.

**One guard.** The library is compiled against the pinned REBOUND
headers and loaded into a process running the wheel's REBOUND, so the
two must be the same source or `struct reb_simulation` is laid out
differently on each side of the call and nothing says so. They are
identical today - the rebound 5.1.1 wheel ships a `rebound.h` byte for
byte the same as the pinned bdfda4bd, on Windows as well as Linux - and
that is luck rather than design, so `check-rebound-match` is an
**order-only prerequisite of the shared library itself**. It therefore
runs before the link rather than beside it, and `check-python` inherits
it; order-only is what lets a phony prerequisite gate a target without
forcing a relink every time. Verified by appending one comment line to
the wheel's `rebound.h`: the target refuses, no library is produced,
and `check-python` refuses too. The line was then removed.

## The Windows rule

A Windows DLL may not carry undefined symbols. Everything else follows
from that one sentence, and it inverts the POSIX design rather than
adjusting it.

**The build links REBOUND; the loader does not.** `make python-lib`
asks `$(PYTHON)` for `rebound.__libpath__` and puts that file on the
link line. Nothing else changes on the C side - `rebound.h` marks its
API `__declspec(dllimport)` unless `BUILDINGLIBREBOUND` is defined, and
`PIC_OBJ` has never defined it, so the objects were already asking to
import these symbols rather than export them. The `.pyd`'s export
table has 494 entries, `reb_integrator_register` among them, so the
link resolves.

**The registration still happens in the constructor.**
`__attribute__((constructor))` in `src/cft_ias15_shared.c` is run by
mingw from `DLL_PROCESS_ATTACH`, so `ctypes.CDLL(path)` registers
`ias15_cft` exactly as `dlopen` does. No source change was needed for
Windows; the whole port is the Makefile and `python/cft_rebound.py`.

**No `RTLD_GLOBAL` dance.** The Windows loader resolves an import
descriptor by base name against the modules already in the process, and
`cft_rebound.load()` imports `rebound` first, so the import binds to
the `.pyd` the caller's simulations live in. That is not an
optimisation - site-packages is on no DLL search path, so if the module
were *not* already loaded the loader would go to disk and fail.

**libcft goes in as the archive.** The correction above - that
`libcft.a` will not enter a shared object because `backend_xrt.o` is
C++ - is a fact about ELF, and about a member that on Windows is not
there: `backend_xrt.o` is built only under `XRT=1`, XRT is Linux-only,
and a Windows `libcft.a` is fifteen C objects. So the DLL swallows it
whole and there is no `cft.dll` to find beside it. **This is the
Windows equivalent of "build libcft.so and link that": don't - link the
archive, and ship one file.**

**What it costs: the import table names one exact file.** After the
link, `objdump -p build/libcft_ias15.dll` reports

```
DLL Name: KERNEL32.dll
DLL Name: librebound.cp312-win_amd64.pyd
DLL Name: msvcrt.dll
```

and that middle name carries the Python ABI tag, so **a DLL built for
CPython 3.12 cannot be used by 3.13**. `make python-lib PYTHON=...`
with the interpreter that will load it is the whole remedy, but the
failure it prevents is nasty enough to be worth naming, because
Windows reports it in two different unhelpful ways and neither of them
mentions Python:

- if no matching librebound is anywhere on the search path,
  `ctypes.CDLL` raises `FileNotFoundError: Could not find module
  '...libcft_ias15.dll' (or one of its dependencies)` - which names our
  library, not the dependency that is missing;
- if one *is* findable - and `ctypes` adds
  `LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR`, so a librebound sitting beside
  our DLL qualifies - the load
  **succeeds silently**, registering the integrator into a second
  REBOUND with its own integrator list, and the only symptom is
  `RuntimeError: Integrator not found.` on the next line. That is the
  "two REBOUNDs" failure the POSIX section warns about, reached by a
  different road.

Both were measured. So `cft_rebound.load()` reads the DLL's import
table itself - `_pe_imported_dlls()`, about fifty lines of `struct`
against the PE headers - and refuses before `CDLL` is called:

```
libcft_ias15.dll was linked against librebound.cp313-win_amd64.pyd and
this interpreter has librebound.cp312-win_amd64.pyd.
...
    make python-lib PYTHON=C:\...\python.exe
```

`cft_rebound.linked_librebound(path)` is that answer on its own, and is
`None` on POSIX by construction - the shared object names no REBOUND at
all, which is the freedom the POSIX rule buys with its loader dance.

**One thing the Windows build does not fix.** The `.pyd` is MSVC-built
and uses the UCRT; MSYS2's `mingw64` gcc links `msvcrt.dll`. Two C
runtimes means two heaps, and `adopt_loaded_state()` in
`src/reb_integrator_cft.c` calls `free()` on blobs REBOUND's
Simulationarchive loader allocated. Nothing on the equivalence path
does that, and nothing here has crashed, but **reading a cft
Simulationarchive from Python on Windows is not gated and should be
treated as untested.** The fix costs no source change: build with
MSYS2's **UCRT64** toolchain instead, and both sides are `ucrtbase`.
That toolchain is not installed on the machine this was measured on, so
the claim is reasoned, not run.

**And one trap that is not ours.** MSYS2 `make` hands recipes a
stripped environment - the same reason `TMP` and `TEMP` are make
variables in this project's Makefile. If your `rebound` is reachable
only through `PYTHONPATH` or a conda activation, pass it the same way:
`make PYTHONPATH=... python-lib`. A `pip install rebound` into the
interpreter you name as `PYTHON=` needs none of this.

## Why there is no pyproject.toml

It was considered when Windows landed and rejected, for one reason in
three registers: **this library's ABI depends on a package, not on a
platform, and no Python packaging tag can say so.**

- **A binary wheel cannot be published.** The tags a wheel carries are
  interpreter, ABI and platform. What this library is actually
  compatible with is *the `rebound.h` in the REBOUND wheel you have* -
  a struct layout, which `check-rebound-match` exists to police - and
  on Windows also *the exact `librebound.<abi>.pyd`*. A wheel that
  claimed `cp312-win_amd64` and nothing else would install cleanly onto
  the wrong REBOUND and be wrong in silence, which is precisely the
  failure this project spends a Makefile target preventing.
- **A source distribution would be a build system in disguise.**
  Building at install time means cloning REBOUND and cft-fp256 at the
  commits in `third_party/MANIFEST` and compiling libcft, from inside
  pip's isolated build environment, with whatever C toolchain happens
  to be there. The Makefile does that in the open, with the pin
  asserted and hard-failing (`tools/fetch_third_party.sh`). Hiding it
  behind `pip install .` would make it no more reproducible and much
  harder to see fail.
- **A pure-Python wheel of `cft_rebound.py` alone would install a
  loader for a file it cannot produce.** The one thing the module
  needs to be told is where the library is; shipping it through pip
  removes the tree that is the only place it could look.

So the unit of distribution stays the repository, the unit of build
stays `make python-lib`, and the Python side stays one file on
`sys.path`:

    import sys; sys.path.insert(0, "path/to/cft-rebound/python")
    import cft_rebound
    cft_rebound.load()          # finds build/libcft_ias15.{so,dylib,dll}

`load()` taking no argument is the concession that replaces the
packaging: `cft_rebound.default_library()` looks beside the module it
was imported from, which is the only path anybody can guess correctly.

## Four sharp edges, all measured

**The name must be lower case.** `__setattr__` calls `value.lower()`
before the string reaches C, and the C lookup is `strcmp`. A name with
an upper-case letter registers fine and can never be selected from
Python: `sim.integrator = "Stub_CFT"` raises
`RuntimeError: Integrator not found.` `ias15_cft` is safe.

**Only one custom integrator per process.** The scan in
`reb_integrator_register` (REBOUND's own src/rebound.c) increments its
index before testing, so the second registration in a process calls
`strcmp(NULL, name)` on the list's `{0}` terminator. At `-O0` that
segfaults; at `-O2` and above - which is how every REBOUND wheel is
built - gcc infers from `strcmp`'s `nonnull` attribute that the
terminator test can never fail, and **the loop never ends**. A second
`reb_integrator_register` in one process hangs, in C and in Python
alike. Registering a name that collides with a *built-in* is fine: it
is refused with `Error! Integrator name must be unique but name already
exists.` and nothing is left broken. This is an upstream bug, not ours,
but it means: load the library once, and never ship two custom
integrators in one library.

**A stock REBOUND reads a cft archive, but the message is wrong and an
error is left behind.** Reading a Simulationarchive whose
`integrator.name` is not registered in that process:
`reb_simulation_set_integrator` fails, so `r->integrator.name` stays
`ias15`, and the next field - `integrator.ias15_cft.<something>` - trips
the prefix check in `binarydata.c`, which sets
`REB_BINARYDATA_WARNING_CORRUPTFILE` and **stops reading the snapshot**.
Because the integrator's fields are written last, the binary64 state -
`t`, `dt`, every particle - is recovered correctly, which is the
behaviour the roadmap wants. But the user is told
`RuntimeWarning: The binary file seems to be corrupted`, which is not
true, and the `Integrator not found.` error stays queued in
`r->messages`: nothing on the load path drains it, so it is raised by
the *next* call that runs `process_messages()`, on an unrelated line.
`REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR` exists in binarydata.h and
`reb_binarydata_process_warnings` prints for it, but nothing in REBOUND
ever sets that bit - the friendly message the roadmap quotes is dead
code at bdfda4bd.

So: **call `cft_rebound.load()` before you open the archive.** Doing it
afterwards, the way binarydata.c's unreachable message advises, does
recover - `sim.process_messages()` to swallow the stale error, then
`sim.integrator = "ias15_cft"` - but `set_integrator` calls your
`create()`, so the integrator's own archived state is not there:
measured, a field the archive held at 8 read back as 0. `t`, `dt` and
the particles are intact.

**Loading twice is safe; loading two copies is not.** Calling
`cft_rebound.load()` on the same path twice is a no-op - `dlopen` and
`LoadLibrary` both return the existing handle without re-running the
constructor. Loading two different files that each register the same
name is the second-registration hang above. Windows has a second way
into it, described in "The Windows rule": a DLL built for the wrong
interpreter, loaded next to a librebound the loader can find, silently
registers into a *second* REBOUND rather than hanging. `load()` refuses
that case before `CDLL` sees it.

## What has not been tested

Two of the four entries here were closed later the same day; they are
kept, struck through, because what a document once could not show is
part of its record. docs/VALIDATION.md entry 28 is the measurement.

- ~~**Parcel A's integrator.**~~ **Closed.** The sections above were
  first measured with a probe integrator written for the purpose (a
  drift-only step with two archivable state fields), because the
  load-bearing question was reachability, not arithmetic, and
  `ias15_cft` did not exist yet. It does: `make check-python` runs
  `ias15_cft` against REBOUND's own `ias15` from Python on the same
  problem and the same fixed step and gets 13 values identical bit for
  bit, and the wide formats reach the arithmetic (the timings above).
- **macOS.** The mechanism is the POSIX one, and the wheel sets
  `-Wl,-install_name,@rpath/librebound<suffix>`, but no run was made.
  The Darwin branch of the Makefile now passes
  `-Wl,-rpath,@loader_path` where it used to pass `$$ORIGIN`, which is
  the Linux spelling and means nothing to `dyld`; that is a correction
  made by reading, not by running.
- **Windows at `epsilon > 0`, and the Windows archive path.** `make
  check-python` and `make check-python-control` are binary64 and
  binary128 at a fixed step; binary256 was run once by hand
  (`--format fp256 --expect-differ --steps 40`: 2 of 13 values move,
  largest difference 8.882e-16, so the widest format reaches the
  arithmetic there too) but is not gated, because 252 steps of software
  binary256 is a minute. The adaptive step has not been run from Python
  on Windows at all. Nor has the archive path - see the CRT note in
  "The Windows rule", which is the one place where the Windows build is
  reasoned rather than measured.
- **Windows on a Python other than 3.12.** The build is asserted to be
  tied to one minor version, and `load()` refuses a mismatch by name -
  demonstrated by relinking against a copy of the `.pyd` whose internal
  export name was byte-patched to `cp313` (renaming the file is not
  enough; `ld` records the name from the PE export directory, not the
  filesystem). No actual 3.13 wheel was installed.
- ~~**The card.**~~ **Partly closed, and not from Python.** The shared
  library links `libcft.so`, which carries the card backend when
  cft-fp256 is built with `XRT=1`, and the shim it contains resolves
  `$CFT_REBOUND_ARTIFACT` like every other entry point - so a Python
  caller *can* reach a tile. Nobody has: the hardware run of entry 27
  was `check_dropin` and `gate_real`, both C. No Python run has opened
  an artifact or a `cft://` server.
- **Simulationarchive round trips with the wide `cft_` state, from
  Python.** Parcel B's fields exist now and `tests/gate_real.c` and
  `tools/check_checkpoint.py` gate the round trip in C, at all three
  formats and across processes. What was shown *here* is only that a
  custom integrator's `field_descriptor_list` survives a write and a
  read through the Python API, including five snapshots of an archive
  read back with `rebound.Simulationarchive` - with the probe's two
  fields, not with the 50 wide blobs.
