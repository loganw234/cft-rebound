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
githash 33549d1d. The numbers and the failures are in
docs/VALIDATION.md.

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
different one does, which is why there is no Windows build. REBOUND
ships librebound as a `.pyd` whose API is `__declspec(dllexport)`, so a
DLL that imports `reb_integrator_register` binds to the module already
loaded in the process, and a plain `ctypes.CDLL(path)` is enough: no
`RTLD_GLOBAL` dance. But a DLL may not carry undefined symbols, so it
has to LINK against `rebound.__libpath__` rather than leave REBOUND's
symbols to the loader, and its import table then names one exact file
(`librebound.cp312-win_amd64.pyd`), tying the build to the Python minor
version. The POSIX library, which links against nothing, is not. See
"Windows is not done" below.

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
   prints the REBOUND it found, the header directory, and the
   integrator names now registered.
4. `python python/example_equivalence.py --library path/to/the/library`
   for the worked example: REBOUND's own `ias15` and `ias15_cft` on the
   same problem, the same fixed step, printed as exact hex floats. At
   binary64 they must be identical, which is parcel A's gate seen from
   Python.

Step 2 is done. `make python-lib` builds
`build/libcft_ias15.so` (`.dylib` on macOS) and `make check-python`
runs the equivalence example against it. Measured on Linux,
2026-09-10, against the rebound 5.1.1 wheel:

```
IDENTICAL: 13 values, bit for bit.
```

That is REBOUND's own `ias15` against `ias15_cft` selected by name from
Python, same problem, same fixed step, at binary64. The wide formats
run too, and the format reaches the arithmetic rather than being
accepted and ignored - 252 steps take 5.2 s at binary64, 15.3 s at
binary128 and 49.4 s at binary256, which is the software backend's
known ratio.

**Four things it needed, and one of them was a correction.**

1. **PIC objects of the four drop-in sources**, which is `PIC_OBJ` in
   the Makefile - the same list as `DROPIN_OBJ` with `-fPIC`, the
   archive added, and `src/cft_ias15_shared.c`.

2. **libcft as a shared object, not the archive.** *This is a
   correction.* An earlier version of this section said `libcft.a`
   links into a shared object as it is, on the evidence that `device.o`
   carries no absolute relocations. One member is not the archive:
   `backend_xrt.o` is C++ and the link stops on `R_X86_64_PC32 against
   symbol _ZSt7nothrow@@GLIBCXX_3.4`. cft-fp256 builds a proper
   `libcft.so` from its own PIC objects, with the card backend when
   `XRT=1`, so the Makefile builds and links that.

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
byte the same as the pinned bdfda4bd - and that is luck rather than
design, so `make python-lib` runs `check-rebound-match` first and
refuses if they ever diverge.

**Windows is not done.** A DLL may not carry undefined symbols, so it
must link against `rebound.__libpath__` rather than leave REBOUND's
symbols to the loader. That is a different rule, and it deserves its own
change and its own test rather than an `ifeq` bolted onto a target that
works.

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
name is the second-registration hang above.

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
