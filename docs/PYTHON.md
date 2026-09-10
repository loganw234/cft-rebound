# Python

**Yes, and it costs a loader.** A REBOUND user in Python reaches a
custom C integrator by name, with no binding layer, no fork of
REBOUND's Python package and no change to it. This file is what was
measured, where the one obstacle is, and what a user installs.

    import rebound
    import cft_rebound                     # python/cft_rebound.py

    cft_rebound.load("build/libias15_cft.so")
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

`reb_simulation_set_integrator` (src/simulation.c) tries the built-ins
and then walks `reb_integrator_configurations_custom`, the global that
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
  `field_descriptor_list` - `sim.integrator.epsilon = 0.0` is how a
  Python user turns off adaptive stepping;
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
    >>> ctypes.CDLL("libias15_cft.so")
    OSError: libias15_cft.so: undefined symbol: reb_integrator_register

The fix is one line. Re-opening the same file with `RTLD_GLOBAL`
promotes the mapping already in the process - it is not a second copy;
the address of `reb_integrator_configurations_custom` seen through both
handles is the same - and our library then resolves against it:

    ctypes.CDLL(rebound.__libpath__, mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL("libias15_cft.so")        # registers "ias15_cft"

That is `cft_rebound.load()`, and it is the entire Python-side cost.

**On Windows there is no obstacle at all.** REBOUND ships librebound as
a `.pyd` whose API is `__declspec(dllexport)`, so a DLL that imports
`reb_integrator_register` binds to the module already loaded in the
process, and a plain `ctypes.CDLL(path)` is enough. The price is that
the DLL's import table names one exact file
(`librebound.cp312-win_amd64.pyd`), so a Windows build is tied to the
Python minor version it was linked against. The POSIX library, which
links against nothing, is not.

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

Step 2 is the part this repository still owes. `libcft` must be built
`-fPIC` and linked in, and today `src/ias15_cft.c` is a program with a
`main`, not a library - so the Makefile has no shared-library target
yet. That is packaging (parcel C), not Python.

## Four sharp edges, all measured

**The name must be lower case.** `__setattr__` calls `value.lower()`
before the string reaches C, and the C lookup is `strcmp`. A name with
an upper-case letter registers fine and can never be selected from
Python: `sim.integrator = "Stub_CFT"` raises
`RuntimeError: Integrator not found.` `ias15_cft` is safe.

**Only one custom integrator per process.** The scan in
`reb_integrator_register` (src/rebound.c) increments its index before
testing, so the second registration in a process calls
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

- **Parcel A's integrator.** It did not exist when this was written.
  Everything above was measured with a probe integrator written for the
  purpose (a drift-only step with two archivable state fields), because
  the load-bearing question is reachability, not arithmetic. The
  equivalence example has been exercised on its identical path
  (`ias15` against `ias15`, 13 values, bit for bit) and on its differing
  path (against the probe), never against `ias15_cft`.
- **macOS.** The mechanism is the POSIX one, and the wheel sets
  `-Wl,-install_name,@rpath/librebound<suffix>`, but no run was made.
- **The card.** Nothing here goes near XRT, an artifact or `cft://`.
- **Simulationarchive round trips with the wide `cft_` state.** Parcel
  B's fields did not exist either. What was shown is that a custom
  integrator's `field_descriptor_list` survives a write and a read
  through the Python API, including five snapshots of an archive read
  back with `rebound.Simulationarchive`.
