# SPDX-License-Identifier: GPL-3.0-or-later
"""Reach a cft-rebound integrator from REBOUND's Python package.

REBOUND's Python layer does not keep its own list of integrators. It
lowercases the string you assign and hands it straight to C
(``rebound/simulation.py``, ``Simulation.__setattr__``)::

    clibrebound.reb_simulation_set_integrator(byref(self), c_char_p(value.encode("ascii")))

and the C side looks the name up in the registry that
``reb_integrator_register`` writes. So an integrator registered from a
shared library is reachable by name from Python, with no binding layer
and no change to REBOUND, provided the registration lands in the *same*
librebound the Python package already loaded. Making that so is the
whole job of this module, and it differs by platform:

* **Windows: not built.** The loading rule would be simpler - REBOUND
  ships librebound as a ``.pyd`` whose API is ``__declspec(dllexport)``,
  so a DLL that imports ``reb_integrator_register`` binds to the module
  already in the process and ``ctypes.CDLL(path)`` is enough - but a DLL
  may not carry undefined symbols, so it would have to LINK against
  ``rebound.__libpath__``, which is a different build rule. The
  Makefile's shared-library target is POSIX only and nothing here
  produces a Windows DLL. See docs/PYTHON.md, "Windows is not done".

* **Linux and macOS.** Not automatic. ``rebound/__init__.py`` loads
  librebound with ``cdll.LoadLibrary``, which is ``RTLD_LOCAL``, so its
  symbols are not in the global scope and cannot resolve ours::

      OSError: libcft_ias15.so: undefined symbol: reb_integrator_register

  Re-opening the same file with ``RTLD_GLOBAL`` promotes it into the
  global scope - it is the same mapping, not a second copy - and our
  library then loads and registers into the very same registry.

Usage::

    import rebound
    import cft_rebound

    cft_rebound.load("build/libcft_ias15.so")   # `make python-lib`
    sim = rebound.Simulation()
    sim.integrator = "ias15_cft"

Run this file to see what your REBOUND accepts::

    python python/cft_rebound.py [path-to-library]
"""

import ctypes
import os

__all__ = ["load", "registered", "include_dir", "library_path"]


def load(path):
    """Load a shared library that registers a REBOUND integrator.

    ``path`` is the library built from this repository: ``make
    python-lib`` writes ``build/libcft_ias15.so`` on Linux and
    ``build/libcft_ias15.dylib`` on macOS. There is no Windows build.
    Registration happens in the
    library's own constructor, so the name is usable as soon as this
    returns. The handle is returned so the caller can keep it alive and
    read any symbols of its own; do not let it be garbage collected
    while a simulation is using the integrator, and never unload it -
    REBOUND stores the registered name and callbacks by pointer.

    Raises ``OSError`` if the library cannot be loaded.
    """
    import rebound

    if os.name != "nt":
        # Promote the already-loaded librebound to the global scope so
        # that our library's undefined REBOUND symbols resolve against
        # it. dlopen on an already-loaded object returns that same
        # object; RTLD_GLOBAL adds it to the global search scope.
        ctypes.CDLL(rebound.__libpath__, mode=ctypes.RTLD_GLOBAL)
    return ctypes.CDLL(os.fspath(path))


def configure(lib, sim, format="fp64", epsilon=0.0, max_iter=0):
    """Set the wide format and step control on a simulation.

    ``sim`` must already be using the integrator (``sim.integrator =
    "ias15_cft"``).

    Why this exists rather than an attribute assignment: the settings
    this integrator needs are not the ones REBOUND generates. REBOUND
    5.1.1 resolves ``sim.integrator.<field>`` against the registered
    ``field_descriptor_list`` (``rebound/integrator.py``
    ``__getattr__``/``__setattr__``), and every name in ours is
    ``cft_``-prefixed, so ``sim.integrator.epsilon`` - the spelling a
    REBOUND user knows - does not resolve. This call takes format,
    epsilon and max_iter by value, validates them together, and applies
    the per-format ``max_iter`` default that a raw field assignment
    would not.

    (An earlier version of this docstring said REBOUND exposes a
    built-in integrator as ``sim.ri_<name>`` and knows nothing of a
    custom one. Neither half is true of the pinned 5.1.1: ``ri_<name>``
    is the 3.x/4.x API and appears nowhere in its Python package except
    one unreachable line of ``citations.py``, and the descriptor-list
    route above makes no distinction between built-in and custom.
    Whether ``sim.integrator.cft_format = 2`` works has not been run.)

    ``lib`` is the handle ``load()`` returned. ``format`` is "fp64",
    "fp128" or "fp256". ``epsilon`` follows REBOUND: 0.0 is a fixed
    step. ``max_iter`` 0 means the default for the format - REBOUND's
    12 at binary64, more above it, because the corrector needs about 6
    passes at binary128 and 18-22 at binary256.

    Raises ValueError if the simulation is not using the integrator or
    a setting is refused.
    """
    code = lib.cft_ias15_format_code(format.encode("ascii"))
    if code < 0:
        raise ValueError("unknown format %r; use fp64, fp128 or fp256" % format)

    lib.cft_ias15_configure.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_int]
    lib.cft_ias15_configure.restype = ctypes.c_int
    rc = lib.cft_ias15_configure(
        ctypes.byref(sim), code, ctypes.c_double(epsilon), int(max_iter))
    if rc:
        raise ValueError(
            {1: "the simulation is not using this integrator; set "
                "sim.integrator = 'ias15_cft' first",
             2: "format %r is not one this integrator has" % format,
             3: "epsilon must not be negative",
             4: "max_iter must not be negative"}.get(rc, "refused (%d)" % rc))


def registered():
    """The integrator names this process's REBOUND will accept.

    Built-ins first, then anything registered by a loaded library.
    Names are matched exactly and ``sim.integrator = ...`` lowercases
    what you assign, so a registered name containing an upper-case
    letter can never be selected from Python.
    """
    import rebound

    f = rebound.clibrebound.reb_integrators_registered
    f.restype = ctypes.POINTER(ctypes.c_char_p)
    p = f()
    names = []
    i = 0
    while p[i]:
        names.append(p[i].decode("ascii"))
        i += 1
    rebound.clibrebound.reb_free(p)  # the list, not the names
    return tuple(names)


def include_dir():
    """Where the installed REBOUND keeps its headers.

    A wheel installs ``librebound.<abi>.so`` and a ``src/`` directory of
    headers side by side, so this is the ``-I`` a library that registers
    an integrator must be compiled against - and compiling against the
    headers of the REBOUND that is actually installed is the only way to
    be sure ``struct reb_simulation`` has the layout the loaded library
    expects.
    """
    import rebound

    return os.path.join(os.path.dirname(rebound.__libpath__), "src")


def library_path():
    """The librebound the Python package loaded. Link against this one."""
    import rebound

    return rebound.__libpath__


if __name__ == "__main__":
    import sys

    import rebound

    print("rebound       ", rebound.__version__, rebound.__githash__)
    print("librebound    ", library_path())
    print("headers       ", include_dir(),
          "(present)" if os.path.isdir(include_dir()) else "(MISSING)")
    if len(sys.argv) > 1:
        print("loading       ", sys.argv[1])
        load(sys.argv[1])
    names = registered()
    print("integrators   ", " ".join(names))
    custom = [n for n in names if n not in (
        "ias15", "whfast", "sei", "leapfrog", "janus", "mercurius",
        "saba", "eos", "bs", "whfast512", "trace", "none")]
    print("custom        ", " ".join(custom) if custom else "(none)")
