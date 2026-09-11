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
whole job of this module, and it differs by platform - the two rules
are opposites, and each one's cost is the other's freedom:

* **Linux and macOS: resolve late.** ``rebound/__init__.py`` loads
  librebound with ``cdll.LoadLibrary``, which is ``RTLD_LOCAL``, so its
  symbols are not in the global scope and cannot resolve ours::

      OSError: libcft_ias15.so: undefined symbol: reb_integrator_register

  Re-opening the same file with ``RTLD_GLOBAL`` promotes it into the
  global scope - it is the same mapping, not a second copy - and our
  library then loads and registers into the very same registry. The
  library itself names no REBOUND at build time, so one build serves
  every Python that ships a compatible ``rebound.h``.

* **Windows: resolve early.** A DLL may not carry undefined symbols,
  so ``make python-lib`` LINKS ``rebound.__libpath__`` and the import
  table then names that one file - ``librebound.cp312-win_amd64.pyd``.
  No ``RTLD_GLOBAL`` dance is needed, because the Windows loader binds
  an import descriptor to a module already in the process by its base
  name and ``import rebound`` has already loaded it. The price is that
  the base name carries the ABI tag: **a DLL built for one Python
  minor version cannot be loaded by another.** The loader's own account
  of that is ``OSError: [WinError 126] The specified module could not
  be found``, which names neither the module nor the fix, so ``load()``
  reads the DLL's import table first and says both.

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
import struct

__all__ = ["load", "registered", "include_dir", "library_path",
           "linked_librebound", "default_library"]


def _pe_imported_dlls(path):
    """The DLL names in a PE image's import table, or None.

    Windows resolves an import descriptor by the base name it carries,
    so this is exactly the list of files the loader will go looking
    for - and for this library one of them is the librebound it was
    linked against, which is the one fact a failed load never tells
    you. Anything unexpected returns None rather than raising: a
    diagnostic that guesses is worse than the OSError it replaces.
    """
    try:
        with open(path, "rb") as f:
            # Two bytes before committing to the whole file: this is
            # called on POSIX too, where the answer is always None and
            # the file is a megabyte of ELF.
            if f.read(2) != b"MZ":
                return None
            f.seek(0)
            data = f.read()
        pe = struct.unpack_from("<I", data, 0x3c)[0]
        if data[pe:pe + 4] != b"PE\0\0":
            return None
        n_sections, = struct.unpack_from("<H", data, pe + 6)
        opt_size, = struct.unpack_from("<H", data, pe + 20)
        opt = pe + 24
        magic, = struct.unpack_from("<H", data, opt)
        if magic == 0x20b:            # PE32+
            data_dirs = opt + 112
        elif magic == 0x10b:          # PE32
            data_dirs = opt + 96
        else:
            return None
        # Data directory entry 1 is the import table.
        import_rva, = struct.unpack_from("<I", data, data_dirs + 8)
        if not import_rva:
            return []

        sections = []
        for i in range(n_sections):
            base = opt + opt_size + 40 * i
            vsize, vaddr, rawsize, rawptr = struct.unpack_from("<IIII", data, base + 8)
            sections.append((vaddr, max(vsize, rawsize), rawptr))

        def file_offset(rva):
            for vaddr, size, rawptr in sections:
                if vaddr <= rva < vaddr + size:
                    return rawptr + (rva - vaddr)
            return None

        names = []
        desc = file_offset(import_rva)
        if desc is None:
            return None
        while True:
            entry = data[desc:desc + 20]
            if len(entry) < 20 or entry == b"\0" * 20:
                break
            at = file_offset(struct.unpack_from("<I", entry, 12)[0])
            if at is not None:
                names.append(data[at:data.index(b"\0", at)].decode("ascii", "replace"))
            desc += 20
        return names
    except Exception:
        return None


def linked_librebound(path):
    """The librebound this library was built against, or None.

    On Windows that is a real answer and a load-bearing one: the DLL
    imports one exact file name and can be loaded by no other Python
    minor version. On POSIX it is None by construction - the shared
    object names no REBOUND at all and resolves against whichever one
    the process already has open, which is why one build serves every
    interpreter there.
    """
    for name in _pe_imported_dlls(path) or ():
        if name.lower().startswith("librebound"):
            return name
    return None


def _refuse_wrong_librebound(path, wanted, have):
    import sys
    raise OSError(
        "%s was linked against %s and this interpreter has %s.\n"
        "A Windows DLL may not carry undefined symbols, so this library links "
        "the librebound it will be loaded into, and the import table then names "
        "one exact file - which carries the Python ABI tag. Build it with the "
        "interpreter that will load it:\n"
        "    make python-lib PYTHON=%s\n"
        "(On Linux and macOS this cannot happen: the shared object names no "
        "REBOUND and resolves against the one already in the process.)"
        % (os.path.basename(path), wanted, have, sys.executable))


def default_library():
    """Where ``make python-lib`` left it, or None.

    Only ``build/`` in the tree this module was imported from. There is
    no system-wide location to guess, deliberately (docs/PYTHON.md,
    "Why there is no pyproject.toml"): the library is built against the
    REBOUND *you* have, so a copy found by search could be one built
    against a different one. ``make install-python-lib`` will put it in
    LIBDIR if you want it there - pass that path to ``load()``.
    """
    names = (("libcft_ias15.dll",) if os.name == "nt"
             else ("libcft_ias15.so", "libcft_ias15.dylib"))
    build = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "build")
    for name in names:
        candidate = os.path.join(build, name)
        if os.path.isfile(candidate):
            return candidate
    return None


def load(path=None):
    """Load a shared library that registers a REBOUND integrator.

    ``path`` is the library built from this repository: ``make
    python-lib`` writes ``build/libcft_ias15.so`` on Linux,
    ``build/libcft_ias15.dylib`` on macOS and
    ``build/libcft_ias15.dll`` on Windows. Omit it to use
    ``default_library()``, which is that file in the tree this module
    came from. Registration happens in the
    library's own constructor - ``__attribute__((constructor))``, which
    mingw runs from ``DLL_PROCESS_ATTACH`` - so the name is usable as
    soon as this returns. The handle is returned so the caller can keep
    it alive and read any symbols of its own; do not let it be garbage
    collected while a simulation is using the integrator, and never
    unload it - REBOUND stores the registered name and callbacks by
    pointer.

    ``import rebound`` happens first on both platforms, and for
    opposite reasons. On POSIX it is what there is to promote to
    ``RTLD_GLOBAL``. On Windows it is what puts
    ``librebound.<abi>.pyd`` in the loaded-module list, so that the
    loader binds our import descriptor to *that* module by base name
    and never goes to disk - where it would not find it, since
    site-packages is on no DLL search path.

    Raises ``OSError`` if the library cannot be loaded.
    """
    import rebound

    if path is None:
        path = default_library()
        if path is None:
            raise OSError(
                "no library given and none built. Run `make python-lib`, or "
                "pass the path to the shared library that registers the "
                "integrator.")
    path = os.fspath(path)
    if os.name == "nt":
        wanted = linked_librebound(path)
        have = os.path.basename(rebound.__libpath__)
        if wanted is not None and wanted.lower() != have.lower():
            _refuse_wrong_librebound(path, wanted, have)
        # No RTLD_GLOBAL dance: the .pyd's API is __declspec(dllexport)
        # and it is already loaded, so the import binds to it. (For a
        # name with a separator ctypes adds
        # LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR to
        # LOAD_LIBRARY_SEARCH_DEFAULT_DIRS, so any OTHER dependency is
        # looked for beside this file - and PATH is not searched at
        # all, which is why an unmatched build cannot be rescued by
        # putting site-packages on PATH.)
        return ctypes.CDLL(path)

    # Promote the already-loaded librebound to the global scope so
    # that our library's undefined REBOUND symbols resolve against
    # it. dlopen on an already-loaded object returns that same
    # object; RTLD_GLOBAL adds it to the global search scope.
    ctypes.CDLL(rebound.__libpath__, mode=ctypes.RTLD_GLOBAL)
    return ctypes.CDLL(path)


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
    """The librebound the Python package loaded.

    On Windows this is what the DLL must be linked against, and
    ``make python-lib`` asks the interpreter for it rather than being
    told. On POSIX nothing links it; it is the file ``load()`` promotes
    to ``RTLD_GLOBAL``.
    """
    import rebound

    return rebound.__libpath__


if __name__ == "__main__":
    import sys

    import rebound

    print("python        ", sys.version.split()[0], "on", os.name)
    print("rebound       ", rebound.__version__, rebound.__githash__)
    print("librebound    ", library_path())
    print("headers       ", include_dir(),
          "(present)" if os.path.isdir(include_dir()) else "(MISSING)")
    lib = sys.argv[1] if len(sys.argv) > 1 else default_library()
    if lib:
        print("loading       ", lib)
        # What it was built against, printed BEFORE the load, because
        # on Windows a mismatch here is the whole reason the load is
        # about to fail and the loader will not say so.
        linked = linked_librebound(lib)
        if linked is not None:
            print("linked against", linked,
                  "(ok)" if linked.lower()
                  == os.path.basename(library_path()).lower() else "(MISMATCH)")
        elif os.name == "nt":
            print("linked against (no librebound in the import table)")
        else:
            print("linked against (nothing - POSIX resolves at load time)")
        load(lib)
    else:
        print("loading        (nothing: give a path, or run `make python-lib`)")
    names = registered()
    print("integrators   ", " ".join(names))
    custom = [n for n in names if n not in (
        "ias15", "whfast", "sei", "leapfrog", "janus", "mercurius",
        "saba", "eos", "bs", "whfast512", "trace", "none")]
    print("custom        ", " ".join(custom) if custom else "(none)")
