# SPDX-License-Identifier: GPL-3.0-or-later
#
# cft-rebound: IAS15 with its arithmetic through libcft.
#
#   make third-party      clone/verify the pinned upstreams (needs git)
#   make libcft           build libcft from the pinned cft-fp256 clone
#   make librebound       build REBOUND's C library from the pinned clone
#   make                  the reference program, the libcft port, the
#                         drop-in gate and the two libraries
#   make constants        re-derive and re-check the IAS15 constants
#   make programs         generate and assemble the sequencer programs
#   make check            every gate at every format. About two and a
#                         half hours on this host, and it NEVER skips:
#                         it is the thing whose meaning must not change.
#   make check-quick      the same gates with the cheaper arguments, and
#                         it SKIPS any leg whose inputs are byte for byte
#                         what they were when that leg last passed
#                         (tools/gate_cache.py). Measured: 42 minutes
#                         cold, 13 seconds when nothing has changed.
#                         Every skip prints itself by name.
#   make check-light      every gate, every case, with most step counts
#                         scaled to a tenth - a fast signal from
#                         everything rather than a thorough one.
#                         Measured: 31 min end to end from clean, and
#                         on check_dropin alone 159 s against 585, the
#                         same 130 cases, when the suite was smaller.
#                         Every banner prints the count it actually ran.
#                         Never a substitute for `make check`.
#   make gate-cache-report
#                         what the cache currently holds, and whether
#                         each stamp came from a full run or a quick one
#   make HIP=1 check-gpu  on a GPU host only: the record-making legs on
#                         hip:0, held to the software backend process by
#                         process (tools/check_gpu.py). Refuses by name
#                         where there is no hip:0, and is in no tier
#   make archive          the Simulationarchive gate programs
#   make check-archive    just those gates (seconds)
#   make dropin           the drop-in library and its gate
#   make python-lib       the shared library for Python; check-python
#                         runs the equivalence gate through it and
#                         check-python-control is its negative control.
#                         Linux, macOS and Windows - but Windows obeys
#                         a different rule and LINKS the wheel's
#                         librebound, so the DLL is tied to one Python
#                         minor version. docs/PYTHON.md, "The Windows
#                         rule".
#   make install-python-lib  that library into LIBDIR. Separate from
#                         install because it needs a REBOUND wheel.
#   make example          build the two worked examples and run them
#   make install          headers and libraries where a REBOUND build
#                         can find them; PREFIX and DESTDIR as usual
#   make uninstall        remove exactly what install put there
#
# Set CFT_REBOUND_ARTIFACT to an .xclbin and every gate that opens the
# engine opens the tile instead; no gate takes a flag for it.
#
# On Windows (MSYS2 mingw64 gcc from Git Bash) the same traps as
# cft-fp256's host/Makefile apply, and they are passed through here:
#
#   PATH="/c/msys64/mingw64/bin:$PATH" make CC=gcc OS=Windows_NT \
#        TMP=C:/Users/you/AppData/Local/Temp TEMP=C:/Users/you/AppData/Local/Temp \
#        PYTHON=C:/path/to/python.exe
#
# TMP/TEMP must be MAKE variables (gcc's subprocesses do not see a
# shell export from Git Bash), OS=Windows_NT selects libcft's dllexport
# branch, and PYTHON must be an interpreter with mpmath.
#
# The same trap catches PYTHONPATH, and python-lib and check-python are
# the targets that care: MSYS2 make strips the environment it hands a
# recipe, so a rebound reachable only through PYTHONPATH or a conda
# activation is invisible to $(PYTHON) here. Pass it the same way -
# `make PYTHONPATH=... python-lib`. A wheel installed into the
# interpreter named by PYTHON= needs none of this.

CC      ?= gcc
PYTHON  ?= python3
CFLAGS  ?= -O2
CSTD    := -std=c99
WARN    := -Wall -Wextra

# make install: the usual conventions. PREFIX is where the package
# belongs on the running system and is compiled into the library (it is
# how libcft_rebound finds the ias15_cft program); DESTDIR is a staging
# root prepended at copy time and is NOT compiled in, so a distribution
# can build once with PREFIX=/usr and install into DESTDIR=/tmp/pkg.
PREFIX     ?= /usr/local
DESTDIR    ?=
BINDIR     ?= $(PREFIX)/bin
LIBDIR     ?= $(PREFIX)/lib
INCLUDEDIR ?= $(PREFIX)/include
INSTALL    ?= install

REB     := third_party/rebound/src
CFT     := third_party/cft-fp256/host
CFTLIB  := $(CFT)/libcft.a
B       := build

# REBOUND's generic (non-Windows-MSVC) flags, plus -O2. -std=c99 also
# pins -ffp-contract=off, so no FMA contraction changes the doubles.
REB_CFLAGS := $(CSTD) $(CFLAGS) -Wpointer-arith -D_GNU_SOURCE -Wno-unknown-pragmas \
              -DGITHASH=$(shell cd $(REB) 2>/dev/null && git rev-parse HEAD || echo unknown)
REB_SRC := $(filter-out $(REB)/integrator_whfast512.c, $(wildcard $(REB)/*.c))
REB_OBJ := $(patsubst $(REB)/%.c, $(B)/rebound/%.o, $(REB_SRC)) $(B)/rebound/whfast512_stub.o

ifeq ($(OS),Windows_NT)
  EXE   := .exe
  LIBS  := -lws2_32
  CFT_MAKEVARS := OS=Windows_NT TMP=$(TMP) TEMP=$(TEMP)
  # rebound.h marks its API __declspec(dllimport) on Windows unless
  # BUILDINGLIBREBOUND is defined. We link REBOUND statically, so both
  # the library objects and the programs are compiled as "building"
  # (dllexport in a static archive is harmless; dllimport would demand
  # a DLL's __imp_ thunks that do not exist).
  REB_CFLAGS += -DBUILDINGLIBREBOUND
  REB_USEFLAGS := -DBUILDINGLIBREBOUND
else
  EXE   :=
  LIBS  := -lm
  CFT_MAKEVARS :=
endif

# The shared library for Python. Three platforms, and Windows obeys a
# different rule from the other two - see docs/PYTHON.md, "The Windows
# rule". In short: a DLL may not carry undefined symbols, so it LINKS
# the wheel's librebound instead of leaving REBOUND's symbols to the
# loader, and libcft goes in as the archive rather than as a second
# shared object.
#
#   SHLIB      what to build
#   SHFLAGS    how to link it
#   SHPIC      -fPIC, where that means something
#   CFT_SHLIB  the libcft this depends on (a prerequisite)
#   CFT_SHLINK how to put that libcft on the link line
#   SHRPATH    how the result finds libcft at run time
ifeq ($(OS),Windows_NT)
  SHLIB      := libcft_ias15.dll
  # Every global is exported. Nothing in our sources is marked
  # __declspec(dllexport) - cft.h's CFT_API is plain unless
  # CFT_BUILD_SHARED is defined and we do not define it - so ld would
  # auto-export anyway; saying it means a header that grows an export
  # mark later cannot silently turn cft_ias15_configure into a symbol
  # ctypes can no longer find.
  SHFLAGS    := -shared -Wl,--export-all-symbols
  SHPIC      :=
  # libcft goes in as the ARCHIVE, and can here for the reason it
  # cannot on POSIX: the member that stops a .so is backend_xrt.o, C++,
  # and it is only built when XRT=1 - XRT is Linux-only, so a Windows
  # libcft.a is fifteen C objects. The DLL is therefore one file with
  # no libcft.dll to find beside it.
  CFT_SHLIB  := $(CFTLIB)
  CFT_SHLINK := $(CFTLIB)
  SHRPATH    :=
  # Nothing extra beside the DLL: the archive is inside it.
  CFT_SHLIB_EXTRA :=
  # Nor a runtime DLL. The trailer's wall clock is clock_gettime()
  # (24900ee, 2026-09-14), which mingw-w64 takes from winpthreads - by
  # default as libwinpthread-1.dll, in MSYS2's bin/ and on no search
  # path a Python process uses, so ctypes refused the library "or one
  # of its dependencies" from that commit until 2026-10-31. Linked
  # static it is code inside the one file.
  SHSTATIC   := -Wl,-Bstatic -lwinpthread -Wl,-Bdynamic
else ifeq ($(shell uname -s 2>/dev/null),Darwin)
  SHLIB      := libcft_ias15.dylib
  SHFLAGS    := -dynamiclib -Wl,-install_name,@rpath/libcft_ias15.dylib
  SHPIC      := -fPIC
  CFT_SHLIB  := $(CFT)/libcft.dylib
  CFT_SHLINK := -L$(CFT) -lcft
  SHRPATH    := -Wl,-rpath,@loader_path -Wl,-rpath,$(abspath $(CFT))
  CFT_SHLIB_EXTRA := $(CFT_SHLIB)
else
  SHLIB      := libcft_ias15.so
  SHFLAGS    := -shared
  SHPIC      := -fPIC
  CFT_SHLIB  := $(CFT)/libcft.so
  CFT_SHLINK := -L$(CFT) -lcft
  SHRPATH    := -Wl,-rpath,'$$ORIGIN' -Wl,-rpath,$(abspath $(CFT))
  CFT_SHLIB_EXTRA := $(CFT_SHLIB)
endif

# What to hand cft_rebound.load(). MSYS/Cygwin make gives recipes POSIX
# paths and a native Windows Python cannot open /c/..., so convert.
ifeq ($(OS),Windows_NT)
  SHLIB_PATH := $(shell cygpath -m '$(CURDIR)' 2>/dev/null || echo '$(CURDIR)')/$(B)/$(SHLIB)
else
  SHLIB_PATH := $(CURDIR)/$(B)/$(SHLIB)
endif

# An XRT build of libcft.a carries backend_xrt.o, which is C++ and needs
# the XRT runtime at link time; without these a plain `make` fails with
# undefined references to xrt::bo and operator new, which names neither
# the cause nor the fix. Sourcing /opt/xilinx/xrt/setup.sh sets
# XILINX_XRT, so that is the signal. XRT=0 turns it off for a
# software-only libcft.a in a shell that happens to have XRT sourced;
# LIBS= still overrides the lot.
ifneq ($(OS),Windows_NT)
  ifneq ($(XILINX_XRT),)
    XRT ?= 1
  endif
  XRT ?= 0
  ifeq ($(XRT),1)
    LIBS += -L$(XILINX_XRT)/lib -lxrt_coreutil -lstdc++ -lpthread -luuid
    # and libcft itself must be built with the backend, which it is not
    # by default. Every hardware run this project has done used an
    # archive built with XRT=1 by hand; this is that, from the signal.
    CFT_MAKEVARS += XRT=1 XRT_ROOT=$(XILINX_XRT)
  endif
endif

# A HIP build of libcft.a (cft-fp256's GPU backend, G1B: the `hip:N`
# device; its docs/HOSTAPI.md, "The GPU backend") carries backend_hip.o,
# which is C++ and needs the GPU runtime at link time: libamdhip64 on
# ROCm, the CUDA runtime under HIP's NVIDIA layer. Nothing in the
# environment says to build it, so it is HIP=1 by hand. HIP_ROOT,
# HIP_PLATFORM, HIP_CUDA and HIP_ARCH pass through to libcft's Makefile,
# whose defaults are /opt/rocm, amd, /usr/local/cuda and gfx1102 (89 120
# on nvidia). One libcft holds a card or a GPU, never both (round 5, P0).
ifneq ($(OS),Windows_NT)
  HIP ?= 0
  ifeq ($(HIP),1)
    ifeq ($(XRT),1)
      $(error HIP=1 with XRT=1: one libcft holds a card (XRT=1) or a GPU (HIP=1), never both - build them in two trees)
    endif
    HIP_ROOT     ?= /opt/rocm
    HIP_PLATFORM ?= amd
    HIP_CUDA     ?= /usr/local/cuda
    ifeq ($(HIP_PLATFORM),nvidia)
      LIBS += -L$(HIP_CUDA)/lib64 -Wl,-rpath,$(HIP_CUDA)/lib64 -lcudart -lstdc++ -lpthread -ldl -lrt
    else
      LIBS += -L$(HIP_ROOT)/lib -Wl,-rpath,$(HIP_ROOT)/lib -lamdhip64 -lstdc++
    endif
    CFT_MAKEVARS += HIP=1 HIP_ROOT=$(HIP_ROOT) HIP_PLATFORM=$(HIP_PLATFORM) HIP_CUDA=$(HIP_CUDA) \
                    $(if $(HIP_ARCH),HIP_ARCH="$(HIP_ARCH)")
  endif
endif

all: $(B)/ias15_ref$(EXE) $(B)/ias15_cft$(EXE) $(B)/check_dropin$(EXE) \
     $(B)/libcft_rebound.a $(B)/libcft_ias15.a

.PHONY: all third-party libcft librebound constants check clean

third-party:
	bash tools/fetch_third_party.sh

$(CFTLIB) libcft:
	$(MAKE) -C $(CFT) CC=$(CC) PYTHON=$(PYTHON) $(CFT_MAKEVARS) libcft.a

$(B)/rebound/%.o: $(REB)/%.c
	@mkdir -p $(B)/rebound
	$(CC) -c $(REB_CFLAGS) -I$(REB) -o $@ $<

$(B)/rebound/whfast512_stub.o: ref/whfast512_stub.c
	@mkdir -p $(B)/rebound
	$(CC) -c $(REB_CFLAGS) -I$(REB) -o $@ $<

$(B)/librebound.a: $(REB_OBJ)
	ar rcs $@ $^

librebound: $(B)/librebound.a

$(B)/ias15_ref$(EXE): ref/ias15_ref.c src/hexfloat.h $(B)/librebound.a
	$(CC) $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -I$(REB) -o $@ ref/ias15_ref.c $(B)/librebound.a $(LIBS)

# The packaging layer: what a REBOUND program links against. It runs
# the ias15_cft program, so it is told where `make install` will put it;
# if that file is not there at run time the search falls through to
# $CFT_REBOUND_IAS15 and then to PATH (src/cft_rebound_run.c).
$(B)/bindir.stamp: FORCE
	@mkdir -p $(B)
	@echo '$(BINDIR)' | cmp -s - $@ || echo '$(BINDIR)' > $@
.PHONY: FORCE
FORCE:

$(B)/cft_rebound_run.o: src/cft_rebound_run.c include/cft_rebound.h src/hexfloat.h src/cft_supported.h $(B)/bindir.stamp
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Iinclude -Isrc -I$(REB) -I$(CFT)/include \
	    -DCFT_REBOUND_IAS15_DEFAULT='"$(BINDIR)/ias15_cft$(EXE)"' -o $@ src/cft_rebound_run.c

$(B)/libcft_rebound.a: $(B)/cft_rebound_run.o $(B)/cft_supported.o
	ar rcs $@ $^

# Gravity's sequencer programs, compiled into both halves (ROUND3.md, P1);
# the header's rule is beside `programs`, below, where $(ASM) is defined.
GRAVITY_CFTA := $(wildcard programs/gravity-*.cfta)
GRAVITY_H    := $(B)/gravity_images.h src/cft_embedded_image.h

$(B)/ias15_cft$(EXE): src/ias15_cft.c src/ias15_zonal.inc src/ias15_constants.h src/hexfloat.h src/ias15_limits.h $(B)/normalabs_images.h $(GRAVITY_H) $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) -Isrc -I$(B) -I$(CFT)/include -o $@ src/ias15_cft.c $(CFTLIB) $(LIBS)

# ---------------------------------------------------------------------
# The drop-in: the same engine compiled as a library (no main) plus the
# struct reb_integrator shim, and the gate that runs it against
# REBOUND's own ias15 inside one program.
#
# src/cft_ias15_fields.c holds the one definition of the archive's
# field descriptors, including the list the integrator registers with.
# It was a placeholder holding only its terminator until 2026-09-10,
# which meant REBOUND could resolve no cft_ field on read - see
# docs/VALIDATION.md entry 25.
DROPIN_OBJ := $(B)/ias15_cft_lib.o $(B)/reb_integrator_cft.o $(B)/cft_ias15_fields.o \
              $(B)/cft_supported.o

$(B)/ias15_cft_lib.o: src/ias15_cft.c src/ias15_zonal.inc src/ias15_constants.h src/ias15_engine.h src/ias15_limits.h $(B)/normalabs_images.h $(GRAVITY_H)
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) -DIAS15_CFT_LIBRARY -Isrc -I$(B) -I$(CFT)/include -o $@ src/ias15_cft.c

$(B)/reb_integrator_cft.o: src/reb_integrator_cft.c src/cft_ias15.h src/ias15_engine.h src/cft_supported.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/reb_integrator_cft.c

# The refusal list, as data: one table both entry points walk, so it
# goes into libcft_ias15.a and libcft_rebound.a alike.
$(B)/cft_supported.o: src/cft_supported.c src/cft_supported.h src/ias15_limits.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/cft_supported.c

# The header carries the macros the lists are generated from, so an
# added blob has to rebuild this object as well as cft_archive.o.
$(B)/cft_ias15_fields.o: src/cft_ias15_fields.c src/cft_ias15_fields.h src/cft_ias15.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/cft_ias15_fields.c

# cft_archive.c compiled for the library rather than for a gate: the
# gates add -Itests for the stub's header and this does not need it.
$(B)/cft_archive.o: src/cft_archive.c src/cft_archive.h src/cft_ias15_fields.h src/cft_ias15_state.h src/cft_ias15.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/cft_archive.c

# The drop-in as a library: the registered integrator, the engine and
# Simulationarchive support. Named for the integrator it registers,
# because libcft_rebound.a - the subprocess API - is a different
# library and one underscore is not enough to tell them apart in a
# lib directory.
$(B)/libcft_ias15.a: $(DROPIN_OBJ) $(B)/cft_archive.o
	ar rcs $@ $^

# ---- the same library as a shared object, for Python ----------------
# The .a and the .so are not the same objects: -fPIC, and the .so also
# carries src/cft_ias15_shared.c, whose constructor registers the
# integrator at dlopen time. A C caller linking the archive calls
# cft_ias15_register() itself and must not get a constructor doing it
# behind their back. (mingw runs the same __attribute__((constructor))
# from DllMain's DLL_PROCESS_ATTACH, so the DLL registers itself at
# LoadLibrary time exactly as the .so does at dlopen time.)
#
# $(REB_USEFLAGS) is deliberately absent, and on Windows that is
# load bearing rather than an omission: rebound.h marks its API
# __declspec(dllimport) unless BUILDINGLIBREBOUND is defined, and
# dllimport is precisely what this library wants - it imports those
# symbols from the librebound already in the process. Everything else
# in this Makefile links REBOUND statically and therefore does define
# it.
PIC_OBJ := $(B)/pic/ias15_cft_lib.o $(B)/pic/reb_integrator_cft.o \
           $(B)/pic/cft_ias15_fields.o $(B)/pic/cft_archive.o \
           $(B)/pic/cft_ias15_shared.o $(B)/pic/cft_supported.o
PICFLAGS := $(CSTD) $(CFLAGS) $(WARN) $(SHPIC) -Isrc -I$(CFT)/include -I$(REB)

$(B)/pic/ias15_cft_lib.o: src/ias15_cft.c src/ias15_zonal.inc src/ias15_constants.h src/ias15_engine.h src/ias15_limits.h $(B)/normalabs_images.h $(GRAVITY_H)
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -DIAS15_CFT_LIBRARY -I$(B) -o $@ src/ias15_cft.c

$(B)/pic/reb_integrator_cft.o: src/reb_integrator_cft.c src/cft_ias15.h src/ias15_engine.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -o $@ src/reb_integrator_cft.c

$(B)/pic/cft_ias15_fields.o: src/cft_ias15_fields.c src/cft_ias15_fields.h src/cft_ias15.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -o $@ src/cft_ias15_fields.c

$(B)/pic/cft_archive.o: src/cft_archive.c src/cft_archive.h src/cft_ias15_fields.h src/cft_ias15_state.h src/cft_ias15.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -o $@ src/cft_archive.c

$(B)/pic/cft_ias15_shared.o: src/cft_ias15_shared.c src/cft_ias15.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -o $@ src/cft_ias15_shared.c

$(B)/pic/cft_supported.o: src/cft_supported.c src/cft_supported.h src/ias15_limits.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -o $@ src/cft_supported.c

# On POSIX, REBOUND's symbols are deliberately NOT linked: the caller's
# process already has librebound loaded and python/cft_rebound.py
# promotes it to RTLD_GLOBAL before opening this. Linking a second copy
# would give a second integrator list and a name that can never be
# selected.
#
# libcft.a will NOT go into a POSIX shared object: backend_xrt.o is C++
# and carries relocations a shared object cannot use (R_X86_64_PC32
# against a GLIBCXX symbol). cft-fp256 builds a proper shared libcft
# from its own PIC objects, so use that. Two rpaths: $$ORIGIN for an
# installed layout where the two sit together, and the tree for running
# it in place without an install.
#
# On Windows both halves of that are different, and the Windows rule is
# below. Only build cft-fp256's shared library where we actually want
# one - on Windows CFT_SHLIB is the archive, which already has a rule.
ifneq ($(CFT_SHLIB),$(CFTLIB))
$(CFT_SHLIB):
	$(MAKE) -C $(CFT) CC=$(CC) PYTHON=$(PYTHON) $(CFT_MAKEVARS) $(notdir $(CFT_SHLIB))
endif

ifeq ($(OS),Windows_NT)
# A Windows DLL may not carry undefined symbols. So this one LINKS the
# librebound it is going to be loaded into - the wheel's .pyd, asked of
# the interpreter that will load it - and the import table then names
# that one exact file. That is the whole cost of the Windows port and
# it is not hideable: a DLL built for cp312 cannot be loaded by cp313,
# and python/cft_rebound.py says so by name rather than letting the
# loader say "The specified module could not be found."
#
# Nothing has to be on the link line for libcft: unlike the .so, the
# DLL swallows the archive whole (CFT_SHLINK above says why).
#
# The import table is asserted, not printed. Printed, it named
# libwinpthread-1.dll on every build for seven weeks and nothing read
# it. MSYS2's runtime DLLs are the lib*.dll family and a Python process
# finds none of them, so such a library is refused and removed here,
# not at ctypes.CDLL with a message that names the wrong file.
$(B)/$(SHLIB): $(PIC_OBJ) $(CFT_SHLIB) | check-rebound-match
	@pyd=`$(PYTHON) -c 'import rebound;print(rebound.__libpath__)' 2>/dev/null`; \
	if [ -z "$$pyd" ]; then \
	    echo "$@: no rebound importable by $(PYTHON)." >&2; \
	    echo "  A Windows DLL may not carry undefined symbols, so this library must" >&2; \
	    echo "  link the librebound it will be loaded into, and there is none to link." >&2; \
	    echo "      $(PYTHON) -m pip install rebound" >&2; \
	    exit 1; \
	fi; \
	echo "$(CC) $(SHFLAGS) -o $@ \$$(PIC_OBJ) $(CFT_SHLINK) $$pyd $(LIBS) $(SHSTATIC)"; \
	$(CC) $(SHFLAGS) -o $@ $(PIC_OBJ) $(CFT_SHLINK) "$$pyd" $(LIBS) $(SHSTATIC)
	@command -v objdump >/dev/null 2>&1 || { echo "$@: no objdump, so its import table is NOT checked"; exit 0; }; \
	imp=`objdump -p $@ | sed -n 's/.*DLL Name: //p' | sort -u | tr '\n' ' '`; \
	echo "$@ imports: $$imp"; \
	bad=`echo "$$imp" | tr ' ' '\n' | grep -i -E '^lib.*\.dll$$' | tr '\n' ' '`; \
	if [ -n "$$bad" ]; then \
	    echo "$@: REFUSING - it imports $$bad- no Python process can find that (docs/PYTHON.md, \"The Windows rule\")" >&2; \
	    rm -f $@; exit 1; \
	fi
else
$(B)/$(SHLIB): $(PIC_OBJ) $(CFT_SHLIB) | check-rebound-match
	$(CC) $(SHFLAGS) -o $@ $(PIC_OBJ) $(CFT_SHLINK) $(SHRPATH) $(LIBS)
endif

.PHONY: python-lib check-python
python-lib: $(B)/$(SHLIB)
	@echo
	@echo "built $(B)/$(SHLIB). From Python:"
	@echo "    import cft_rebound"
	@echo "    cft_rebound.load(r'$(SHLIB_PATH)')"

# The library is compiled against the pinned REBOUND headers and is
# loaded into a process running the WHEEL's REBOUND, so the two must be
# the same source or struct reb_simulation is laid out differently on
# each side of the call and nothing says so. They are identical for
# rebound 5.1.1 and the pinned bdfda4bd - on Windows too, where the
# wheel is MSVC-built and the library is not - and this refuses to
# build a library that would be quietly wrong if that ever stops being
# true. It is an order-only prerequisite of $(B)/$(SHLIB), so it now
# runs BEFORE the link rather than beside it, and check-python gets it
# too; order-only means a phony prerequisite cannot force a relink.
#
# chr(92) rather than a backslash literal: the path the wheel reports
# on Windows is C:\..., and MSYS cmp is handed it as a shell word.
.PHONY: check-rebound-match
check-rebound-match:
	@w=`$(PYTHON) -c 'import rebound,os;print(os.path.dirname(rebound.__libpath__).replace(chr(92),"/"))' 2>/dev/null`; \
	if [ -z "$$w" ]; then \
	    echo "check-rebound-match: no rebound wheel importable by $(PYTHON); skipping the comparison"; \
	elif [ ! -f "$$w/src/rebound.h" ]; then \
	    echo "check-rebound-match: the wheel at $$w ships no src/rebound.h; cannot compare"; \
	elif cmp -s $(REB)/rebound.h "$$w/src/rebound.h"; then \
	    echo "check-rebound-match: the pinned REBOUND and the installed wheel are the same headers"; \
	else \
	    echo "check-rebound-match: REFUSING - $(REB)/rebound.h differs from the wheel's"; \
	    echo "  $$w/src/rebound.h"; \
	    echo "  A library built against one and loaded into the other shares a struct"; \
	    echo "  layout it may not have. Re-pin third_party/MANIFEST to the wheel's"; \
	    echo "  REBOUND, or install the wheel matching the pin."; \
	    exit 1; \
	fi

# The Python gate: REBOUND's own ias15 against the registered one,
# from Python, on the same problem and the same fixed step. The same 13
# values on every platform, and the same claim: at binary64 they are
# identical or the gate has failed. $(SHLIB_PATH) rather than
# $(CURDIR)/... because MSYS make hands recipes /c/... and a native
# Windows Python cannot open that.
check-python: $(B)/$(SHLIB)
	$(PYTHON) python/example_equivalence.py --library '$(SHLIB_PATH)'

# The control for it. The same run at binary128 MUST differ: if a
# format change leaves 13 values identical then the format is being
# accepted and ignored, and check-python above is comparing binary64
# against binary64 and passing against itself. --expect-differ inverts
# the exit status, so this fails loudly when nothing moved.
.PHONY: check-python-control
check-python-control: $(B)/$(SHLIB)
	$(PYTHON) python/example_equivalence.py --library '$(SHLIB_PATH)' \
	    --format fp128 --expect-differ

.PHONY: dropin
dropin: $(B)/libcft_ias15.a $(B)/check_dropin$(EXE)

# The gate is a file per topic now: main() and the refusals in
# check_dropin.c, the scaffolding in dropin_common.c, the cases in
# cases_*.c. A parcel adds its file here and nowhere else -
# tools/dropin_cases.h says how.
CASES_SRC := tools/check_dropin.c tools/dropin_common.c \
             tools/cases_core.c tools/cases_wide.c \
             tools/cases_forces.c tools/cases_modes.c \
             tools/cases_pairs.c tools/cases_seam.c tools/cases_corrector.c \
             tools/cases_gravity.c tools/cases_scalar.c

$(B)/check_dropin$(EXE): $(CASES_SRC) tools/dropin_cases.h src/cft_supported.h $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) \
	      -o $@ $(CASES_SRC) $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB) $(LIBS)

constants:
	$(PYTHON) tools/gen_constants.py

# --- Simulationarchive support and its gates -------------------------
# The archive module is compiled against BOTH upstreams: REBOUND for the
# field descriptors and the loader, libcft for the wide formats.
ARCHIVE_SRC := src/cft_archive.c src/cft_archive.h src/cft_ias15_state.h src/cft_ias15_fields.c src/cft_ias15_fields.h
ARCH_CFLAGS := $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -I$(REB) -I$(CFT)/include -Isrc -Itests
ARCHIVE_GATES := $(B)/gate_restart$(EXE) $(B)/gate_write$(EXE) \
                 $(B)/gate_stock$(EXE) $(B)/gate_promote$(EXE) \
                 $(B)/gate_real$(EXE) $(B)/gate_subprocess$(EXE)

# The subprocess API's gate. It links libcft_rebound.a rather than the
# sources, so what it exercises is what a caller linking the shipped
# archive would get, and it needs the standalone program at run time
# because spawning it is the whole of what that API does.
$(B)/gate_subprocess$(EXE): tests/gate_subprocess.c src/cft_supported.h \
                            $(B)/libcft_rebound.a $(B)/librebound.a $(CFTLIB) \
                            $(B)/ias15_cft$(EXE)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Iinclude -Isrc -I$(REB) -I$(CFT)/include \
	      -o $@ tests/gate_subprocess.c $(B)/libcft_rebound.a $(B)/librebound.a $(CFTLIB) $(LIBS)

$(B)/gate_restart$(EXE): tests/gate_restart.c tests/cft_shim_stub.c tests/cft_shim_stub.h $(ARCHIVE_SRC) $(B)/librebound.a $(CFTLIB)
	$(CC) $(ARCH_CFLAGS) -o $@ tests/gate_restart.c tests/cft_shim_stub.c src/cft_archive.c src/cft_ias15_fields.c $(B)/librebound.a $(CFTLIB) $(LIBS)

$(B)/gate_write$(EXE): tests/gate_write.c tests/cft_shim_stub.c tests/cft_shim_stub.h $(ARCHIVE_SRC) $(B)/librebound.a $(CFTLIB)
	$(CC) $(ARCH_CFLAGS) -o $@ tests/gate_write.c tests/cft_shim_stub.c src/cft_archive.c src/cft_ias15_fields.c $(B)/librebound.a $(CFTLIB) $(LIBS)

$(B)/gate_promote$(EXE): tests/gate_promote.c tests/cft_shim_stub.c tests/cft_shim_stub.h $(ARCHIVE_SRC) $(B)/librebound.a $(CFTLIB)
	$(CC) $(ARCH_CFLAGS) -o $@ tests/gate_promote.c tests/cft_shim_stub.c src/cft_archive.c src/cft_ias15_fields.c $(B)/librebound.a $(CFTLIB) $(LIBS)

# Deliberately links the pinned upstream librebound and NOTHING of
# cft-rebound's: this program is a stock REBOUND reader.
$(B)/gate_stock$(EXE): tests/gate_stock.c $(B)/librebound.a
	$(CC) $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -I$(REB) -o $@ tests/gate_stock.c $(B)/librebound.a $(LIBS)

# The only target that links the REAL shim and the archive together.
# The four gates above link tests/cft_shim_stub.c, whose own header says
# to re-point them once the real integrator exists; this gate is that
# re-pointing, kept separate so the stub gates go on proving descriptor
# completeness - their step touches 48 of the 50 blobs and their
# snapshot compares all 50, and 20 steps of a real IAS15 need not.
# tests/cft_shim_stub.h has the exact division.
$(B)/gate_real$(EXE): tests/gate_real.c $(ARCHIVE_SRC) $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB)
	$(CC) $(ARCH_CFLAGS) -o $@ tests/gate_real.c src/cft_archive.c $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB) $(LIBS)

.PHONY: archive check-archive
archive: $(ARCHIVE_GATES)

check-archive: $(ARCHIVE_GATES)
	$(PYTHON) tools/check_archive.py --build $(B)

# The sequencer programs: generated as text, assembled by the pinned
# clone's own assembler (built here from its sources).
ASM := $(CFT)/cft-asm$(EXE)
# After libcft.a: both are sub-makes in the same clone, and under -j two of
# them building libcft.a at once is a race (the lead's 2026-11-02 trap).
# `all` reaches this rule now, through gravity's header.
$(ASM): $(CFTLIB)
	$(MAKE) -C $(CFT) CC=$(CC) PYTHON=$(PYTHON) $(CFT_MAKEVARS) cft-asm$(EXE)

# Gravity's sequencer programs, compiled into both halves (ROUND3.md, P1):
# the library half loads nothing from disk. Assembled by the pinned
# cft-asm from the checked-in programs/gravity-*.cfta, which `make
# programs` writes. The header is rewritten only when its bytes change, so
# that `programs` rewriting every .cfta does not relink the engine on each
# suite run (tools/embed_programs.py says why that matters on Windows).
$(B)/gravity_images.h: tools/embed_programs.py $(GRAVITY_CFTA) $(ASM)
	@mkdir -p $(B)
	$(PYTHON) tools/embed_programs.py --asm $(ASM) --name gravity --out $@ $(GRAVITY_CFTA)

.PHONY: programs
programs: $(ASM)
	$(PYTHON) tools/gen_programs.py
	@mkdir -p programs/out
	@for f in programs/*.cfta; do b=$$(basename $$f .cfta); $(ASM) $$f -o programs/out/$$b.cftp || exit 1; done
	@echo "assembled $$(ls programs/out/*.cftp | wc -l) programs"

# Program images compiled INTO the engine. The library half loads nothing
# from disk - load_one() is the standalone program's alone - so a route
# whose new half is a program run carries its image in the binary, one
# generated header a family, which both halves compile in
# (tools/embed_programs.py, src/cft_embedded_image.h).
#
# This family is cft-fp256's normalabs-<fmt>, the corrector's exit test on
# the device (ROUND3.md, P2; pc_error() in src/ias15_cft.c). It is
# assembled from the PINNED tree by the pinned assembler, so the pin
# governs it and no copy drifts, and tools/check_corrector.py holds the
# image the engine actually loads (cft_program_digest) to the SHA-256 the
# pinned tree's own programs/MANIFEST gives it.
CFTPROGS      := third_party/cft-fp256/programs
NORMALABS_SRC := $(foreach f,fp64 fp128 fp256,$(CFTPROGS)/normalabs-$(f).cfta)
$(B)/normalabs_images.h: tools/embed_programs.py $(NORMALABS_SRC) $(ASM)
	@mkdir -p $(B)
	$(PYTHON) tools/embed_programs.py --asm $(ASM) --name normalabs --out $@ $(NORMALABS_SRC)

# ---------------------------------------------------------------------
# Suite legs go through tools/gate_cache.py. `check` runs them with
# --mode full: it never skips and it writes the stamps. `check-quick`
# runs them with --mode quick, which skips a leg whose inputs are
# byte-for-byte what they were when that leg last passed.
#
# The key is over the ARTIFACTS a leg executes, its own script and the
# data it reads - never over source files, because a source list has to
# enumerate every header and flag and the cost of missing one is
# skipping a leg that would have failed. The header of gate_cache.py
# has the reasoning.
#
# $(ART) is in every key: a pass on the software backend must never
# satisfy a run against a card. So is every other CFT_REBOUND_
# variable, by prefix and not by list. The door's route variables
# choose which route a leg exercises (src/ias15_cft.c, engine_routes()),
# and a pass with every route hidden must not satisfy a default run.
GATE_CACHE  := $(PYTHON) tools/gate_cache.py --cache-dir $(B)/.gate-cache \
               --env CFT_REBOUND_ARTIFACT --env-prefix CFT_REBOUND_
PROBLEMS    := $(wildcard data/problems/*)
# The program images the legs run, named from their sources. `programs`
# assembles each programs/*.cfta into programs/out/<name>.cftp, but only
# after make has read this file, so a wildcard over programs/out alone
# was empty in a fresh clone, and the legs that list the images
# ($(IN_PROGS)) keyed on none of them (measured in round 5's first
# suite). Named from the sources, the list is the same whether
# programs/out exists yet or not. Whatever programs/out already holds is
# listed too, as that wildcard listed it: check_gravity holds every
# gravity image it finds there, so an image whose source has gone is part
# of its verdict. A source that `programs` writes for the first time (a
# program new in tools/gen_programs.py) is listed from the next make on.
PROGRAMS    := $(sort $(patsubst programs/%.cfta,programs/out/%.cftp,$(wildcard programs/*.cfta)) \
                      $(wildcard programs/out/*.cftp))
RECORDS     := $(wildcard results/raw/*.rec)

IN_DROPIN   := --input $(B)/check_dropin$(EXE)
IN_PROGRAM  := --input $(B)/ias15_cft$(EXE)
IN_REF      := --input $(B)/ias15_ref$(EXE)
IN_PROB     := $(patsubst %,--input %,$(PROBLEMS))
IN_PROGS    := $(patsubst %,--input %,$(PROGRAMS))
IN_RECS     := $(patsubst %,--input %,$(RECORDS))
IN_ARCHIVE  := --input $(B)/gate_restart$(EXE) --input $(B)/gate_write$(EXE) \
               --input $(B)/gate_stock$(EXE) --input $(B)/gate_promote$(EXE)
IN_REAL     := --input $(B)/gate_real$(EXE)
IN_SUBPROC  := --input $(B)/gate_subprocess$(EXE) --input $(B)/ias15_cft$(EXE)
# the corrector's gate runs the program on every problem, on the program
# engine too, generates an ensemble, and reads the pinned MANIFEST
IN_CORRECTOR := --input tools/check_corrector.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) \
                --input tools/make_ensemble.py --input $(CFTPROGS)/MANIFEST
# tools/check_replay.py: the program, the program images its program-engine
# cases run, and the tool it imports to cut a test particle out of its
# wave; it writes its own waves (ROUND4.md, P3)
IN_REPLAY   := --input tools/check_replay.py --input tools/extract_particle.py $(IN_PROGRAM) $(IN_PROGS)
# tools/check_gravity.py: the program, the problems, the program images it
# holds the binary to, the two generators it runs and the header it reads
# the body cap from; --library adds the drop-in gate it runs.
IN_GRAVITY  := --input tools/check_gravity.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) \
               --input tools/make_nbody.py --input tools/make_ensemble.py --input src/ias15_limits.h
# tools/check_scalar.py: the program on every problem, on both engines, the
# ensembles it generates, the source it recounts the sites from (with the
# zonal term's include, which the recount reads too), and with --library the
# drop-in gate and gate_real, the route forced on.
IN_SCALAR   := --input tools/check_scalar.py --input src/ias15_cft.c --input src/ias15_zonal.inc $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) \
               --input tools/make_ensemble.py $(IN_DROPIN) $(IN_REAL)
# tools/check_zonal.py (round 4's P1): the program, its selftest included, on
# the zonal problem it derives and holds, on the program engine too, and the
# golden model it runs from the pinned cft-fp256 clone (whose pin is MANIFEST).
GOLDEN_PY   := $(wildcard third_party/cft-fp256/python/cft_golden/*.py)
IN_ZONAL    := --input tools/check_zonal.py $(IN_PROGRAM) --input data/problems/zonal_saturn.txt $(IN_PROGS) \
               $(patsubst %,--input %,$(GOLDEN_PY)) --input third_party/MANIFEST
# tools/check_seams.py: gate_real through check_checkpoint's runner, and
# the program on the problems, the program images and an ensemble
IN_SEAMS    := --input tools/check_seams.py --input tools/check_checkpoint.py $(IN_REAL) \
               $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) --input tools/make_ensemble.py
# tools/check_ring.py (round 4's P2): the program on the rings the generator
# writes, on both engines and so the program images, REBOUND's own IAS15 on
# one wave, and the helper it imports; it reads no file of data/problems/
IN_RING     := --input tools/check_ring.py --input tools/make_ring.py --input tools/check_equivalence.py \
               $(IN_PROGRAM) $(IN_REF) $(IN_PROGS)
# tools/check_ring_seams.py (round 4's seams, the lead's): the program on the
# rings the generator writes, with and without their zonal lines, on both
# engines and so the program images; a particle replayed out of a wave
# through P3's tool; P4's driver, tools/compare_run.py, on a ring; and the
# modules those import; it reads no file of data/problems/
IN_RSEAMS   := --input tools/check_ring_seams.py --input tools/make_ring.py --input tools/check_ensemble.py \
               --input tools/check_equivalence.py --input tools/extract_particle.py --input tools/check_replay.py \
               --input tools/compare_run.py --input tools/check_records.py $(IN_PROGRAM) $(IN_PROGS)
# tools/check_compare_run.py: the driver it holds and every script the
# driver imports or runs, the program the driver reruns, REBOUND's own
# program (whose record the driver must refuse) on the problem it reads, and
# the program images for the program-engine case (ROUND4.md, P4). P3's two
# files are inputs whenever they exist: the driver's verdict turns on their
# presence (its particle place), so a pass from before they landed must not
# satisfy a run after.
IN_COMPARE  := --input tools/check_compare_run.py --input tools/compare_run.py --input tools/check_ensemble.py \
               --input tools/check_records.py --input tools/check_equivalence.py --input tools/compare_formats.py \
               --input tools/oracle.py $(IN_PROGRAM) $(IN_REF) --input data/problems/kepler.txt $(IN_PROGS) \
               $(patsubst %,--input %,$(wildcard tools/extract_particle.py tools/check_replay.py))
# tools/check_reversal.py (VALIDATION entry 52): tools/reversal.py and the
# two files it runs or imports, which the gate also copies and mutates; the
# program on zonal_saturn, whose `active 2` the gate makes `active 1`; and
# the committed reversal records and problems the gate rebuilds
IN_REVERSAL := --input tools/check_reversal.py --input tools/reversal.py --input tools/state_to_problem.py \
               --input tools/check_equivalence.py $(IN_PROGRAM) --input data/problems/zonal_saturn.txt \
               $(patsubst %,--input %,$(wildcard results/horizon/reversal/*))

# The pins, checked at build time and not only at fetch time: a clone
# under third_party/ that is not at its MANIFEST commit fails every tier
# before a single gate runs, and so does a README whose one cft-fp256
# row does not begin its commit with the one MANIFEST pins. Not through the gate cache
# - it is instant, and a cached pass would be exactly the wrong thing.
.PHONY: check-pins
check-pins:
	@bash tools/check_pins.sh

check: check-pins all programs $(ARCHIVE_GATES)
	$(GATE_CACHE) --leg constants --mode full --input tools/gen_constants.py --input third_party/rebound/src/integrator_ias15.c -- $(PYTHON) tools/gen_constants.py --no-write
	$(GATE_CACHE) --leg check_dropin --mode full $(IN_DROPIN) -- $(B)/check_dropin$(EXE)
	$(GATE_CACHE) --leg check_dropin_wide --mode full $(IN_DROPIN) -- $(B)/check_dropin$(EXE) --wide
	$(GATE_CACHE) --leg check_equivalence --mode full --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_equivalence.py --build $(B)
	$(GATE_CACHE) --leg check_program_engine --mode full --input tools/check_program_engine.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_program_engine.py --build $(B)
	$(GATE_CACHE) --leg check_door --mode full --input tools/check_door.py $(IN_PROGRAM) $(IN_PROB) -- $(PYTHON) tools/check_door.py --build $(B)
	$(GATE_CACHE) --leg check_corrector --mode full $(IN_CORRECTOR) $(IN_DROPIN) -- $(PYTHON) tools/check_corrector.py --build $(B) --library
	$(GATE_CACHE) --leg check_replay --mode full $(IN_REPLAY) -- $(PYTHON) tools/check_replay.py --build $(B)
	$(GATE_CACHE) --leg check_gravity --mode full $(IN_GRAVITY) $(IN_DROPIN) -- $(PYTHON) tools/check_gravity.py --build $(B) --library
	$(GATE_CACHE) --leg check_zonal --mode full $(IN_ZONAL) -- $(PYTHON) tools/check_zonal.py --build $(B)
	$(GATE_CACHE) --leg check_scalar --mode full $(IN_SCALAR) -- $(PYTHON) tools/check_scalar.py --build $(B) --library
	$(GATE_CACHE) --leg check_records --mode full --input tools/check_records.py --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_PROB) $(IN_RECS) -- $(PYTHON) tools/check_records.py --build $(B)
	$(GATE_CACHE) --leg check_ensemble --mode full --input tools/check_ensemble.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_ensemble.py --build $(B)
	$(GATE_CACHE) --leg check_ring --mode full $(IN_RING) -- $(PYTHON) tools/check_ring.py --build $(B)
	$(GATE_CACHE) --leg check_compare_run --mode full $(IN_COMPARE) -- $(PYTHON) tools/check_compare_run.py --build $(B)
	$(GATE_CACHE) --leg check_archive --mode full --input tools/check_archive.py $(IN_ARCHIVE) -- $(PYTHON) tools/check_archive.py --build $(B)
	$(GATE_CACHE) --leg check_bodycount --mode full --input tools/check_bodycount.py --input tools/check_equivalence.py --input tools/make_nbody.py --input src/ias15_limits.h $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_bodycount.py --build $(B)
	$(GATE_CACHE) --leg gate_real_fp64 --mode full $(IN_REAL) -- $(B)/gate_real$(EXE) --fp64
	$(GATE_CACHE) --leg gate_real_fp128 --mode full $(IN_REAL) -- $(B)/gate_real$(EXE) --fp128
	$(GATE_CACHE) --leg gate_real_fp256 --mode full $(IN_REAL) -- $(B)/gate_real$(EXE) --fp256
	$(GATE_CACHE) --leg gate_subprocess --mode full $(IN_SUBPROC) -- $(B)/gate_subprocess$(EXE) --build $(B)
	$(GATE_CACHE) --leg check_checkpoint --mode full --input tools/check_checkpoint.py $(IN_REAL) -- $(PYTHON) tools/check_checkpoint.py --build $(B)
	$(GATE_CACHE) --leg check_seams --mode full $(IN_SEAMS) -- $(PYTHON) tools/check_seams.py --build $(B)
	$(GATE_CACHE) --leg check_ring_seams --mode full $(IN_RSEAMS) -- $(PYTHON) tools/check_ring_seams.py --build $(B)
	$(GATE_CACHE) --leg check_reversal --mode full $(IN_REVERSAL) -- $(PYTHON) tools/check_reversal.py --build $(B)

# the same gates at binary64 only, in a few minutes
.PHONY: check-quick
# What the cache holds, and whether each stamp came from a full run or
# a quick one. Reads nothing else and changes nothing.
.PHONY: gate-cache-report
gate-cache-report:
	@$(PYTHON) tools/gate_cache.py --report --cache-dir $(B)/.gate-cache

check-quick: check-pins all programs $(ARCHIVE_GATES)
	$(GATE_CACHE) --leg constants --mode quick --input tools/gen_constants.py --input third_party/rebound/src/integrator_ias15.c -- $(PYTHON) tools/gen_constants.py --no-write
	$(GATE_CACHE) --leg check_dropin --mode quick $(IN_DROPIN) -- $(B)/check_dropin$(EXE)
	$(GATE_CACHE) --leg check_dropin_wide --mode quick $(IN_DROPIN) -- $(B)/check_dropin$(EXE) --wide
	$(GATE_CACHE) --leg check_equivalence --mode quick --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_equivalence.py --build $(B) --quick
	$(GATE_CACHE) --leg check_program_engine --mode quick --input tools/check_program_engine.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_program_engine.py --build $(B) --formats fp64
	$(GATE_CACHE) --leg check_door --mode quick --input tools/check_door.py $(IN_PROGRAM) $(IN_PROB) -- $(PYTHON) tools/check_door.py --build $(B) --formats fp64
	$(GATE_CACHE) --leg check_corrector --mode quick $(IN_CORRECTOR) $(IN_DROPIN) -- $(PYTHON) tools/check_corrector.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_replay --mode quick $(IN_REPLAY) -- $(PYTHON) tools/check_replay.py --build $(B) --quick
	$(GATE_CACHE) --leg check_gravity --mode quick $(IN_GRAVITY) $(IN_DROPIN) -- $(PYTHON) tools/check_gravity.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_zonal --mode quick $(IN_ZONAL) -- $(PYTHON) tools/check_zonal.py --build $(B) --quick
	$(GATE_CACHE) --leg check_scalar --mode quick $(IN_SCALAR) -- $(PYTHON) tools/check_scalar.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_records --mode quick --input tools/check_records.py --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_PROB) $(IN_RECS) -- $(PYTHON) tools/check_records.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ensemble --mode quick --input tools/check_ensemble.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_ensemble.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ring --mode quick $(IN_RING) -- $(PYTHON) tools/check_ring.py --build $(B) --quick
	$(GATE_CACHE) --leg check_compare_run --mode quick $(IN_COMPARE) -- $(PYTHON) tools/check_compare_run.py --build $(B) --quick
	$(GATE_CACHE) --leg check_archive --mode quick --input tools/check_archive.py $(IN_ARCHIVE) -- $(PYTHON) tools/check_archive.py --build $(B)
	$(GATE_CACHE) --leg check_bodycount --mode quick --input tools/check_bodycount.py --input tools/check_equivalence.py --input tools/make_nbody.py --input src/ias15_limits.h $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_bodycount.py --build $(B) --quick
	$(GATE_CACHE) --leg gate_real_fp64 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp64
	$(GATE_CACHE) --leg gate_real_fp128 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp128
	$(GATE_CACHE) --leg gate_real_fp256 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp256
	$(GATE_CACHE) --leg gate_subprocess --mode quick $(IN_SUBPROC) -- $(B)/gate_subprocess$(EXE) --build $(B)
	$(GATE_CACHE) --leg check_checkpoint --mode quick --input tools/check_checkpoint.py $(IN_REAL) -- $(PYTHON) tools/check_checkpoint.py --build $(B) --quick
	$(GATE_CACHE) --leg check_seams --mode quick $(IN_SEAMS) -- $(PYTHON) tools/check_seams.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ring_seams --mode quick $(IN_RSEAMS) -- $(PYTHON) tools/check_ring_seams.py --build $(B) --quick
	$(GATE_CACHE) --leg check_reversal --mode quick $(IN_REVERSAL) -- $(PYTHON) tools/check_reversal.py --build $(B) --quick
	@$(PYTHON) tools/gate_cache.py --report --cache-dir $(B)/.gate-cache

# Every case touched, none of them thoroughly - the step counts are
# scaled (tools/dropin_cases.h, light_steps). Measured on check_dropin:
# 159 s against 585 s, the same 130 cases, zero failures. The one case
# whose step count IS the case opts out by name; see case_collision.
#
# Its check_dropin legs are named apart from check-quick's on purpose.
# The cache keys on the command, so the two could never cross-satisfy -
# but sharing a stamp file would make each overwrite the other's and
# neither would ever skip.
.PHONY: check-light
check-light: check-pins all programs $(ARCHIVE_GATES)
	$(GATE_CACHE) --leg constants --mode quick --input tools/gen_constants.py --input third_party/rebound/src/integrator_ias15.c -- $(PYTHON) tools/gen_constants.py --no-write
	$(GATE_CACHE) --leg check_dropin_light --mode quick $(IN_DROPIN) -- $(B)/check_dropin$(EXE) --light
	$(GATE_CACHE) --leg check_dropin_wide_light --mode quick $(IN_DROPIN) -- $(B)/check_dropin$(EXE) --wide --light
	$(GATE_CACHE) --leg check_equivalence --mode quick --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_equivalence.py --build $(B) --quick
	$(GATE_CACHE) --leg check_program_engine --mode quick --input tools/check_program_engine.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_program_engine.py --build $(B) --formats fp64
	$(GATE_CACHE) --leg check_door --mode quick --input tools/check_door.py $(IN_PROGRAM) $(IN_PROB) -- $(PYTHON) tools/check_door.py --build $(B) --formats fp64
	$(GATE_CACHE) --leg check_corrector --mode quick $(IN_CORRECTOR) $(IN_DROPIN) -- $(PYTHON) tools/check_corrector.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_replay --mode quick $(IN_REPLAY) -- $(PYTHON) tools/check_replay.py --build $(B) --quick
	$(GATE_CACHE) --leg check_gravity --mode quick $(IN_GRAVITY) $(IN_DROPIN) -- $(PYTHON) tools/check_gravity.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_zonal --mode quick $(IN_ZONAL) -- $(PYTHON) tools/check_zonal.py --build $(B) --quick
	$(GATE_CACHE) --leg check_scalar --mode quick $(IN_SCALAR) -- $(PYTHON) tools/check_scalar.py --build $(B) --quick --library
	$(GATE_CACHE) --leg check_records --mode quick --input tools/check_records.py --input tools/check_equivalence.py $(IN_PROGRAM) $(IN_PROB) $(IN_RECS) -- $(PYTHON) tools/check_records.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ensemble --mode quick --input tools/check_ensemble.py --input tools/check_equivalence.py --input tools/make_ensemble.py $(IN_PROGRAM) $(IN_PROB) $(IN_PROGS) -- $(PYTHON) tools/check_ensemble.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ring --mode quick $(IN_RING) -- $(PYTHON) tools/check_ring.py --build $(B) --quick
	$(GATE_CACHE) --leg check_compare_run --mode quick $(IN_COMPARE) -- $(PYTHON) tools/check_compare_run.py --build $(B) --quick
	$(GATE_CACHE) --leg check_archive --mode quick --input tools/check_archive.py $(IN_ARCHIVE) -- $(PYTHON) tools/check_archive.py --build $(B)
	$(GATE_CACHE) --leg check_bodycount --mode quick --input tools/check_bodycount.py --input tools/check_equivalence.py --input tools/make_nbody.py --input src/ias15_limits.h $(IN_PROGRAM) $(IN_REF) $(IN_PROB) -- $(PYTHON) tools/check_bodycount.py --build $(B) --quick
	$(GATE_CACHE) --leg gate_real_fp64 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp64
	$(GATE_CACHE) --leg gate_real_fp128 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp128
	$(GATE_CACHE) --leg gate_real_fp256 --mode quick $(IN_REAL) -- $(B)/gate_real$(EXE) --fp256
	$(GATE_CACHE) --leg gate_subprocess --mode quick $(IN_SUBPROC) -- $(B)/gate_subprocess$(EXE) --build $(B)
	$(GATE_CACHE) --leg check_checkpoint --mode quick --input tools/check_checkpoint.py $(IN_REAL) -- $(PYTHON) tools/check_checkpoint.py --build $(B) --quick
	$(GATE_CACHE) --leg check_seams --mode quick $(IN_SEAMS) -- $(PYTHON) tools/check_seams.py --build $(B) --quick
	$(GATE_CACHE) --leg check_ring_seams --mode quick $(IN_RSEAMS) -- $(PYTHON) tools/check_ring_seams.py --build $(B) --quick
	$(GATE_CACHE) --leg check_reversal --mode quick $(IN_REVERSAL) -- $(PYTHON) tools/check_reversal.py --build $(B) --quick
	@$(PYTHON) tools/gate_cache.py --report --cache-dir $(B)/.gate-cache

# This tree on a GPU as it is (ROUND5.md, P1; ROADMAP's GPU half, step 2).
# A GPU host's target, made with HIP=1, and in none of the three tiers
# above: on a host with no hip:0, the desktop among them, it refuses by
# name. tools/check_gpu.py runs the witness's own controls, the routes on
# the device and the GPU seam test, then the record-making legs in their
# quick forms, their commands read from check-quick's recipe above
# (check_dropin from check-light's), with CFT_REBOUND_ARTIFACT=$(GPU_ARTIFACT);
# every process a leg starts goes through its witness, which holds the
# device line and the software backend's record from the same arguments.
# Not through the gate cache: a GPU run is its own evidence, and a cached
# pass would carry no device line. GPU_LEGS names a subset (and then the
# controls, routes and seam do not run), GPU_SKIP leaves legs out,
# and GPU_EXPECT is text the device line must contain or nothing runs (in
# WSL cft-gpu, "RTX 5060 Ti": hip:0 must never be the 4070 there).
GPU_ARTIFACT ?= hip:0
GPU_LEGS     ?=
GPU_SKIP     ?=
GPU_EXPECT   ?=
.PHONY: check-gpu
check-gpu: check-pins all programs
	$(PYTHON) tools/check_gpu.py --build $(B) --artifact $(GPU_ARTIFACT) --hip '$(HIP)' --expect '$(GPU_EXPECT)' \
	    $(if $(GPU_LEGS),--legs $(GPU_LEGS)) $(if $(GPU_SKIP),--skip $(GPU_SKIP))

# The worked round trip, built in the tree and run against the
# programs here rather than an installed copy.
.PHONY: example
example: $(B)/librebound.a
	$(MAKE) install PREFIX=$(CURDIR)/$(B)/stage
	$(MAKE) -C examples clean
	$(MAKE) -C examples PREFIX=$(CURDIR)/$(B)/stage \
	    REBOUND_INCLUDE=../$(REB) REBOUND_LIB=../$(B)/librebound.a
	examples/roundtrip$(EXE)
	examples/dropin$(EXE)

# A header and a library where a REBOUND build can find them, plus the
# ias15_cft program the library runs. libcft is installed too (it is
# Apache-2.0; see NOTICE) because libcft_rebound has undefined
# references into it and `-lcft` has to resolve.
.PHONY: install uninstall
install: all $(CFTLIB)
	$(INSTALL) -d $(DESTDIR)$(INCLUDEDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 644 include/cft_rebound.h $(DESTDIR)$(INCLUDEDIR)/
	$(INSTALL) -m 644 $(B)/libcft_rebound.a $(DESTDIR)$(LIBDIR)/
	$(INSTALL) -m 644 $(B)/libcft_ias15.a $(DESTDIR)$(LIBDIR)/
	$(INSTALL) -m 644 src/cft_ias15.h src/cft_ias15_state.h src/cft_archive.h \
	    $(DESTDIR)$(INCLUDEDIR)/
	$(INSTALL) -m 755 $(B)/ias15_cft$(EXE) $(DESTDIR)$(BINDIR)/
	$(INSTALL) -m 644 $(CFT)/include/cft.h $(CFT)/include/cft_config.h $(CFT)/include/cft.hpp $(DESTDIR)$(INCLUDEDIR)/
	$(INSTALL) -m 644 $(CFTLIB) $(DESTDIR)$(LIBDIR)/
	@echo
	@echo "installed into $(DESTDIR)$(PREFIX). Two lines in your own build:"
	@echo "    CFLAGS  += -I$(INCLUDEDIR)"
	@echo "    LDLIBS  += -L$(LIBDIR) -lcft_rebound -lcft"
	@echo
	@echo "The Python integrator is NOT part of this: it is built against the"
	@echo "REBOUND wheel a particular interpreter has, so it cannot be a"
	@echo "system-wide artifact the way these are. \`make install-python-lib\`"
	@echo "if you want it in $(LIBDIR) anyway; docs/PYTHON.md says why there"
	@echo "is no pip package."

# Deliberately separate from install, and not a prerequisite of it: it
# needs a REBOUND wheel, which on Windows it must LINK, so folding it in
# would make `make install` fail for every C-only caller without one.
.PHONY: install-python-lib
install-python-lib: $(B)/$(SHLIB)
	$(INSTALL) -d $(DESTDIR)$(LIBDIR)
	$(INSTALL) -m 755 $(B)/$(SHLIB) $(DESTDIR)$(LIBDIR)/
	@test -z "$(CFT_SHLIB_EXTRA)" || \
	    $(INSTALL) -m 755 $(CFT_SHLIB_EXTRA) $(DESTDIR)$(LIBDIR)/
	@echo
	@echo "    cft_rebound.load('$(LIBDIR)/$(SHLIB)')"

uninstall:
	rm -f $(DESTDIR)$(LIBDIR)/$(SHLIB) \
	      $(if $(CFT_SHLIB_EXTRA),$(DESTDIR)$(LIBDIR)/$(notdir $(CFT_SHLIB_EXTRA)))
	rm -f $(DESTDIR)$(INCLUDEDIR)/cft_rebound.h
	rm -f $(DESTDIR)$(INCLUDEDIR)/cft.h $(DESTDIR)$(INCLUDEDIR)/cft_config.h $(DESTDIR)$(INCLUDEDIR)/cft.hpp
	rm -f $(DESTDIR)$(LIBDIR)/libcft_rebound.a $(DESTDIR)$(LIBDIR)/libcft_ias15.a \
	      $(DESTDIR)$(LIBDIR)/libcft.a
	rm -f $(DESTDIR)$(INCLUDEDIR)/cft_ias15.h $(DESTDIR)$(INCLUDEDIR)/cft_ias15_state.h \
	      $(DESTDIR)$(INCLUDEDIR)/cft_archive.h
	rm -f $(DESTDIR)$(BINDIR)/ias15_cft$(EXE)

clean:
	rm -rf $(B)
	# tests/gate_real.c writes its checkpoint beside the binary when it
	# is run by hand rather than through tools/check_checkpoint.py,
	# which uses a temporary directory.
	rm -f gate_real.bin
	$(MAKE) -C examples clean
