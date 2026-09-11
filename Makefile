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
#   make check            every gate at every format (about half an hour)
#   make check-quick      the same at binary64 only (a few minutes)
#   make archive          the Simulationarchive gate programs
#   make check-archive    just those gates (seconds)
#   make dropin           the drop-in library and its gate
#   make python-lib       the shared library for Python; check-python
#                         runs the equivalence gate through it
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

# The shared library for Python. POSIX only: a Windows DLL may not
# carry undefined symbols, so it would have to link against the
# wheel's librebound rather than leave REBOUND's symbols to the
# loader, and that is a different rule - see docs/PYTHON.md.
ifeq ($(shell uname -s 2>/dev/null),Darwin)
  SHLIB     := libcft_ias15.dylib
  SHFLAGS   := -dynamiclib -Wl,-install_name,@rpath/libcft_ias15.dylib
  CFT_SHLIB := $(CFT)/libcft.dylib
else
  SHLIB     := libcft_ias15.so
  SHFLAGS   := -shared
  CFT_SHLIB := $(CFT)/libcft.so
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

$(B)/ias15_cft$(EXE): src/ias15_cft.c src/ias15_constants.h src/hexfloat.h src/ias15_limits.h $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) -I$(CFT)/include -o $@ src/ias15_cft.c $(CFTLIB) $(LIBS)

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

$(B)/ias15_cft_lib.o: src/ias15_cft.c src/ias15_constants.h src/ias15_engine.h src/ias15_limits.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) -DIAS15_CFT_LIBRARY -Isrc -I$(CFT)/include -o $@ src/ias15_cft.c

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
# behind their back.
PIC_OBJ := $(B)/pic/ias15_cft_lib.o $(B)/pic/reb_integrator_cft.o \
           $(B)/pic/cft_ias15_fields.o $(B)/pic/cft_archive.o \
           $(B)/pic/cft_ias15_shared.o $(B)/pic/cft_supported.o
PICFLAGS := $(CSTD) $(CFLAGS) $(WARN) -fPIC -Isrc -I$(CFT)/include -I$(REB)

$(B)/pic/ias15_cft_lib.o: src/ias15_cft.c src/ias15_constants.h src/ias15_engine.h src/ias15_limits.h
	@mkdir -p $(B)/pic
	$(CC) -c $(PICFLAGS) -DIAS15_CFT_LIBRARY -o $@ src/ias15_cft.c

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

# REBOUND's symbols are deliberately NOT linked: the caller's process
# already has librebound loaded and python/cft_rebound.py promotes it
# to RTLD_GLOBAL before opening this. Linking a second copy would give
# a second integrator list and a name that can never be selected.
# libcft.a will NOT go into a shared object: backend_xrt.o is C++ and
# carries relocations a shared object cannot use (R_X86_64_PC32 against
# a GLIBCXX symbol). cft-fp256 builds a proper shared libcft from its
# own PIC objects, so use that. Two rpaths: $$ORIGIN for an installed
# layout where the two sit together, and the tree for running it in
# place without an install.
$(CFT_SHLIB):
	$(MAKE) -C $(CFT) CC=$(CC) PYTHON=$(PYTHON) $(CFT_MAKEVARS) $(notdir $(CFT_SHLIB))

$(B)/$(SHLIB): $(PIC_OBJ) $(CFT_SHLIB)
	$(CC) $(SHFLAGS) -o $@ $(PIC_OBJ) -L$(CFT) -lcft \
	    -Wl,-rpath,'$$ORIGIN' -Wl,-rpath,$(abspath $(CFT)) $(LIBS)

.PHONY: python-lib check-python
python-lib: $(B)/$(SHLIB) check-rebound-match
	@echo
	@echo "built $(B)/$(SHLIB). From Python:"
	@echo "    import cft_rebound"
	@echo "    cft_rebound.load('$(CURDIR)/$(B)/$(SHLIB)')"

# The library is compiled against the pinned REBOUND headers and is
# loaded into a process running the WHEEL's REBOUND, so the two must be
# the same source or struct reb_simulation is laid out differently on
# each side of the call and nothing says so. They are identical for
# rebound 5.1.1 and the pinned bdfda4bd; this refuses to build a
# library that would be quietly wrong if that ever stops being true.
.PHONY: check-rebound-match
check-rebound-match:
	@w=`$(PYTHON) -c 'import rebound,os;print(os.path.dirname(rebound.__libpath__))' 2>/dev/null`; \
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
# from Python, on the same problem and the same fixed step.
check-python: $(B)/$(SHLIB)
	$(PYTHON) python/example_equivalence.py --library $(CURDIR)/$(B)/$(SHLIB)

.PHONY: dropin
dropin: $(B)/libcft_ias15.a $(B)/check_dropin$(EXE)

# The gate is a file per topic now: main() and the refusals in
# check_dropin.c, the scaffolding in dropin_common.c, the cases in
# cases_*.c. A parcel adds its file here and nowhere else -
# tools/dropin_cases.h says how.
CASES_SRC := tools/check_dropin.c tools/dropin_common.c \
             tools/cases_core.c tools/cases_wide.c \
             tools/cases_forces.c tools/cases_modes.c

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
                 $(B)/gate_real$(EXE)

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
$(ASM):
	$(MAKE) -C $(CFT) CC=$(CC) PYTHON=$(PYTHON) $(CFT_MAKEVARS) cft-asm$(EXE)

.PHONY: programs
programs: $(ASM)
	$(PYTHON) tools/gen_programs.py
	@mkdir -p programs/out
	@for f in programs/*.cfta; do b=$$(basename $$f .cfta); $(ASM) $$f -o programs/out/$$b.cftp || exit 1; done
	@echo "assembled $$(ls programs/out/*.cftp | wc -l) programs"

check: all programs $(ARCHIVE_GATES)
	$(PYTHON) tools/gen_constants.py --no-write
	$(B)/check_dropin$(EXE)
	$(B)/check_dropin$(EXE) --wide
	$(PYTHON) tools/check_equivalence.py --build $(B)
	$(PYTHON) tools/check_program_engine.py --build $(B)
	$(PYTHON) tools/check_records.py --build $(B)
	$(PYTHON) tools/check_ensemble.py --build $(B)
	$(PYTHON) tools/check_archive.py --build $(B)
	$(PYTHON) tools/check_bodycount.py --build $(B)
	$(B)/gate_real$(EXE) --fp64
	$(B)/gate_real$(EXE) --fp128
	$(B)/gate_real$(EXE) --fp256
	$(PYTHON) tools/check_checkpoint.py --build $(B)

# the same gates at binary64 only, in a few minutes
.PHONY: check-quick
check-quick: all programs $(ARCHIVE_GATES)
	$(PYTHON) tools/gen_constants.py --no-write
	$(B)/check_dropin$(EXE)
	$(B)/check_dropin$(EXE) --wide
	$(PYTHON) tools/check_equivalence.py --build $(B) --quick
	$(PYTHON) tools/check_program_engine.py --build $(B) --formats fp64
	$(PYTHON) tools/check_records.py --build $(B) --quick
	$(PYTHON) tools/check_ensemble.py --build $(B) --quick
	$(PYTHON) tools/check_archive.py --build $(B)
	$(PYTHON) tools/check_bodycount.py --build $(B) --quick
	$(B)/gate_real$(EXE) --fp64
	$(B)/gate_real$(EXE) --fp128
	$(B)/gate_real$(EXE) --fp256
	$(PYTHON) tools/check_checkpoint.py --build $(B) --quick

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

uninstall:
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
