# SPDX-License-Identifier: GPL-3.0-or-later
#
# cft-rebound: IAS15 with its arithmetic through libcft.
#
#   make third-party      clone/verify the pinned upstreams (needs git)
#   make libcft           build libcft from the pinned cft-fp256 clone
#   make librebound       build REBOUND's C library from the pinned clone
#   make                  the reference program and the libcft port
#   make constants        re-derive and re-check the IAS15 constants
#   make check            the gates: constants, and the port at binary64
#                         against REBOUND's own IAS15, bit for bit
#   make archive          the Simulationarchive gate programs
#   make check-archive    just those gates (seconds)
#   make example          build the worked round trip and run it
#   make install          a header and a library where a REBOUND build
#                         can find them; PREFIX and DESTDIR as usual
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
  endif
endif

all: $(B)/ias15_ref$(EXE) $(B)/ias15_cft$(EXE) $(B)/check_dropin$(EXE) $(B)/libcft_rebound.a

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

$(B)/cft_rebound_run.o: src/cft_rebound_run.c include/cft_rebound.h src/hexfloat.h $(B)/bindir.stamp
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Iinclude -Isrc -I$(REB) -I$(CFT)/include \
	    -DCFT_REBOUND_IAS15_DEFAULT='"$(BINDIR)/ias15_cft$(EXE)"' -o $@ src/cft_rebound_run.c

$(B)/libcft_rebound.a: $(B)/cft_rebound_run.o
	ar rcs $@ $^

$(B)/ias15_cft$(EXE): src/ias15_cft.c src/ias15_constants.h src/hexfloat.h $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) -I$(CFT)/include -o $@ src/ias15_cft.c $(CFTLIB) $(LIBS)

# ---------------------------------------------------------------------
# The drop-in: the same engine compiled as a library (no main) plus the
# struct reb_integrator shim, and the gate that runs it against
# REBOUND's own ias15 inside one program.
#
# src/cft_ias15_fields.c is PARCEL B's file: it defines
# cft_ias15_field_descriptor_list and today holds only the terminator.
DROPIN_OBJ := $(B)/ias15_cft_lib.o $(B)/reb_integrator_cft.o $(B)/cft_ias15_fields.o

$(B)/ias15_cft_lib.o: src/ias15_cft.c src/ias15_constants.h src/ias15_engine.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) -DIAS15_CFT_LIBRARY -Isrc -I$(CFT)/include -o $@ src/ias15_cft.c

$(B)/reb_integrator_cft.o: src/reb_integrator_cft.c src/cft_ias15.h src/ias15_engine.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/reb_integrator_cft.c

$(B)/cft_ias15_fields.o: src/cft_ias15_fields.c src/cft_ias15.h
	@mkdir -p $(B)
	$(CC) -c $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) -o $@ src/cft_ias15_fields.c

$(B)/libcftrebound.a: $(DROPIN_OBJ)
	ar rcs $@ $^

.PHONY: dropin
dropin: $(B)/libcftrebound.a $(B)/check_dropin$(EXE)

$(B)/check_dropin$(EXE): tools/check_dropin.c $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) $(REB_USEFLAGS) -Isrc -I$(CFT)/include -I$(REB) 	      -o $@ tools/check_dropin.c $(DROPIN_OBJ) $(B)/librebound.a $(CFTLIB) $(LIBS)

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
# completeness - their step touches all 48 blobs, and 20 steps of a real
# IAS15 need not.
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

# The worked round trip, built in the tree and run against the
# programs here rather than an installed copy.
.PHONY: example
example: $(B)/librebound.a
	$(MAKE) install PREFIX=$(CURDIR)/$(B)/stage
	$(MAKE) -C examples clean
	$(MAKE) -C examples PREFIX=$(CURDIR)/$(B)/stage \
	    REBOUND_INCLUDE=../$(REB) REBOUND_LIB=../$(B)/librebound.a
	examples/roundtrip$(EXE)

# A header and a library where a REBOUND build can find them, plus the
# ias15_cft program the library runs. libcft is installed too (it is
# Apache-2.0; see NOTICE) because libcft_rebound has undefined
# references into it and `-lcft` has to resolve.
.PHONY: install uninstall
install: all $(CFTLIB)
	$(INSTALL) -d $(DESTDIR)$(INCLUDEDIR) $(DESTDIR)$(LIBDIR) $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 644 include/cft_rebound.h $(DESTDIR)$(INCLUDEDIR)/
	$(INSTALL) -m 644 $(B)/libcft_rebound.a $(DESTDIR)$(LIBDIR)/
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
	rm -f $(DESTDIR)$(LIBDIR)/libcft_rebound.a $(DESTDIR)$(LIBDIR)/libcft.a
	rm -f $(DESTDIR)$(BINDIR)/ias15_cft$(EXE)

clean:
	rm -rf $(B)
	$(MAKE) -C examples clean
