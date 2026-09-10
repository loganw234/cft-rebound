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

all: $(B)/ias15_ref$(EXE) $(B)/ias15_cft$(EXE)

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

$(B)/ias15_cft$(EXE): src/ias15_cft.c src/ias15_constants.h src/hexfloat.h $(CFTLIB)
	$(CC) $(CSTD) $(CFLAGS) $(WARN) -I$(CFT)/include -o $@ src/ias15_cft.c $(CFTLIB) $(LIBS)

constants:
	$(PYTHON) tools/gen_constants.py

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

check: all programs
	$(PYTHON) tools/gen_constants.py --no-write
	$(PYTHON) tools/check_equivalence.py --build $(B)
	$(PYTHON) tools/check_program_engine.py --build $(B)
	$(PYTHON) tools/check_records.py --build $(B)
	$(PYTHON) tools/check_ensemble.py --build $(B)

# the same gates at binary64 only, in a few minutes
.PHONY: check-quick
check-quick: all programs
	$(PYTHON) tools/gen_constants.py --no-write
	$(PYTHON) tools/check_equivalence.py --build $(B) --quick
	$(PYTHON) tools/check_program_engine.py --build $(B) --formats fp64
	$(PYTHON) tools/check_records.py --build $(B) --quick
	$(PYTHON) tools/check_ensemble.py --build $(B) --quick

clean:
	rm -rf $(B)
