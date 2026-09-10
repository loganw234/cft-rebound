# SPDX-License-Identifier: GPL-3.0-or-later
#
# Building a REBOUND program against an INSTALLED cft-rebound.
#
# This is a normal REBOUND program's makefile with the cft-rebound lines
# added, and they are marked below. There are two ways in and they are
# different libraries, so copy the pair that matches the one you want:
# CFLAGS plus CFTLIBS for the subprocess API (roundtrip.c), CFLAGS plus
# DROPINLIBS for the registered integrator (dropin.c).
#
#   make PREFIX=/usr/local
#   ./roundtrip
#   ./dropin
#
# If REBOUND itself is not installed system-wide, say where it is:
#
#   make PREFIX=/usr/local \
#        REBOUND_INCLUDE=../third_party/rebound/src \
#        REBOUND_LIB=../build/librebound.a
#
# (The repository root's `make example` does exactly that, and runs it.)

CC      ?= gcc
PREFIX  ?= /usr/local
CFLAGS  ?= -O2
CSTD    := -std=c99

# Your own REBOUND. Leave REBOUND_INCLUDE empty and REBOUND_LIB at
# -lrebound if it is installed where the compiler already looks.
REBOUND_INCLUDE ?=
REBOUND_LIB     ?= -lrebound

# ---- the two lines a cft-rebound user adds ------------------------
CFLAGS  += -I$(PREFIX)/include
CFTLIBS += -L$(PREFIX)/lib -lcft_rebound -lcft
# -------------------------------------------------------------------
# The other way in: the registered integrator, in this process. Same
# shape, a different library.
DROPINLIBS += -L$(PREFIX)/lib -lcft_ias15 -lcft

ifeq ($(OS),Windows_NT)
  EXE     := .exe
  # librebound is linked statically here, so rebound.h must not mark
  # its API dllimport. Your own build already handles this its own way.
  CFLAGS  += -DBUILDINGLIBREBOUND
  SYSLIBS := -lws2_32
else
  EXE     :=
  SYSLIBS := -lm
endif

INC := $(if $(REBOUND_INCLUDE),-I$(REBOUND_INCLUDE),)

.PHONY: all clean
all: roundtrip$(EXE) dropin$(EXE)

roundtrip$(EXE): roundtrip.c
	$(CC) $(CSTD) $(CFLAGS) $(INC) -o $@ $< $(CFTLIBS) $(REBOUND_LIB) $(SYSLIBS)

dropin$(EXE): dropin.c
	$(CC) $(CSTD) $(CFLAGS) $(INC) -o $@ $< $(DROPINLIBS) $(REBOUND_LIB) $(SYSLIBS)

clean:
	rm -f roundtrip roundtrip.exe dropin dropin.exe dropin_checkpoint.bin
