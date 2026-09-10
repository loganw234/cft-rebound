/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * A stand-in for REBOUND's integrator_whfast512.c on MinGW-w64.
 *
 * REBOUND's WHFast512 file defines its own `__m512d` and `__mmask8`
 * types so that it need not be compiled with AVX-512 flags. On
 * MinGW-w64 (gcc 16, msys2) rebound_internal.h pulls in <winsock2.h>
 * -> <windows.h> -> <winnt.h> -> <x86intrin.h>, which already defines
 * `__m512d`, and the two definitions conflict:
 *
 *   integrator_whfast512.c:60:28: error: conflicting types for '__m512d'
 *
 * REBOUND's own Windows build assumes MSVC (cl.exe), where windows.h
 * does not drag the intrinsics in. This repository builds REBOUND's
 * library without that one file and links this stub instead, so that
 * the built-in integrator table in simulation.c still resolves. IAS15
 * does not touch WHFast512; selecting "whfast512" from this build is
 * refused, not silently mis-run. */
#include "rebound.h"
#include <stdio.h>

static void* whfast512_stub_create(void){
    fprintf(stderr, "cft-rebound: whfast512 is not built in this tree (see ref/whfast512_stub.c)\n");
    return NULL;
}
static void whfast512_stub_free(void* p){ (void)p; }
static void whfast512_stub_step(struct reb_simulation* r, void* p){
    (void)p;
    reb_simulation_error(r, "whfast512 is not built in this tree (see ref/whfast512_stub.c)");
}

const struct reb_integrator reb_integrator_whfast512 = {
    .documentation = "Not built in cft-rebound's tree; see ref/whfast512_stub.c.",
    .step = whfast512_stub_step,
    .synchronize = NULL,
    .create = whfast512_stub_create,
    .free = whfast512_stub_free,
    .did_add_particle = NULL,
    .will_remove_particle = NULL,
    .field_descriptor_list = NULL,
};
