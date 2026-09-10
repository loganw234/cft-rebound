/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * What makes the shared library usable from Python: it registers
 * itself when it is loaded.
 *
 * python/cft_rebound.py's load() says so - "registration happens in the
 * library's own constructor, so the name is usable as soon as this
 * returns" - because there is nowhere else for it to happen. REBOUND's
 * Python layer selects an integrator by assigning a string to
 * sim.integrator, which reaches C as a strcmp against the registered
 * list; there is no call a Python caller could make first without a
 * binding layer of its own, which is exactly what this avoids needing.
 *
 * The order that makes it work is the loader's: it promotes the
 * already-loaded librebound to RTLD_GLOBAL before opening this library,
 * so reb_integrator_register resolves against the same REBOUND the
 * caller's simulations live in. Opened any other way - RTLD_LOCAL, or
 * before `import rebound` - the constructor runs against an unresolved
 * symbol and dlopen fails, which is the right failure: a second
 * REBOUND would be a second integrator list and a silently unusable
 * name.
 *
 * This file is compiled ONLY into the shared library. The static
 * libcft_ias15.a deliberately does not carry it: a C caller links the
 * archive and calls cft_ias15_register() itself, at a point of its own
 * choosing, which is what examples/dropin.c does.
 */
#include "cft_ias15.h"

#if defined(__GNUC__) || defined(__clang__)

/* Registering twice is not merely wasteful here: reb_integrator_register
 * increments its scan index before testing, so a second registration in
 * one process reaches the list's {0} terminator and hands a NULL name to
 * strcmp - see the note in src/cft_ias15.h. cft_ias15_register() returns
 * early if the same name is already registered, so loading this library
 * twice is safe; loading it AND a second library that registers the same
 * name is not, and nothing here can prevent that. */
__attribute__((constructor))
static void cft_ias15_autoregister(void){
    cft_ias15_register(CFT_IAS15_INTEGRATOR_NAME);
}

#else
# error "the shared library needs a constructor; add one for this compiler"
#endif

/* Exported so a caller that loaded the library some other way can do it
 * by hand, and so `nm -D` shows something recognisable. */
void cft_ias15_shared_register(void){
    cft_ias15_register(CFT_IAS15_INTEGRATOR_NAME);
}
