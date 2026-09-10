/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * A stand-in for parcel A's REBOUND integrator shim, used ONLY by
 * parcel B's archive gates. Its step is NOT IAS15: it is a deterministic
 * mixing over the wide state, at a real libcft format, arranged so that
 * every one of 48 archived blobs feeds the next step. That is what the
 * archive gate needs - if one blob were missing from the descriptor
 * list, or written at the wrong width, a restart would not continue
 * bit-identically - and it is deliberately independent of whether the
 * real IAS15 step is finished.
 *
 * 48, not CFT_N_BLOBS (50): shim_step() never reads or writes the live
 * x and v, the two blobs added on 2026-09-10. Those two are still
 * covered, but by a weaker mechanism - cft_shim_setup() seeds every
 * blob from index 2 up and cft_shim_snapshot() compares all
 * CFT_N_BLOBS of them - so a descriptor dropped for x or v shows up in
 * the snapshot comparison rather than in the step's forward
 * dependency. tests/gate_real.c is the gate that exercises x and v as
 * live state.
 *
 * This header used to say "when parcel A lands, the gates should be
 * re-pointed at the real integrator". Parcel A landed and they were
 * deliberately NOT re-pointed: tests/gate_real.c is a separate gate
 * that drives the real shim through an archive, and these four keep
 * driving the stub because its step touches 48 of the 50 blobs (see
 * above) and 20 steps of a real IAS15 need not. Leaving
 * that sentence standing for a whole round is what let the seam go
 * untested - docs/VALIDATION.md entry 25.
 *
 * The blob walker and the blob count are NOT duplicated here any more:
 * cft_archive_state_blob() and CFT_N_BLOBS are the definitions, and a
 * private copy of the walker in cft_shim_stub.c returned NULL for the
 * two blobs added on 2026-09-10 and crashed three gates.
 */
#ifndef CFT_SHIM_STUB_H
#define CFT_SHIM_STUB_H

#include "rebound.h"
#include "cft_archive.h"

/* Register under CFT_IAS15_INTEGRATOR_NAME. Safe to call twice. */
void cft_shim_register(void);

/* Set the format, allocate and seed the wide state for r->N particles.
 * Must be called after the particles are added. */
int  cft_shim_setup(struct reb_simulation *r, int format);

/* Bytes of the whole wide state plus the particles, for bit comparison.
 * Caller frees. */
unsigned char *cft_shim_snapshot(struct reb_simulation *r, size_t *len);

/* How many of the CFT_N_BLOBS blobs hold at least one non-zero byte. A
 * gate that
 * compares all-zero states proves nothing. */
int  cft_shim_blobs_nonzero(struct reb_simulation *r);

#endif
