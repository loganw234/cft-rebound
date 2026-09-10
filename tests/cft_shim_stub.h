/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * A stand-in for parcel A's REBOUND integrator shim, used ONLY by
 * parcel B's archive gates. Its step is NOT IAS15: it is a deterministic
 * mixing over the wide state, at a real libcft format, arranged so that
 * every one of the 48 archived blobs feeds the next step. That is what
 * the archive gate needs - if one blob were missing from the descriptor
 * list, or written at the wrong width, a restart would not continue
 * bit-identically - and it is deliberately independent of whether the
 * real IAS15 step is finished.
 *
 * When parcel A lands, the gates should be re-pointed at the real
 * integrator; the archive module itself does not change.
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

/* How many of the 48 blobs hold at least one non-zero byte. A gate that
 * compares all-zero states proves nothing. */
int  cft_shim_blobs_nonzero(struct reb_simulation *r);

#endif
