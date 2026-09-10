/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft_ias15_state.h - kept as the name the archive includes.
 *
 * The struct itself lives in cft_ias15.h, which the integrator shim
 * owns. This file held a second, byte-identical copy while the shim and
 * the archive were written in parallel; at integration the copies were
 * diffed, found identical, and this one removed. What remains is the
 * one helper the archive added and the include that reaches the real
 * definition.
 */
#ifndef CFT_IAS15_STATE_H
#define CFT_IAS15_STATE_H

#include <stddef.h>
#include <stdint.h>
#include "cft.h"
#include "cft_ias15.h"

#define CFT_IAS15_STATE_ABI 1

static inline size_t cft_ias15_state_width(const struct cft_ias15_state *s){
    return cft_format_size((cft_format)s->format);
}

#endif /* CFT_IAS15_STATE_H */
