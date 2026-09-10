/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * PLACEHOLDER - THIS FILE IS PARCEL B'S.
 *
 * The Simulationarchive field descriptors for struct cft_ias15_state.
 * REBOUND writes a registered integrator's state by walking this list
 * (binarydata.c, output_fields_from_list) with `r->integrator.state` as
 * the base address and "integrator.<registered name>." as the prefix, so
 * the archived names come out as
 *
 *     integrator.ias15_cft.cft_x0
 *
 * and a reader that does not know them warns and seeks past
 * (REB_BINARYDATA_WARNING_FIELD_UNKNOWN). The list is terminated by an
 * entry whose name is empty; what is here is that terminator and nothing
 * else, which REBOUND reads as "this integrator adds no fields" - a cft
 * archive written today is therefore a plain, complete binary64 REBOUND
 * archive, and parcel B is what adds the wide state to it.
 *
 * The shape parcel B fills in, from the ROADMAP:
 *
 *   { "", REB_DOUBLE,  "cft_epsilon", offsetof(struct cft_ias15_state, epsilon), 0, 0, 0 },
 *   { "", REB_POINTER, "cft_x0",      offsetof(struct cft_ias15_state, x0),
 *                                     offsetof(struct cft_ias15_state, n_elem),
 *                                     W, 0 },        // W = cft_format_size(format)
 *   ...
 *
 * with W the wide element width. Note that element_size is a compile-time
 * constant in every REBOUND descriptor and W is not: it is 8, 16 or 32
 * according to state->format. That is parcel B's first real decision and
 * it is written up in the handover note at the end of docs/VALIDATION.md.
 */
#include "rebound.h"
#include "cft_ias15.h"

const struct reb_binarydata_field_descriptor cft_ias15_field_descriptor_list[] = {
    { 0 },   /* null terminated list */
};
