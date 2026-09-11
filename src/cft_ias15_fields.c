/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The one definition of the Simulationarchive field descriptors for
 * struct cft_ias15_state, and the list the integrator registers with.
 *
 * REBOUND writes a registered integrator's state by walking
 * r->integrator.callbacks.field_descriptor_list with
 * r->integrator.state as the base and "integrator.<name>." as the
 * prefix, so the archived names come out as
 *
 *     integrator.ias15_cft.cft_x0
 *
 * On READ it resolves each name against the simulation's list, then the
 * built-in integrators, then the registered custom ones
 * (binarydata.c, reb_binarydata_field_descriptor_for_name). All three of
 * those are whatever registration installed, so this list has to be a
 * real one or nothing is read back - which is what tests/gate_real.c
 * caught when this file held only a terminator.
 *
 * WIDTH. element_size is a compile-time constant in a REBOUND
 * descriptor and the wide element is 8, 16 or 32 bytes according to
 * state->format, so there are three lists rather than one.
 * cft_archive_bind() installs the one matching the state before a write,
 * which is what makes size_data right on disk.
 *
 * For a READ the registered list is cft_fd_fp64 whatever the archive's
 * format, and that is correct rather than merely tolerable: REBOUND
 * allocates field.size_data bytes from the file and reads them, using
 * element_size only for *pointer_N = size_data/element_size.
 * cft_archive_finish_load() overwrites n_elem from the file for exactly
 * that reason (its NOTE 3), and 8 divides 16 and 32, so the
 * "Inconsistent size_data" warning cannot fire for a wide archive read
 * through the narrow list.
 */
#include "cft_ias15_fields.h"

const struct reb_binarydata_field_descriptor cft_fd_fp64[]  = CFT_FD_LIST(8);
const struct reb_binarydata_field_descriptor cft_fd_fp128[] = CFT_FD_LIST(16);
const struct reb_binarydata_field_descriptor cft_fd_fp256[] = CFT_FD_LIST(32);

/* What struct reb_integrator carries. See WIDTH above for why the
 * narrow list is the right one to register. */
const struct reb_binarydata_field_descriptor cft_ias15_field_descriptor_list[] = CFT_FD_LIST(8);

/* cft_ias15_state.provenance in words. Here rather than in
 * src/cft_archive.c because a caller holding a state may want it and
 * not every target that holds a state links the archive module:
 * DROPIN_OBJ in the Makefile is the engine, the shim, this file and
 * the refusal table, and build/check_dropin links exactly that. */
const char *cft_ias15_provenance_str(int provenance){
    switch (provenance){
        case CFT_PROV_NONE:     return "created here";
        case CFT_PROV_EXACT:    return "restored exactly";
        case CFT_PROV_PROMOTED: return "promoted from binary64";
    }
    return "?";
}
