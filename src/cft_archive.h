/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft_archive.h - Simulationarchive support for the cft IAS15 integrator.
 *
 * REBOUND identifies archive fields by NAME, not by index, and an
 * integrator's own state is written from the field_descriptor_list its
 * `struct reb_integrator` carries, under the prefix
 * "integrator.<integrator name>.". Every name added here begins with
 * `cft_`, so it cannot collide with any present or future REBOUND field.
 *
 * Two things a reader of ROADMAP.md should know before reading further,
 * because the source says something narrower than the roadmap does:
 *
 *  1. binarydata.c's REB_FIELD_NOT_FOUND case warns, seeks past the
 *     field, and then `goto finish_fields` - it STOPS reading the
 *     snapshot; it does not continue to the next field. The cft_ fields
 *     are therefore only safe for a stock reader because REBOUND writes
 *     every simulation field (particles included) BEFORE any integrator
 *     field. What a stock reader loses is whatever follows: the
 *     "functionpointers" flag.
 *
 *  2. Integrator fields are additionally gated on the integrator NAME:
 *     an archive naming an integrator the reader has not registered
 *     produces "Integrator not found." from reb_simulation_set_integrator
 *     and then REB_BINARYDATA_WARNING_CORRUPTFILE at the first
 *     integrator-prefixed field - not the FIELD_UNKNOWN warning. The
 *     binary64 state still arrives intact. docs/VALIDATION.md records
 *     the verbatim output.
 */
#ifndef CFT_ARCHIVE_H
#define CFT_ARCHIVE_H

#include <stdint.h>
#include "rebound.h"
#include "cft.h"
#include "cft_ias15.h"        /* CFT_IAS15_INTEGRATOR_NAME, the state */
#include "cft_ias15_state.h"

/* The name this project registers its integrator under. The archive
 * writes its fields as "integrator.<CFT_IAS15_INTEGRATOR_NAME>.cft_*".
 * The name itself is the integrator's, and is defined in cft_ias15.h
 * beside the call that registers it. */

/* What a load did. */
enum cft_archive_status {
    CFT_ARCHIVE_EXACT           =  0, /* cft_ fields present and consistent: the wide state was restored bit for bit */
    CFT_ARCHIVE_PROMOTED        =  1, /* no cft_ fields: REBOUND's binary64 state was promoted */
    CFT_ARCHIVE_FORMAT_MISMATCH = -1, /* the archive's cft_format is not the run's format: refused */
    CFT_ARCHIVE_COUNT_MISMATCH  = -2, /* the element count and the blob lengths disagree: refused */
    CFT_ARCHIVE_NO_FILE         = -3,
    CFT_ARCHIVE_BAD             = -4  /* unreadable or corrupt */
};

const char *cft_archive_status_str(enum cft_archive_status s);

/* Check the shape of the three descriptor lists (names, widths, order).
 * 0 if sound, otherwise the number of complaints, each on stderr. */
int cft_archive_selftest(void);

/* What one snapshot of an archive says about itself, read from the file
 * itself rather than inferred from a loaded simulation. For a snapshot
 * other than 0 the answer is snapshot 0 overlaid with that snapshot's
 * diff, which is what REBOUND's own loader reconstructs. */
struct cft_archive_info {
    char     integrator[256];  /* the archive's integrator.name, "" if absent */
    int      has_cft;          /* at least one integrator.<name>.cft_* field */
    int      n_cft_fields;
    int      format;           /* cft_format, -1 if absent */
    uint64_t n_elem;           /* cft_n_elem, 0 if absent */
    uint64_t E;                /* cft_E, 0 if absent */
    uint64_t blob_bytes;       /* size_data of cft_x0, 0 if absent */
    int      n_blobs_seen;     /* how many of the CFT_N_BLOBS wide blobs are present */
    int      blob_lengths_agree;
    char     abi[17];          /* cft_abi, NUL terminated */
    uint64_t constants_digest;
    int64_t  n_snapshots;

    /* The alias family, counted apart from the blobs above because its
     * length is a different member (see CFT_FD_A in
     * src/cft_ias15_fields.h). Every one of them is 0 for an archive
     * written by a run that never removed a particle - and so for every
     * archive this project wrote before 2026-09-11, which is why the
     * family is optional on read rather than required. */
    uint64_t hiwater_n_elem;   /* cft_hiwater_n_elem: 3*N_allocated */
    int      has_hiwater;      /* the field is NAMED, whatever its value.
                                * A zero value and an absent field are
                                * different facts: an appended snapshot is
                                * a diff, so a mark that has risen back to
                                * the live count writes a 0 there while the
                                * base snapshot's blobs stay in the file. */
    uint64_t alias_blob_bytes; /* size_data of cft_alias_g0 */
    int      n_alias_blobs_seen;
    int      alias_lengths_agree;

    /* What the WRITING run's state said about itself. provenance is a
     * cft_ias15_provenance; -1 when the archive names no such field,
     * which every archive written before 2026-09-11 does. */
    int      provenance;
    uint64_t iterations_max_exceeded;
};

/* The high-water body count an archive implies: hiwater_n_elem/3 when
 * it carries the alias family, and otherwise the live count, which is
 * what it means for the family to be absent. One place, because the
 * shim and the gates both need it and REBOUND recovers its own
 * N_allocated by exactly this division. */
uint64_t cft_archive_hiwater_N(const struct cft_archive_info *info);

/* Read one snapshot's field table. Returns 0 on success. */
int cft_archive_probe(const char *filename, int64_t snapshot,
                      struct cft_archive_info *out);

/* The field_descriptor_list for a given libcft format. The list's
 * REB_POINTER entries carry element_size = cft_format_size(format), so
 * the list in force must match the state's format whenever a simulation
 * is WRITTEN. Returns NULL for a format libcft does not have. */
const struct reb_binarydata_field_descriptor *cft_archive_descriptor_list(int format);

/* Point r's integrator at the descriptor list for the state's current
 * format. Call after changing state->format and after a load. Returns 0
 * on success. */
int cft_archive_bind(struct reb_simulation *r);

/* Allocate / free the wide arrays for n_elem elements at the state's
 * format. Zero-filled (+0 in every format). The archive owns these two
 * only so that the promotion path can allocate; the integrator shim may
 * use them or its own. */
/* The i'th wide blob of a state, by the order of CFT_FD_BLOBS, or NULL
 * past the end. The one walker: anything that iterates the blobs uses
 * this, so a list that grows cannot leave a second copy behind. */
unsigned char **cft_archive_state_blob(struct cft_ias15_state *s, int i);

/* The i'th blob of the ALIAS family, by the order of CFT_FD_ALIAS, or
 * NULL past the end: alias[0][0..6] through alias[5][0..6], then
 * alias_csx and alias_csv. These are sized by hiwater_n_elem rather
 * than by n_elem, which is the whole reason they are a second family. */
unsigned char **cft_archive_state_alias(struct cft_ias15_state *s, int i);

int  cft_archive_state_alloc(struct cft_ias15_state *s, size_t n_elem);
void cft_archive_state_free(struct cft_ias15_state *s);

/* Promote n binary64 values into fmt, exactly (binary64 is a subset of
 * every wider format). dev may be NULL for the software backend. */
int cft_archive_promote_doubles(cft_device *dev, int fmt, const double *src,
                                unsigned char *dst, size_t n);

/* Round n wide values back to binary64, once, to nearest-even. */
int cft_archive_demote_doubles(cft_device *dev, int fmt, const unsigned char *src,
                               double *dst, size_t n);

/* Load a simulation from an archive.
 *
 *   want_format  one of CFT_FP64/FP128/FP256, or -1 to adopt whatever
 *                the archive says (and CFT_FP64 when it says nothing).
 *   dev          a device for the promotion path; NULL opens and closes
 *                a software one.
 *
 * If the archive carries cft_ fields and the format and element count
 * agree, the wide state is restored exactly (CFT_ARCHIVE_EXACT).
 * If it carries none, REBOUND's binary64 state - the particles, and
 * REBOUND's own IAS15 arrays when the archive was written by "ias15" -
 * is promoted into the wide arrays (CFT_ARCHIVE_PROMOTED).
 * A disagreement is refused: a message on stderr and NULL, never a
 * silently reinterpreted state.
 *
 * Returns NULL on every negative status. *status is always written. */
struct reb_simulation *cft_archive_load(const char *filename, int64_t snapshot,
                                        int want_format, cft_device *dev,
                                        enum cft_archive_status *status);

/* The same decision applied to a simulation the caller has already
 * loaded (with reb_simulation_create_from_file, say). `info` must come
 * from cft_archive_probe on the same file and snapshot. */
enum cft_archive_status cft_archive_finish_load(struct reb_simulation *r,
                                                const struct cft_archive_info *info,
                                                int want_format, cft_device *dev);

#endif /* CFT_ARCHIVE_H */
