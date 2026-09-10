/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft_archive.c - Simulationarchive support for the cft IAS15 integrator.
 *
 * The field names mirror REBOUND's own IAS15 descriptor list in
 * integrator_ias15.c so the mapping can be checked line by line:
 *
 *   REBOUND                                  here
 *   epsilon         REB_DOUBLE               cft_epsilon      REB_DOUBLE
 *   min_dt          REB_DOUBLE               cft_min_dt       REB_DOUBLE
 *   adaptive_mode   REB_UINT                 cft_adaptive_mode REB_INT   (int in the shared struct)
 *   at              REB_POINTER  8           - see NOTE 1
 *   x0              REB_POINTER  8           cft_x0           REB_POINTER  W
 *   v0              REB_POINTER  8           cft_v0           REB_POINTER  W
 *   a0              REB_POINTER  8           cft_a0           REB_POINTER  W
 *   csx             REB_POINTER  8           cft_csx          REB_POINTER  W
 *   csv             REB_POINTER  8           cft_csv          REB_POINTER  W
 *   csa0            REB_POINTER  8           cft_csa0         REB_POINTER  W
 *   g               REB_POINTER  7*8         cft_g0 .. cft_g6 REB_POINTER  W   NOTE 2
 *   b               REB_POINTER  7*8         cft_b0 .. cft_b6
 *   csb             REB_POINTER  7*8         cft_csb0 .. cft_csb6
 *   e               REB_POINTER  7*8         cft_e0 .. cft_e6
 *   br              REB_POINTER  7*8         cft_br0 .. cft_br6
 *   er              REB_POINTER  7*8         cft_er0 .. cft_er6
 *
 * NOTE 1. REBOUND archives `at`, the acceleration at the current
 * substep. ROADMAP.md's mirror list names it too, but the struct it
 * defines has no `at` member, so there is nothing here to point a
 * descriptor at. `at` is written before it is read in every substep of
 * REBOUND's own step_try, so a restart does not need it; if the
 * integrator shim's state grows an `at`, add one line to CFT_FD_BLOBS.
 *
 * NOTE 2. REBOUND keeps each seven-level array as ONE allocation of
 * 7*N3 doubles and archives it as one field with element_size 7*8
 * (dpcast() in integrator_ias15.c slices it level-major). ROADMAP.md's
 * struct keeps seven separate pointers, which cannot be one REB_POINTER
 * field, so each level gets its own name. The disk layout is the same
 * bytes in the same order; only the field boundaries differ.
 *
 * NOTE 3. n_elem cannot be trusted to survive a load, and this is the
 * one place where ROADMAP.md's "REB_POINTER with element_size set to
 * the wide width" meets the code and loses. REBOUND's REB_POINTER read
 * path writes size_data/element_size into whatever offset_N names -
 * here that is n_elem - and element_size is fixed in the descriptor at
 * compile time while W is a property of the run. On load the list in
 * force is whichever one the integrator was REGISTERED with (the
 * binary64 one), because reb_simulation_set_integrator installs it when
 * "integrator.name" is read and nothing re-points it before the blobs
 * arrive. So every blob read leaves n_elem = size_data/8 behind.
 *
 * Putting cft_n_elem LAST in the list repairs that for a whole
 * snapshot - it is written last and read last, so the true count
 * overwrites the wrong ones - but NOT for an appended snapshot, which
 * REBOUND stores as a diff and which therefore omits cft_n_elem
 * precisely because it did not change. That failure is real and was
 * caught by gate 1: binary128 and binary256 restarts from snapshot 1
 * came back with n_elem = 18 and 36 for a 9-element state, while
 * binary64 passed because 72/8 happens to be 9. docs/VALIDATION.md has
 * the verbatim output.
 *
 * The fix is not to lean on the mechanism at all: cft_archive_probe
 * reads cft_n_elem and the blob lengths FROM THE FILE, checks they
 * agree with each other, with cft_format and with 3N (or 3NE), and
 * cft_archive_finish_load then writes the true count into the state.
 * The blob CONTENTS are never at risk - the read path freads exactly
 * size_data bytes into a buffer it reallocates to that size - so this
 * costs nothing but the call. It does mean a caller who loads an
 * archive with plain reb_simulation_create_from_file and never calls
 * cft_archive_finish_load holds a state whose n_elem is wrong, and the
 * next save would write blobs of the wrong length. Both entry points
 * here do the repair; nothing else should be used.
 *
 * The list keeps element_size = W anyway, because that is what makes
 * the WRITE side correct (size_data = n_elem*W), and cft_n_elem stays
 * last so that a whole-snapshot load is self-consistent even before the
 * repair runs.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <inttypes.h>

#include "rebound.h"
#include "binarydata.h"
#include "simulationarchive.h"

#include "cft_archive.h"

/* ------------------------------------------------------------------ */
/* the descriptor lists                                               */
/* ------------------------------------------------------------------ */

#include "cft_ias15_fields.h"

/* the number of wide blobs: 6 flat arrays + 6 seven-level arrays */

const struct reb_binarydata_field_descriptor *cft_archive_descriptor_list(int format){
    switch (format){
        case CFT_FP64:  return cft_fd_fp64;
        case CFT_FP128: return cft_fd_fp128;
        case CFT_FP256: return cft_fd_fp256;
        default:        return NULL;
    }
}

/* A list that has drifted writes the wrong bytes silently, so check its
 * shape rather than trusting it: 48 blobs then 11 scalars, every name
 * cft_-prefixed and unique, every blob's element_size the format's
 * width and its offset_N the one member the read path is allowed to
 * scribble on, and cft_n_elem last (NOTE 3). Returns 0 if all three
 * lists are sound, otherwise the number of complaints, each printed. */
int cft_archive_selftest(void){
    int bad = 0;
    const int fmts[3] = { CFT_FP64, CFT_FP128, CFT_FP256 };
    for (int fi = 0; fi < 3; fi++){
        const struct reb_binarydata_field_descriptor *l = cft_archive_descriptor_list(fmts[fi]);
        size_t w = cft_format_size((cft_format)fmts[fi]);
        int n = 0, blobs = 0, scalars = 0;
        for (; l[n].name[0]; n++){
            const struct reb_binarydata_field_descriptor *f = &l[n];
            if (strncmp(f->name, "cft_", 4)){
                fprintf(stderr, "cft_archive_selftest: \"%s\" is not cft_-prefixed\n", f->name); bad++;
            }
            for (int m = 0; m < n; m++)
                if (!strcmp(l[m].name, f->name)){
                    fprintf(stderr, "cft_archive_selftest: duplicate name \"%s\"\n", f->name); bad++;
                }
            if (f->dtype == REB_POINTER){
                blobs++;
                if (f->element_size != w){
                    fprintf(stderr, "cft_archive_selftest: %s element_size %zu, want %zu\n",
                            f->name, f->element_size, w); bad++;
                }
                if (f->offset_N != CFT_OFF(n_elem)){
                    fprintf(stderr, "cft_archive_selftest: %s offset_N is not n_elem\n", f->name); bad++;
                }
            }else scalars++;
        }
        if (blobs != CFT_N_BLOBS){
            fprintf(stderr, "cft_archive_selftest: %d blobs at %s, want %d\n",
                    blobs, cft_format_name((cft_format)fmts[fi]), CFT_N_BLOBS); bad++;
        }
        if (scalars != 11){
            fprintf(stderr, "cft_archive_selftest: %d scalars, want 11\n", scalars); bad++;
        }
        if (n < 1 || strcmp(l[n-1].name, "cft_n_elem")){
            fprintf(stderr, "cft_archive_selftest: last field is \"%s\", must be cft_n_elem\n",
                    n ? l[n-1].name : "(empty)"); bad++;
        }
    }
    return bad;
}

const char *cft_archive_status_str(enum cft_archive_status s){
    switch (s){
        case CFT_ARCHIVE_EXACT:           return "restored exactly";
        case CFT_ARCHIVE_PROMOTED:        return "promoted from binary64";
        case CFT_ARCHIVE_FORMAT_MISMATCH: return "refused: format mismatch";
        case CFT_ARCHIVE_COUNT_MISMATCH:  return "refused: element count mismatch";
        case CFT_ARCHIVE_NO_FILE:         return "no such archive";
        case CFT_ARCHIVE_BAD:             return "refused: unreadable or inconsistent archive";
    }
    return "?";
}

int cft_archive_bind(struct reb_simulation *r){
    if (!r || !r->integrator.state) return -1;
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    const struct reb_binarydata_field_descriptor *l = cft_archive_descriptor_list(s->format);
    if (!l) return -1;
    r->integrator.callbacks.field_descriptor_list = l;
    return 0;
}

/* ------------------------------------------------------------------ */
/* allocation                                                         */
/* ------------------------------------------------------------------ */

/* Is this the tail of a wide blob's name? Answered from the descriptor
 * list rather than from a second list of names kept in step by hand:
 * a blob is exactly a REB_POINTER field called cft_<tag>. */
static int is_blob_tag(const char *tag){
    char want[64];
    const struct reb_binarydata_field_descriptor *f;
    if (snprintf(want, sizeof want, "cft_%s", tag) >= (int)sizeof want) return 0;
    for (f = cft_fd_fp64; f->name[0]; f++)
        if (f->dtype == REB_POINTER && strcmp(f->name, want) == 0) return 1;
    return 0;
}

unsigned char **cft_archive_state_blob(struct cft_ias15_state *s, int i){
    /* the same order as CFT_FD_BLOBS, so a caller can walk them all */
    switch (i){
        case 0: return &s->x0;
        case 1: return &s->v0;
        case 2: return &s->a0;
        case 3: return &s->csx;
        case 4: return &s->csv;
        case 5: return &s->csa0;
        default: break;
    }
    if (i == 48) return &s->x;      /* appended: see CFT_FD_BLOBS */
    if (i == 49) return &s->v;
    i -= 6;
    if (i < 0 || i >= 42) return NULL;
    switch (i / 7){
        case 0: return &s->g  [i % 7];
        case 1: return &s->b  [i % 7];
        case 2: return &s->csb[i % 7];
        case 3: return &s->e  [i % 7];
        case 4: return &s->br [i % 7];
        case 5: return &s->er [i % 7];
    }
    return NULL;
}

int cft_archive_state_alloc(struct cft_ias15_state *s, size_t n_elem){
    size_t w = cft_ias15_state_width(s);
    if (!w) return -1;
    for (int i = 0; i < CFT_N_BLOBS; i++){
        unsigned char **p = cft_archive_state_blob(s, i);
        free(*p);
        *p = calloc(n_elem ? n_elem : 1, w);   /* +0 in every format */
        if (!*p) return -1;
    }
    s->n_elem = n_elem;
    return 0;
}

void cft_archive_state_free(struct cft_ias15_state *s){
    for (int i = 0; i < CFT_N_BLOBS; i++){
        unsigned char **p = cft_archive_state_blob(s, i);
        free(*p);
        *p = NULL;
    }
    s->n_elem = 0;
}

/* ------------------------------------------------------------------ */
/* binary64 <-> wide                                                  */
/* ------------------------------------------------------------------ */

int cft_archive_promote_doubles(cft_device *dev, int fmt, const double *src,
                                unsigned char *dst, size_t n){
    if (n == 0) return 0;
    if (fmt == CFT_FP64){ memcpy(dst, src, n * sizeof(double)); return 0; }
    cft_device *own = NULL;
    if (!dev){
        if (cft_open(NULL, 0, &own) != CFT_OK) return -1;
        dev = own;
    }
    cft_status st = cft_convert(dev, CFT_FP64, (cft_format)fmt, CFT_RNE,
                                src, dst, n, NULL);
    if (own) cft_close(own);
    return st == CFT_OK ? 0 : -1;
}

int cft_archive_demote_doubles(cft_device *dev, int fmt, const unsigned char *src,
                               double *dst, size_t n){
    if (n == 0) return 0;
    if (fmt == CFT_FP64){ memcpy(dst, src, n * sizeof(double)); return 0; }
    cft_device *own = NULL;
    if (!dev){
        if (cft_open(NULL, 0, &own) != CFT_OK) return -1;
        dev = own;
    }
    cft_status st = cft_convert(dev, (cft_format)fmt, CFT_FP64, CFT_RNE,
                                src, dst, n, NULL);
    if (own) cft_close(own);
    return st == CFT_OK ? 0 : -1;
}

/* ------------------------------------------------------------------ */
/* the probe: what one snapshot says about itself                     */
/* ------------------------------------------------------------------ */

/* The prefix an integrator field carries: "integrator.<name>." */
static int cft_field_suffix(const char *name, const char *integrator, const char **suffix){
    size_t l = strlen(integrator);
    if (strncmp(name, "integrator.", 11) != 0) return 0;
    if (!l || strncmp(name + 11, integrator, l) != 0) return 0;
    if (name[11 + l] != '.') return 0;
    *suffix = name + 11 + l + 1;
    return 1;
}

static int cft_probe_one(FILE *f, uint64_t offset, struct cft_archive_info *out){
    struct reb_binarydata_field field;
    char name[REB_STRING_SIZE_MAX];
    if (fseek(f, (long)offset, SEEK_SET)) return -1;
    while (1){
        if (fread(&field, sizeof field, 1, f) != 1) return -1;
        if (field.size_name == reb_binarydata_header){
            if (fseek(f, 64 - (long)sizeof field, SEEK_CUR)) return -1;
            continue;
        }
        if (field.size_name == 0 || field.size_name > REB_STRING_SIZE_MAX) return -1;
        if (fread(name, (size_t)field.size_name, 1, f) != 1) return -1;
        if (strcmp(name, "end") == 0) return 0;

        if (strcmp(name, "integrator.name") == 0){
            size_t n = (size_t)field.size_data;
            if (n == 0 || n > sizeof out->integrator) return -1;
            if (fread(out->integrator, n, 1, f) != 1) return -1;
            out->integrator[sizeof out->integrator - 1] = '\0';
            continue;
        }

        const char *suffix = NULL;
        if (out->integrator[0] && cft_field_suffix(name, out->integrator, &suffix)
                && strncmp(suffix, "cft_", 4) == 0){
            const char *tag = suffix + 4;
            out->has_cft = 1;
            out->n_cft_fields++;
            if (strcmp(tag, "format") == 0){
                int v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                out->format = v;
                continue;
            }
            if (strcmp(tag, "n_elem") == 0 || strcmp(tag, "E") == 0){
                size_t v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                if (tag[0] == 'n') out->n_elem = (uint64_t)v; else out->E = (uint64_t)v;
                continue;
            }
            if (strcmp(tag, "abi_0") == 0 || strcmp(tag, "abi_1") == 0){
                if (field.size_data != 8) return -1;
                if (fread(out->abi + (tag[4] == '0' ? 0 : 8), 8, 1, f) != 1) return -1;
                out->abi[16] = '\0';
                continue;
            }
            if (strcmp(tag, "constants_digest") == 0){
                uint64_t v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                out->constants_digest = v;
                continue;
            }
            /* a wide blob: only its length matters here */
            if (is_blob_tag(tag)){
                if (out->blob_bytes == 0) out->blob_bytes = (uint64_t)field.size_data;
                else if (out->blob_bytes != (uint64_t)field.size_data) out->blob_lengths_agree = 0;
                out->n_blobs_seen++;
            }
        }
        if (fseek(f, (long)field.size_data, SEEK_CUR)) return -1;
    }
}

int cft_archive_probe(const char *filename, int64_t snapshot,
                      struct cft_archive_info *out){
    memset(out, 0, sizeof *out);
    out->format = -1;
    out->blob_lengths_agree = 1;

    struct reb_simulationarchive *sa = reb_simulationarchive_create_from_file(filename);
    if (!sa) return -1;
    if (!sa->inf){ reb_simulationarchive_free(sa); return -1; }
    out->n_snapshots = sa->nblobs;
    if (snapshot < 0) snapshot += sa->nblobs;
    if (snapshot < 0 || snapshot >= sa->nblobs){ reb_simulationarchive_free(sa); return -1; }

    FILE *f = fopen(filename, "rb");
    if (!f){ reb_simulationarchive_free(sa); return -1; }
    int err = cft_probe_one(f, sa->offset[0], out);
    /* A snapshot other than 0 is a DIFF over snapshot 0 (see
     * reb_binarydata_diff in simulationarchive.c's append path), so it
     * is read as an overlay, exactly as REBOUND's loader applies it. */
    if (!err && snapshot != 0){
        struct cft_archive_info over;
        memset(&over, 0, sizeof over);
        over.format = -1;
        over.blob_lengths_agree = 1;
        memcpy(over.integrator, out->integrator, sizeof over.integrator);
        err = cft_probe_one(f, sa->offset[snapshot], &over);
        if (!err){
            if (over.integrator[0]) memcpy(out->integrator, over.integrator, sizeof out->integrator);
            if (over.format >= 0)   out->format = over.format;
            if (over.n_elem)        out->n_elem = over.n_elem;
            if (over.E)             out->E = over.E;
            if (over.abi[0])        memcpy(out->abi, over.abi, sizeof out->abi);
            if (over.constants_digest) out->constants_digest = over.constants_digest;
            if (over.blob_bytes){
                if (out->blob_bytes && out->blob_bytes != over.blob_bytes) out->blob_lengths_agree = 0;
                out->blob_bytes = over.blob_bytes;
            }
            if (!over.blob_lengths_agree) out->blob_lengths_agree = 0;
            out->has_cft |= over.has_cft;
            out->n_cft_fields += over.n_cft_fields;
        }
    }
    fclose(f);
    reb_simulationarchive_free(sa);
    return err;
}

/* ------------------------------------------------------------------ */
/* loading                                                            */
/* ------------------------------------------------------------------ */

static const char *fmt_name(int f){
    const char *n = cft_format_name((cft_format)f);
    return n ? n : "(not a libcft format)";
}

static void refuse(const char *filename, int64_t snapshot, const char *what,
                   const char *detail){
    fprintf(stderr,
        "cft-rebound: refusing to load \"%s\" snapshot %" PRId64 ": %s\n"
        "  %s\n"
        "  The wide state is not reinterpretable across this difference and\n"
        "  will not be silently narrowed or reshaped. Load the archive with\n"
        "  the settings it was written with, or start a new run.\n",
        filename, snapshot, what, detail);
}

/* Promote REBOUND's binary64 state into the wide arrays. */
static enum cft_archive_status promote(struct reb_simulation *r, int want, cft_device *dev){
    size_t N  = r->N;
    size_t N3 = 3 * N;

    /* Take a copy of REBOUND's own IAS15 series before set_integrator
     * frees it. The particles are the authority for x0/v0: REBOUND's
     * ias15->x0 holds the START of the last completed step. */
    double eps = 1e-9, min_dt = 0.0;
    int adaptive_mode = 2;                 /* PRS23, REBOUND's default since 01/2024 */
    double *keep[12];                      /* a0 csx csv csa0 then g b csb e br er */
    int have_ias15 = 0;
    for (int i = 0; i < 12; i++) keep[i] = NULL;

    if (r->integrator.name && strcmp(r->integrator.name, "ias15") == 0 && r->integrator.state){
        struct reb_integrator_ias15_state *ia =
            (struct reb_integrator_ias15_state*)r->integrator.state;
        eps = ia->epsilon; min_dt = ia->min_dt; adaptive_mode = (int)ia->adaptive_mode;
        if (ia->N_allocated >= N3){
            have_ias15 = 1;
            double *flat[4] = { ia->a0, ia->csx, ia->csv, ia->csa0 };
            for (int i = 0; i < 4; i++){
                keep[i] = malloc(N3 * sizeof(double));
                if (keep[i] && flat[i]) memcpy(keep[i], flat[i], N3 * sizeof(double));
                else if (keep[i]) memset(keep[i], 0, N3 * sizeof(double));
            }
            double *lvl[6] = { ia->g, ia->b, ia->csb, ia->e, ia->br, ia->er };
            for (int i = 0; i < 6; i++){
                keep[4 + i] = malloc(7 * N3 * sizeof(double));
                if (keep[4 + i] && lvl[i]) memcpy(keep[4 + i], lvl[i], 7 * N3 * sizeof(double));
                else if (keep[4 + i]) memset(keep[4 + i], 0, 7 * N3 * sizeof(double));
            }
        }
    }

    void *sp = reb_simulation_set_integrator(r, CFT_IAS15_INTEGRATOR_NAME);
    if (!sp){
        for (int i = 0; i < 12; i++) free(keep[i]);
        return CFT_ARCHIVE_BAD;
    }
    struct cft_ias15_state *s = (struct cft_ias15_state*)sp;
    s->format        = want;
    s->E             = 1;
    s->epsilon       = eps;
    s->min_dt        = min_dt;
    s->adaptive_mode = adaptive_mode;
    if (s->max_iter <= 0) s->max_iter = 12;
    snprintf(s->cft_abi, sizeof s->cft_abi, "%d.%d",
             CFT_ABI_VERSION_MAJOR, CFT_ABI_VERSION_MINOR);
    if (cft_archive_state_alloc(s, N3)){
        for (int i = 0; i < 12; i++) free(keep[i]);
        return CFT_ARCHIVE_BAD;
    }

    double *tmp = malloc(N3 * sizeof(double));
    if (!tmp){
        for (int i = 0; i < 12; i++) free(keep[i]);
        return CFT_ARCHIVE_BAD;
    }
    for (size_t i = 0; i < N; i++){ tmp[3*i] = r->particles[i].x; tmp[3*i+1] = r->particles[i].y; tmp[3*i+2] = r->particles[i].z; }
    cft_archive_promote_doubles(dev, want, tmp, s->x0, N3);
    for (size_t i = 0; i < N; i++){ tmp[3*i] = r->particles[i].vx; tmp[3*i+1] = r->particles[i].vy; tmp[3*i+2] = r->particles[i].vz; }
    cft_archive_promote_doubles(dev, want, tmp, s->v0, N3);

    if (have_ias15){
        unsigned char *flat[4] = { s->a0, s->csx, s->csv, s->csa0 };
        for (int i = 0; i < 4; i++)
            cft_archive_promote_doubles(dev, want, keep[i], flat[i], N3);
        unsigned char **lvl[6] = { s->g, s->b, s->csb, s->e, s->br, s->er };
        for (int i = 0; i < 6; i++)
            for (int k = 0; k < 7; k++)
                cft_archive_promote_doubles(dev, want, keep[4 + i] + (size_t)k * N3,
                                            lvl[i][k], N3);
    }else{
        for (size_t i = 0; i < N; i++){ tmp[3*i] = r->particles[i].ax; tmp[3*i+1] = r->particles[i].ay; tmp[3*i+2] = r->particles[i].az; }
        cft_archive_promote_doubles(dev, want, tmp, s->a0, N3);
    }

    free(tmp);
    for (int i = 0; i < 12; i++) free(keep[i]);
    cft_archive_bind(r);
    return CFT_ARCHIVE_PROMOTED;
}

enum cft_archive_status cft_archive_finish_load(struct reb_simulation *r,
                                                const struct cft_archive_info *info,
                                                int want_format, cft_device *dev){
    if (!r) return CFT_ARCHIVE_BAD;

    if (!info->has_cft){
        int want = (want_format >= 0) ? want_format : CFT_FP64;
        if (!cft_format_size((cft_format)want)) return CFT_ARCHIVE_BAD;
        return promote(r, want, dev);
    }

    if (!r->integrator.name || strcmp(r->integrator.name, CFT_IAS15_INTEGRATOR_NAME) != 0)
        return CFT_ARCHIVE_BAD;   /* the shim was not registered before the load */
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    if (!s) return CFT_ARCHIVE_BAD;

    if (!cft_format_size((cft_format)info->format)) return CFT_ARCHIVE_BAD;
    if (want_format >= 0 && info->format != want_format) return CFT_ARCHIVE_FORMAT_MISMATCH;
    if (s->format != info->format) return CFT_ARCHIVE_BAD;

    size_t w = cft_format_size((cft_format)info->format);
    if (info->n_blobs_seen != CFT_N_BLOBS)            return CFT_ARCHIVE_COUNT_MISMATCH;
    if (!info->blob_lengths_agree)                    return CFT_ARCHIVE_COUNT_MISMATCH;
    if (!info->n_elem)                                return CFT_ARCHIVE_COUNT_MISMATCH;
    if (info->blob_bytes != info->n_elem * (uint64_t)w) return CFT_ARCHIVE_COUNT_MISMATCH;
    {
        uint64_t E = info->E ? info->E : 1;
        if (info->n_elem != 3ull * (uint64_t)r->N * E) return CFT_ARCHIVE_COUNT_MISMATCH;
        s->E = (size_t)E;
    }

    /* The file is the authority for the element count; the loaded state
     * is not, because the reader derived it with the wrong width. See
     * NOTE 3 at the top of this file. */
    s->n_elem = (size_t)info->n_elem;

    /* Did the load actually happen? Everything above was read out of
     * the file by cft_archive_probe() and says nothing about what
     * reached memory. REBOUND allocates a REB_POINTER field's buffer
     * only when it reads that field, so a NULL blob here means the
     * loader never saw it - a descriptor list that does not match the
     * file. That is exactly the state this function used to certify
     * as "restored exactly" (docs/VALIDATION.md entry 25), and a
     * wrong answer reported as a right one is the one outcome this
     * repository refuses. */
    {
        int i;
        for (i = 0; i < CFT_N_BLOBS; i++){
            unsigned char **b = cft_archive_state_blob(s, i);
            if (!b || !*b){
                fprintf(stderr,
                        "cft-rebound: the archive names %d wide blobs and blob %d did not "
                        "reach memory. The integrator was registered with a field "
                        "descriptor list that does not match this file, so the state "
                        "was not restored; refusing rather than continuing from "
                        "defaults.\n", info->n_blobs_seen, i);
                return CFT_ARCHIVE_BAD;
            }
        }
    }

    if (cft_archive_bind(r)) return CFT_ARCHIVE_BAD;
    return CFT_ARCHIVE_EXACT;
}

struct reb_simulation *cft_archive_load(const char *filename, int64_t snapshot,
                                        int want_format, cft_device *dev,
                                        enum cft_archive_status *status){
    enum cft_archive_status st;
    struct cft_archive_info info;

    if (cft_archive_probe(filename, snapshot, &info)){
        FILE *f = fopen(filename, "rb");
        if (f){ fclose(f); st = CFT_ARCHIVE_BAD; }
        else     st = CFT_ARCHIVE_NO_FILE;
        if (st == CFT_ARCHIVE_BAD)
            refuse(filename, snapshot, "its field table could not be read",
                   "The file is not a Simulationarchive, or the snapshot is corrupt.");
        if (status) *status = st;
        return NULL;
    }

    /* Refuse a format disagreement before anything is loaded: the state
     * on disk is the run's, and reinterpreting binary256 bytes as
     * binary64 is a different simulation, not a lossy one. */
    if (info.has_cft && want_format >= 0 && info.format != want_format){
        char detail[512];
        snprintf(detail, sizeof detail,
                 "the archive's cft_format is %s (%zu bytes an element); this run is "
                 "configured for %s (%zu bytes an element).",
                 fmt_name(info.format), cft_format_size((cft_format)info.format),
                 fmt_name(want_format), cft_format_size((cft_format)want_format));
        refuse(filename, snapshot, "it was written at a different precision", detail);
        if (status) *status = CFT_ARCHIVE_FORMAT_MISMATCH;
        return NULL;
    }

    struct reb_simulation *r = reb_simulation_create_from_file((char*)filename, snapshot);
    if (!r){
        if (status) *status = CFT_ARCHIVE_NO_FILE;
        return NULL;
    }

    st = cft_archive_finish_load(r, &info, want_format, dev);
    if (st < 0){
        char detail[512];
        switch (st){
            case CFT_ARCHIVE_COUNT_MISMATCH:
                snprintf(detail, sizeof detail,
                    "cft_n_elem is %" PRIu64 " and cft_E is %" PRIu64 ", but the archive holds "
                    "%d of %d wide blobs and each is %" PRIu64 " bytes, against %zu N particles "
                    "at %zu bytes an element.",
                    info.n_elem, info.E ? info.E : 1, info.n_blobs_seen, CFT_N_BLOBS,
                    info.blob_bytes, r->N,
                    cft_format_size((cft_format)info.format));
                refuse(filename, snapshot, "its element count does not add up", detail);
                break;
            case CFT_ARCHIVE_FORMAT_MISMATCH:
                snprintf(detail, sizeof detail,
                    "the archive's cft_format is %s; this run is configured for %s.",
                    fmt_name(info.format), fmt_name(want_format));
                refuse(filename, snapshot, "it was written at a different precision", detail);
                break;
            default:
                snprintf(detail, sizeof detail,
                    "the archive names integrator \"%s\" and carries %d cft_ fields; "
                    "the loaded simulation's integrator is \"%s\".",
                    info.integrator, info.n_cft_fields,
                    r->integrator.name ? r->integrator.name : "(none)");
                refuse(filename, snapshot, "its integrator state could not be adopted", detail);
                break;
        }
        reb_simulation_free(r);
        if (status) *status = st;
        return NULL;
    }

    if (status) *status = st;
    return r;
}
