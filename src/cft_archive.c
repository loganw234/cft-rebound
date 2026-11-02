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
 *   -                                       cft_x            REB_POINTER  W   NOTE 4
 *   -                                       cft_v            REB_POINTER  W   NOTE 4
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
 *
 * NOTE 5. cft_alias_g0 .. cft_alias_br6, cft_alias_csx and
 * cft_alias_csv are the arrays REBOUND holds at N_allocated rather than
 * at the live 3N:
 *
 *   REBOUND                                  here
 *   x0 v0 a0 csx csv csa0  at N_allocated    cft_alias_csx, cft_alias_csv
 *   g b csb e br er        at N_allocated    cft_alias_<family><level>
 *
 * REBOUND's reb_integrator_ias15_alloc() reallocates only when 3N
 * exceeds N_allocated, so below the mark nothing moves and nothing is
 * zeroed. Two consequences, and this port had to reproduce both:
 *
 *  - the six coefficient families share one allocation apiece that
 *    dpcast() re-slices at the CURRENT 3N, so a shrink strands a tail
 *    no stride reaches and a regrow under the mark reads it straight
 *    back. ias15_engine_alias_resize() keeps that in a flat shadow;
 *    cft_alias_<family><level> is that shadow sliced at the MARK's
 *    stride, which is what REBOUND's buffer holds.
 *
 *  - csx and csv keep their values above the live region too, and the
 *    predictor reads them before writing them. Of REBOUND's flat
 *    arrays they are the only two that do: x0, v0, a0 and csa0 are all
 *    written at the top of every attempt, and this port re-promotes x
 *    and v from r->particles whenever the body count changes. Measured
 *    rather than argued - without these two, gate_real's second
 *    sequence differs in the last bit of exactly the coordinates a
 *    stale carry touches.
 *
 * That is the whole of ROADMAP.md's "What the collision work leaves
 * open, across a checkpoint". Nothing carried the tail across a save,
 * and nothing carried N_allocated either: the shim set its high-water
 * mark from r->N when it bound, so a resumed run's mark was the SMALLER
 * count and a regrow zeroed where REBOUND aliases.
 *
 * REBOUND archives no N_allocated field - every one of its pointers
 * names N_allocated as offset_N and the reader recovers it as
 * size_data/element_size. Same here: the mark is
 * cft_hiwater_n_elem/3.
 *
 * They are written only when the mark is ABOVE the live count, which is
 * the only time they hold anything a live blob does not. REBOUND's
 * writer skips a REB_POINTER field of zero length, so an archive from a
 * run that never shrank carries none of them.
 *
 * It is not, however, byte for byte what this project wrote before they
 * existed - the three new SCALARS are unconditional, and measuring it
 * gives 10590 bytes against 10783 and 62 named cft_ fields against 65.
 * The compatibility that does hold is the useful one: the shared names
 * keep their order and their bytes, and an older archive resumes here
 * bit-identically.
 *
 * NOTE 4. cft_x and cft_v have no counterpart in REBOUND's list because
 * REBOUND's live coordinates are in r->particles, which it archives as
 * simulation fields. Here they are the integrator's own: step_attempt()
 * copies x into x0 at the top of a step and leaves the advanced position
 * in x, so after a completed step x0 is the PREVIOUS step's start and x
 * is where the bodies are. They were missing from the first version of
 * this list, which at binary128 and binary256 truncated a checkpoint's
 * position and velocity to binary64 while every gate passed; that is
 * what took the blob count from 48 to 50. docs/VALIDATION.md entry 25
 * has the measurements.
 *
 * NOTE 6. cft_t, cft_t_hi, cft_t_lo, cft_dt and cft_dt_last are the
 * engine's wide clock and steps, and cft_identity is what the producer
 * says it ran on. All six were added on 2026-11-02 and have no
 * counterpart in REBOUND's list, whose t, dt and dt_last_done are the
 * binary64 views of the first, fourth and fifth (and are simulation
 * fields, not integrator ones).
 *
 *   - Their length members count BYTES (element_size 1), so they come
 *     back right through the binary64 list and need none of NOTE 3's
 *     repair. src/cft_ias15_fields.h, CFT_FD_B.
 *   - They are OPTIONAL on read, like NOTE 5's family. Every archive
 *     written before 2026-11-02 lacks them. Its load restores the wide
 *     state exactly and rebuilds the five from binary64, as every load
 *     did then, and it says so: CFT_ARCHIVE_EXACT_CLOCK_STEPS64, and
 *     CFT_PROV_EXACT_CLOCK_STEPS64 in the state.
 *   - They come LAST, before cft_n_elem. An older build stops reading a
 *     snapshot at the first name it does not know, so it loses only
 *     cft_n_elem, which its finish_load repairs from the file. Measured,
 *     in both directions, in tools/check_archive.py's docstring.
 *   - cft_identity differs between builds by design. Anything that
 *     compares whole archives across builds or pins leaves it out.
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

/* The lists themselves are defined once, in src/cft_ias15_fields.c,
 * from the macros in the header above; this file only selects among
 * them. CFT_N_BLOBS is 50: eight flat arrays (x0, v0, a0, csx, csv,
 * csa0 and the live x, v) and six seven-level ones. They were
 * file-scope statics here until 2026-09-10, which meant the registered
 * integrator could not name them - see docs/VALIDATION.md entry 25.
 *
 * CFT_N_ALIAS_BLOBS is a second family, sized by hiwater_n_elem rather
 * than by n_elem and optional on disk - NOTE 5. Every count comes from
 * src/cft_ias15_fields.h; nothing here writes one down. */

const struct reb_binarydata_field_descriptor *cft_archive_descriptor_list(int format){
    switch (format){
        case CFT_FP64:  return cft_fd_fp64;
        case CFT_FP128: return cft_fd_fp128;
        case CFT_FP256: return cft_fd_fp256;
        default:        return NULL;
    }
}

/* A list that has drifted writes the wrong bytes silently, so check its
 * shape rather than trusting it: CFT_N_BLOBS blobs of n_elem,
 * CFT_N_ALIAS_BLOBS of hiwater_n_elem, the clock family and the
 * identity as byte fields (element_size 1, each its own length member),
 * and CFT_N_SCALARS scalars; every name cft_-prefixed and unique, every
 * blob's element_size the format's width and its offset_N one of the
 * length members the read path is allowed to scribble on, the alias and
 * clock descriptors landing in their walkers' order, every field
 * appended on 2026-11-02 after every older one (NOTE 6), and cft_n_elem
 * last (NOTE 3). Returns 0 if all three lists are sound, otherwise the
 * number of complaints, each printed. */
int cft_archive_selftest(void){
    int bad = 0;
    const int fmts[3] = { CFT_FP64, CFT_FP128, CFT_FP256 };
    const size_t apfx = sizeof CFT_ALIAS_PREFIX - 1;
    /* Addresses in a dummy state, so nothing is allocated and nothing
     * is read; only offsets within it are taken. */
    struct cft_ias15_state probe;
    memset(&probe, 0, sizeof probe);
    for (int fi = 0; fi < 3; fi++){
        const struct reb_binarydata_field_descriptor *l = cft_archive_descriptor_list(fmts[fi]);
        size_t w = cft_format_size((cft_format)fmts[fi]);
        int n = 0, blobs = 0, alias = 0, clock = 0, ident = 0, scalars = 0;
        int first_new = -1;          /* the first field appended on 2026-11-02 */
        for (; l[n].name[0]; n++){
            const struct reb_binarydata_field_descriptor *f = &l[n];
            int is_alias_name = !strncmp(f->name, CFT_ALIAS_PREFIX, apfx);
            if (strncmp(f->name, "cft_", 4)){
                fprintf(stderr, "cft_archive_selftest: \"%s\" is not cft_-prefixed\n", f->name); bad++;
            }
            for (int m = 0; m < n; m++)
                if (!strcmp(l[m].name, f->name)){
                    fprintf(stderr, "cft_archive_selftest: duplicate name \"%s\"\n", f->name); bad++;
                }
            /* An older build's reader abandons a snapshot at the first
             * name it does not know, so a field it does know must never
             * follow one it does not - cft_n_elem, which every build
             * repairs from the file, the one exception (NOTE 6). */
            if (f->offset >= CFT_OFF_APPENDED_2026_11_02){
                if (first_new < 0) first_new = n;
            }else if (first_new >= 0 && strcmp(f->name, "cft_n_elem")){
                fprintf(stderr, "cft_archive_selftest: %s predates 2026-11-02 and follows %s, "
                                "which does not; a build from before then would stop "
                                "reading at %s and never see %s\n",
                        f->name, l[first_new].name, l[first_new].name, f->name);
                bad++;
            }
            if (f->dtype == REB_POINTER
                    && (f->offset_N == CFT_OFF(clock_bytes) || f->offset_N == CFT_OFF(identity_bytes))){
                /* The byte fields: their length is counted in bytes, so
                 * element_size is 1 at every format (CFT_FD_B). */
                if (f->element_size != 1){
                    fprintf(stderr, "cft_archive_selftest: %s is a byte field and its "
                                    "element_size is %zu, not 1\n", f->name, f->element_size);
                    bad++;
                }
                if (f->offset_N == CFT_OFF(clock_bytes)){
                    /* In cft_archive_state_clock()'s order, for the reason
                     * the alias family is held to its walker's. */
                    unsigned char **m = cft_archive_state_clock(&probe, clock);
                    size_t want = m ? (size_t)((char*)m - (char*)&probe) : (size_t)-1;
                    if (f->offset != want){
                        fprintf(stderr, "cft_archive_selftest: clock field %d (%s) is at offset "
                                        "%zu, want %zu; CFT_FD_CLOCK and "
                                        "cft_archive_state_clock() are out of step\n",
                                clock, f->name, f->offset, want);
                        bad++;
                    }
                    clock++;
                }else{
                    if (f->offset != CFT_OFF(identity)){
                        fprintf(stderr, "cft_archive_selftest: %s is sized by identity_bytes "
                                        "and is not the identity\n", f->name);
                        bad++;
                    }
                    ident++;
                }
            }else if (f->dtype == REB_POINTER){
                if (f->element_size != w){
                    fprintf(stderr, "cft_archive_selftest: %s element_size %zu, want %zu\n",
                            f->name, f->element_size, w); bad++;
                }
                if (f->offset_N == CFT_OFF(n_elem)){
                    blobs++;
                    if (is_alias_name){
                        fprintf(stderr, "cft_archive_selftest: %s is named for the alias "
                                        "family but its length is n_elem\n", f->name); bad++;
                    }
                }else if (f->offset_N == CFT_OFF(hiwater_n_elem)){
                    /* Each must land on the member cft_archive_state_alias()
                     * hands back at the same index: nothing else checks
                     * that the descriptor order and the member order
                     * agree, and two families swapped would read each
                     * other's coefficients back without a word. */
                    unsigned char **m = cft_archive_state_alias(&probe, alias);
                    size_t want = m ? (size_t)((char*)m - (char*)&probe) : (size_t)-1;
                    if (!is_alias_name){
                        fprintf(stderr, "cft_archive_selftest: %s is sized by hiwater_n_elem "
                                        "and is not named \"%s...\"\n", f->name, CFT_ALIAS_PREFIX);
                        bad++;
                    }
                    if (f->offset != want){
                        fprintf(stderr, "cft_archive_selftest: alias blob %d (%s) is at offset "
                                        "%zu, want %zu; CFT_FD_ALIAS and "
                                        "cft_archive_state_alias() are out of step\n",
                                alias, f->name, f->offset, want);
                        bad++;
                    }
                    alias++;
                }else{
                    fprintf(stderr, "cft_archive_selftest: %s offset_N names none of n_elem, "
                                    "hiwater_n_elem, clock_bytes and identity_bytes\n", f->name);
                    bad++;
                }
            }else scalars++;
        }
        if (blobs != CFT_N_BLOBS){
            fprintf(stderr, "cft_archive_selftest: %d blobs at %s, want %d\n",
                    blobs, cft_format_name((cft_format)fmts[fi]), CFT_N_BLOBS); bad++;
        }
        if (alias != CFT_N_ALIAS_BLOBS){
            fprintf(stderr, "cft_archive_selftest: %d alias blobs at %s, want %d\n",
                    alias, cft_format_name((cft_format)fmts[fi]), CFT_N_ALIAS_BLOBS); bad++;
        }
        if (clock != CFT_N_CLOCK_BLOBS || ident != CFT_N_IDENTITY){
            fprintf(stderr, "cft_archive_selftest: %d clock fields and %d identity at %s, "
                            "want %d and %d\n", clock, ident,
                    cft_format_name((cft_format)fmts[fi]), CFT_N_CLOCK_BLOBS, CFT_N_IDENTITY);
            bad++;
        }
        if (scalars != CFT_N_SCALARS){
            fprintf(stderr, "cft_archive_selftest: %d scalars, want %d\n",
                    scalars, CFT_N_SCALARS); bad++;
        }
        if (n < 1 || strcmp(l[n-1].name, "cft_n_elem")){
            fprintf(stderr, "cft_archive_selftest: last field is \"%s\", must be cft_n_elem\n",
                    n ? l[n-1].name : "(empty)"); bad++;
        }
    }

    /* The walker has to reach every blob the descriptor list names, and
     * reach a different one each time. cft_archive_state_blob() maps an
     * index to a member by position, so a blob added to CFT_FD_BLOBS
     * without a matching case would return NULL here - which is how a
     * stale copy of this walker took three gates down with an access
     * violation rather than a message. The same for the alias family
     * and cft_archive_state_alias(), and for the clock family and
     * cft_archive_state_clock(). */
    {
        enum { NWALK = CFT_N_BLOBS + CFT_N_ALIAS_BLOBS + CFT_N_CLOCK_BLOBS };
        unsigned char **seen[NWALK];
        int k = 0;
        for (int i = 0; i < NWALK; i++){
            int in_alias = i >= CFT_N_BLOBS && i < CFT_N_BLOBS + CFT_N_ALIAS_BLOBS;
            int in_clock = i >= CFT_N_BLOBS + CFT_N_ALIAS_BLOBS;
            int idx = in_clock ? i - CFT_N_BLOBS - CFT_N_ALIAS_BLOBS
                    : in_alias ? i - CFT_N_BLOBS : i;
            unsigned char **b = in_clock ? cft_archive_state_clock(&probe, idx)
                              : in_alias ? cft_archive_state_alias(&probe, idx)
                                         : cft_archive_state_blob(&probe, idx);
            if (!b){
                fprintf(stderr, "cft_archive_selftest: %s blob %d of %d has no member; "
                                "the macro list and the walker disagree\n",
                        in_clock ? "clock" : in_alias ? "alias" : "state", idx,
                        in_clock ? CFT_N_CLOCK_BLOBS : in_alias ? CFT_N_ALIAS_BLOBS : CFT_N_BLOBS);
                bad++; seen[k++] = NULL; continue;
            }
            for (int j = 0; j < k; j++){
                if (seen[j] == b){
                    fprintf(stderr, "cft_archive_selftest: blobs %d and %d are the same "
                                    "member\n", j, i); bad++; break;
                }
            }
            seen[k++] = b;
        }
        if (cft_archive_state_blob(&probe, CFT_N_BLOBS)){
            fprintf(stderr, "cft_archive_selftest: index %d is past the end and still "
                            "returns a member; CFT_N_BLOBS is too small\n", CFT_N_BLOBS);
            bad++;
        }
        if (cft_archive_state_alias(&probe, CFT_N_ALIAS_BLOBS)){
            fprintf(stderr, "cft_archive_selftest: alias index %d is past the end and "
                            "still returns a member; CFT_N_ALIAS_BLOBS is too small\n",
                    CFT_N_ALIAS_BLOBS);
            bad++;
        }
        if (cft_archive_state_clock(&probe, CFT_N_CLOCK_BLOBS)){
            fprintf(stderr, "cft_archive_selftest: clock index %d is past the end and "
                            "still returns a member; CFT_FD_CLOCK is too short\n",
                    CFT_N_CLOCK_BLOBS);
            bad++;
        }
    }
    return bad;
}

const char *cft_archive_status_str(enum cft_archive_status s){
    switch (s){
        case CFT_ARCHIVE_EXACT:           return "restored exactly";
        case CFT_ARCHIVE_PROMOTED:        return "promoted from binary64";
        case CFT_ARCHIVE_EXACT_CLOCK_STEPS64:
            return "restored exactly, its clock and steps from binary64";
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

/* Is this the tail of a wide blob's name, and which family is it in?
 * Answered from the descriptor list rather than from a second list of
 * names kept in step by hand: a blob is exactly a REB_POINTER field
 * called cft_<tag>, and its family is which length member it is sized
 * by - which is the difference that matters and cannot drift from the
 * list the way a name test could. */
#define CFT_BLOB_NONE  0
#define CFT_BLOB_STATE 1    /* n_elem elements: CFT_FD_BLOBS */
#define CFT_BLOB_ALIAS 2    /* hiwater_n_elem elements: CFT_FD_ALIAS */
#define CFT_BLOB_CLOCK 3    /* clock_bytes bytes: CFT_FD_CLOCK (NOTE 6) */
#define CFT_BLOB_IDENT 4    /* identity_bytes bytes: CFT_FD_IDENTITY (NOTE 6) */
static int blob_class(const char *tag){
    char want[64];
    const struct reb_binarydata_field_descriptor *f;
    if (snprintf(want, sizeof want, "cft_%s", tag) >= (int)sizeof want) return CFT_BLOB_NONE;
    for (f = cft_fd_fp64; f->name[0]; f++){
        if (f->dtype != REB_POINTER || strcmp(f->name, want)) continue;
        if (f->offset_N == CFT_OFF(hiwater_n_elem)) return CFT_BLOB_ALIAS;
        if (f->offset_N == CFT_OFF(clock_bytes))    return CFT_BLOB_CLOCK;
        if (f->offset_N == CFT_OFF(identity_bytes)) return CFT_BLOB_IDENT;
        return CFT_BLOB_STATE;
    }
    return CFT_BLOB_NONE;
}

uint64_t cft_archive_hiwater_N(const struct cft_archive_info *info){
    /* Three coordinates a body, which is how REBOUND recovers its own
     * N_allocated from a field length too - see NOTE 5. No alias
     * family means the mark was the live count, which is what its
     * absence records. */
    if (info->hiwater_n_elem) return info->hiwater_n_elem / 3u;
    return info->n_elem / 3u;
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
    /* The two appended at the end of CFT_FD_BLOBS. Positional, and
     * checked by cft_archive_selftest() rather than trusted. */
    if (i == CFT_N_BLOBS - 2) return &s->x;
    if (i == CFT_N_BLOBS - 1) return &s->v;
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

unsigned char **cft_archive_state_alias(struct cft_ias15_state *s, int i){
    /* the same order as CFT_FD_ALIAS: six families of seven levels,
     * then the two flat carries */
    if (i < 0) return NULL;
    if (i < 6 * 7) return &s->alias[i / 7][i % 7];
    if (i == 6 * 7)     return &s->alias_csx;
    if (i == 6 * 7 + 1) return &s->alias_csv;
    return NULL;
}

unsigned char **cft_archive_state_clock(struct cft_ias15_state *s, int i){
    /* the same order as CFT_FD_CLOCK, which cft_archive_selftest() holds
     * this to */
    switch (i){
        case 0: return &s->t;
        case 1: return &s->t_hi;
        case 2: return &s->t_lo;
        case 3: return &s->dt;
        case 4: return &s->dt_last;
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
    /* A freshly allocated state is at its own high-water mark by
     * definition - nothing has been removed from it - so the alias
     * family is empty and its absence is what records that. */
    for (int f = 0; f < CFT_N_ALIAS_BLOBS; f++){
        unsigned char **p = cft_archive_state_alias(s, f);
        if (p){ free(*p); *p = NULL; }
    }
    s->hiwater_n_elem = 0;
    s->n_elem = n_elem;
    return 0;
}

void cft_archive_state_free(struct cft_ias15_state *s){
    for (int i = 0; i < CFT_N_BLOBS; i++){
        unsigned char **p = cft_archive_state_blob(s, i);
        free(*p);
        *p = NULL;
    }
    for (int f = 0; f < CFT_N_ALIAS_BLOBS; f++){
        unsigned char **p = cft_archive_state_alias(s, f);
        if (p){ free(*p); *p = NULL; }
    }
    /* And NOTE 6's members, which a state this function is used on owns
     * whenever it has them: a loaded one's are the loader's. */
    for (int c = 0; c < CFT_N_CLOCK_BLOBS; c++){
        unsigned char **p = cft_archive_state_clock(s, c);
        if (p){ free(*p); *p = NULL; }
    }
    s->clock_bytes = 0;
    free(s->identity);
    s->identity = NULL;
    s->identity_bytes = 0;
    s->hiwater_n_elem = 0;
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
            if (strcmp(tag, "n_elem") == 0 || strcmp(tag, "E") == 0
                    || strcmp(tag, "hiwater_n_elem") == 0){
                size_t v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                if      (!strcmp(tag, "E"))              out->E = (uint64_t)v;
                else if (!strcmp(tag, "hiwater_n_elem")){ out->hiwater_n_elem = (uint64_t)v;
                                                          out->has_hiwater = 1; }
                else                                     out->n_elem = (uint64_t)v;
                continue;
            }
            if (strcmp(tag, "provenance") == 0){
                int v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                out->provenance = v;
                continue;
            }
            if (strcmp(tag, "iterations_max_exceeded") == 0){
                uint64_t v = 0;
                if (field.size_data != sizeof v) return -1;
                if (fread(&v, sizeof v, 1, f) != 1) return -1;
                out->iterations_max_exceeded = v;
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
            /* a wide blob: only its length matters here. The families
             * are counted apart because their lengths are different
             * members - CFT_FD_A and CFT_FD_B in src/cft_ias15_fields.h. */
            switch (blob_class(tag)){
                case CFT_BLOB_STATE:
                    if (out->blob_bytes == 0) out->blob_bytes = (uint64_t)field.size_data;
                    else if (out->blob_bytes != (uint64_t)field.size_data) out->blob_lengths_agree = 0;
                    out->n_blobs_seen++;
                    break;
                case CFT_BLOB_ALIAS:
                    if (out->alias_blob_bytes == 0) out->alias_blob_bytes = (uint64_t)field.size_data;
                    else if (out->alias_blob_bytes != (uint64_t)field.size_data) out->alias_lengths_agree = 0;
                    out->n_alias_blobs_seen++;
                    break;
                case CFT_BLOB_CLOCK:
                    if (out->clock_bytes == 0) out->clock_bytes = (uint64_t)field.size_data;
                    else if (out->clock_bytes != (uint64_t)field.size_data) out->clock_lengths_agree = 0;
                    out->n_clock_seen++;
                    break;
                case CFT_BLOB_IDENT: {
                    /* The one field whose contents the probe reports: what
                     * the writing run said it ran on. Read whole when it
                     * fits, which everything this project writes does. */
                    uint64_t n = field.size_data;
                    size_t take = n < (uint64_t)sizeof out->identity ? (size_t)n
                                                                      : sizeof out->identity - 1;
                    memset(out->identity, 0, sizeof out->identity);
                    if (take && fread(out->identity, take, 1, f) != 1) return -1;
                    out->identity[sizeof out->identity - 1] = '\0';
                    out->identity_bytes = n;
                    out->identity_cut = n > (uint64_t)take;
                    if (fseek(f, (long)(n - (uint64_t)take), SEEK_CUR)) return -1;
                    continue;                  /* consumed: skip the seek below */
                }
                default: break;
            }
        }
        if (fseek(f, (long)field.size_data, SEEK_CUR)) return -1;
    }
}

int cft_archive_probe(const char *filename, int64_t snapshot,
                      struct cft_archive_info *out){
    memset(out, 0, sizeof *out);
    out->format = -1;
    out->provenance = -1;              /* no cft_provenance field at all */
    out->blob_lengths_agree = 1;
    out->alias_lengths_agree = 1;
    out->clock_lengths_agree = 1;

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
        over.provenance = -1;
        over.blob_lengths_agree = 1;
        over.alias_lengths_agree = 1;
        over.clock_lengths_agree = 1;
        memcpy(over.integrator, out->integrator, sizeof over.integrator);
        err = cft_probe_one(f, sa->offset[snapshot], &over);
        if (!err){
            if (over.integrator[0]) memcpy(out->integrator, over.integrator, sizeof out->integrator);
            if (over.format >= 0)   out->format = over.format;
            if (over.provenance >= 0) out->provenance = over.provenance;
            if (over.n_elem)        out->n_elem = over.n_elem;
            if (over.E)             out->E = over.E;
            if (over.abi[0])        memcpy(out->abi, over.abi, sizeof out->abi);
            if (over.constants_digest) out->constants_digest = over.constants_digest;
            if (over.iterations_max_exceeded) out->iterations_max_exceeded = over.iterations_max_exceeded;
            /* A blob's effective length is the overlay's where the
             * overlay carries it and the base's where it does not, so
             * the two lengths have to agree only when the overlay
             * carries SOME of them. When it carries all of them - which
             * is what a particle count change produces, because
             * reb_binarydata_diff emits every field whose size_data
             * differs - the base's lengths are all superseded and a
             * disagreement with them is not one. Reading it as one
             * refused every appended snapshot taken after a removal. */
            if (over.n_blobs_seen == CFT_N_BLOBS){
                out->blob_bytes = over.blob_bytes;
                out->blob_lengths_agree = over.blob_lengths_agree;
            }else if (over.blob_bytes){
                if (out->blob_bytes && out->blob_bytes != over.blob_bytes) out->blob_lengths_agree = 0;
                out->blob_bytes = over.blob_bytes;
                if (!over.blob_lengths_agree) out->blob_lengths_agree = 0;
            }
            /* The alias family. cft_hiwater_n_elem is a scalar like any
             * other and a diff omits it when it has not changed, so an
             * ABSENT field means "unchanged" while a field carrying 0
             * means the mark has risen back to the live count and the
             * base snapshot's blobs are now dead weight. Those are two
             * different facts and has_hiwater is what tells them apart;
             * reading a present 0 as "unchanged" would refuse a
             * perfectly good append. */
            if (over.has_hiwater){
                out->hiwater_n_elem = over.hiwater_n_elem;
                out->has_hiwater = 1;
            }
            if (over.n_alias_blobs_seen == CFT_N_ALIAS_BLOBS){
                out->alias_blob_bytes = over.alias_blob_bytes;
                out->alias_lengths_agree = over.alias_lengths_agree;
            }else if (over.alias_blob_bytes){
                if (out->alias_blob_bytes && out->alias_blob_bytes != over.alias_blob_bytes)
                    out->alias_lengths_agree = 0;
                out->alias_blob_bytes = over.alias_blob_bytes;
                if (!over.alias_lengths_agree) out->alias_lengths_agree = 0;
            }
            if (over.n_alias_blobs_seen > out->n_alias_blobs_seen)
                out->n_alias_blobs_seen = over.n_alias_blobs_seen;
            /* The clock family (NOTE 6). A diff carries only the members
             * that changed - t moves every step, t_lo may not - and the
             * shim publishes all five together and never withdraws them,
             * so the snapshot holds the larger of the two counts, and
             * every length either side carries must agree. A file that
             * breaks that is refused in finish_load, which counts. */
            if (over.n_clock_seen){
                if (out->clock_bytes && out->clock_bytes != over.clock_bytes)
                    out->clock_lengths_agree = 0;
                out->clock_bytes = over.clock_bytes;
                if (!over.clock_lengths_agree) out->clock_lengths_agree = 0;
                if (over.n_clock_seen > out->n_clock_seen) out->n_clock_seen = over.n_clock_seen;
            }
            if (over.identity_bytes){
                memcpy(out->identity, over.identity, sizeof out->identity);
                out->identity_bytes = over.identity_bytes;
                out->identity_cut = over.identity_cut;
            }
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
        enum cft_archive_status st;
        if (!cft_format_size((cft_format)want)) return CFT_ARCHIVE_BAD;
        st = promote(r, want, dev);
        if (st == CFT_ARCHIVE_PROMOTED){
            struct cft_ias15_state *ps = (struct cft_ias15_state*)r->integrator.state;
            if (ps) ps->provenance = CFT_PROV_PROMOTED;
        }
        return st;
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

    /* The alias family, on the same terms and for the same reason. It
     * is OPTIONAL: an archive written by a run whose mark was its live
     * count carries none of it, and neither does any archive written
     * before 2026-09-11, so its absence means a mark of n_elem/3 and
     * an empty shadow - which is what a load did before these fields
     * existed, and still the right answer. NOTE 5. */
    if (info->hiwater_n_elem){
        if (info->n_alias_blobs_seen != CFT_N_ALIAS_BLOBS) return CFT_ARCHIVE_COUNT_MISMATCH;
        if (!info->alias_lengths_agree)                    return CFT_ARCHIVE_COUNT_MISMATCH;
        if (info->alias_blob_bytes != info->hiwater_n_elem * (uint64_t)w)
            return CFT_ARCHIVE_COUNT_MISMATCH;
        /* three coordinates a body, and a high-water mark that cannot
         * be BELOW the live count - it only ever rises, and REBOUND's
         * N_allocated is the same. */
        if (info->hiwater_n_elem % 3u)                     return CFT_ARCHIVE_COUNT_MISMATCH;
        if (info->hiwater_n_elem < info->n_elem)           return CFT_ARCHIVE_COUNT_MISMATCH;
        s->hiwater_n_elem = (size_t)info->hiwater_n_elem;
    }else{
        /* Blobs with no field naming their length is a file that
         * disagrees with itself. Blobs with a length of ZERO are not:
         * that is an appended snapshot whose mark has risen back to the
         * live count, so what the base snapshot wrote is dead weight
         * and is dropped rather than refused. */
        if (info->n_alias_blobs_seen && !info->has_hiwater)
            return CFT_ARCHIVE_COUNT_MISMATCH;
        s->hiwater_n_elem = 0;
    }

    /* The clock family, NOTE 6: OPTIONAL, like the alias family, because
     * every archive written before 2026-11-02 lacks it - but whole when
     * present, and one element of this format's width each. */
    const int has_clock = info->n_clock_seen != 0;
    if (has_clock){
        if (info->n_clock_seen != CFT_N_CLOCK_BLOBS) return CFT_ARCHIVE_COUNT_MISMATCH;
        if (!info->clock_lengths_agree)              return CFT_ARCHIVE_COUNT_MISMATCH;
        if (info->clock_bytes != (uint64_t)w)        return CFT_ARCHIVE_COUNT_MISMATCH;
    }

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
        /* The alias family, by the same argument: the file names six of
         * them, so six must have reached memory. */
        if (s->hiwater_n_elem){
            for (i = 0; i < CFT_N_ALIAS_BLOBS; i++){
                unsigned char **b = cft_archive_state_alias(s, i);
                if (!b || !*b){
                    fprintf(stderr,
                            "cft-rebound: the archive names %d high-water coefficient blobs "
                            "and blob %d did not reach memory. The integrator was registered "
                            "with a field descriptor list that does not match this file, so "
                            "the state was not restored; refusing rather than continuing "
                            "from defaults.\n", info->n_alias_blobs_seen, i);
                    return CFT_ARCHIVE_BAD;
                }
            }
        }
        /* And the clock family. REBOUND's reader leaves each member's
         * byte count in clock_bytes, so that must be the width too. */
        if (has_clock){
            for (i = 0; i < CFT_N_CLOCK_BLOBS; i++){
                unsigned char **b = cft_archive_state_clock(s, i);
                if (!b || !*b || s->clock_bytes != w){
                    fprintf(stderr,
                            "cft-rebound: the archive names %d wide clock and step fields "
                            "and field %d did not reach memory at %zu bytes. The integrator "
                            "was registered with a field descriptor list that does not "
                            "match this file, so the state was not restored; refusing "
                            "rather than continuing from a binary64 clock.\n",
                            info->n_clock_seen, i, w);
                    return CFT_ARCHIVE_BAD;
                }
            }
        }else{
            for (i = 0; i < CFT_N_CLOCK_BLOBS; i++){
                unsigned char **b = cft_archive_state_clock(s, i);
                if (b){ free(*b); *b = NULL; }
            }
            s->clock_bytes = 0;
        }
    }

    /* The identity is the writer's, stated and not checked, and nothing
     * is refused over it. It must still be a string: REBOUND allocated
     * exactly its size_data, and a file is not trusted to have ended it
     * with a NUL. */
    if (s->identity && s->identity_bytes) s->identity[s->identity_bytes - 1] = '\0';
    else{ free(s->identity); s->identity = NULL; s->identity_bytes = 0; }

    if (cft_archive_bind(r)) return CFT_ARCHIVE_BAD;
    /* Last, so that it is about THIS load and not about the run that
     * wrote the file - whose own provenance arrived in the cft_ field
     * and has just been overwritten. cft_archive_probe() is where to
     * read the writer's. An archive with no clock family restores the
     * wide state exactly and rebuilds the clock and steps from binary64
     * (adopt_loaded_state() in src/reb_integrator_cft.c), as every load
     * did before 2026-11-02, and says so rather than "restored exactly". */
    if (!has_clock){
        s->provenance = CFT_PROV_EXACT_CLOCK_STEPS64;
        return CFT_ARCHIVE_EXACT_CLOCK_STEPS64;
    }
    s->provenance = CFT_PROV_EXACT;
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
                    "at %zu bytes an element; cft_hiwater_n_elem is %" PRIu64 " against "
                    "%d of %d high-water blobs of %" PRIu64 " bytes; and it holds %d of %d "
                    "clock and step fields of %" PRIu64 " bytes.",
                    info.n_elem, info.E ? info.E : 1, info.n_blobs_seen, CFT_N_BLOBS,
                    info.blob_bytes, r->N,
                    cft_format_size((cft_format)info.format),
                    info.hiwater_n_elem, info.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS,
                    info.alias_blob_bytes, info.n_clock_seen, CFT_N_CLOCK_BLOBS,
                    info.clock_bytes);
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
