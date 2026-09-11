/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The body-count cap, once.
 *
 * It was written twice - CFT_MAX_BODIES in src/ias15_cft.c and
 * CFT_REBOUND_MAX_BODIES in src/cft_rebound_run.c - with a comment on
 * the second saying "kept in step with" the first. A number kept in
 * step by a comment is a number that will drift; this repository has
 * the receipts (ROADMAP.md, the last section, on the blob count).
 *
 * This header exists rather than the number living in cft_supported.h
 * because the engine must not include rebound.h: src/ias15_cft.c builds
 * both as the standalone program and, with -DIAS15_CFT_LIBRARY, as a
 * library half that knows nothing about REBOUND. So the shared fact
 * needs a home neither side has to reach through the other for, and it
 * depends on nothing.
 *
 * The cap is set where the worst case, binary256, still fits an
 * ordinary workstation: 1024 bodies want 4.4 GB and 2048 would want
 * 17.7 GB. Time stops you well before memory does - README has the
 * table. The engine is the authority and refuses on its own; the two
 * callers check first so that the refusal arrives before any work, and
 * on the subprocess path before a temporary file is written.
 */
#ifndef IAS15_LIMITS_H
#define IAS15_LIMITS_H

#define CFT_MAX_BODIES 1024

/* The seed of the FNV-STYLE hash this project digests wide bytes with -
 * the constants at ias15_cft.c's digest_constants(), the state at
 * gate_real.c's wide_digest(). FNV-style, not FNV-1a: the prime is
 * FNV's, the offset basis is FNV-1a's 14695981039346656037 with its
 * last decimal digit missing. Kept rather than corrected, because the
 * constants digest is ARCHIVED and correcting the seed would change
 * every stored value to no purpose - nothing compares it against an
 * independent implementation. It is here rather than in either file
 * because it was in both, three hours after this repository wrote down
 * the rule against that. */
#define CFT_DIGEST_SEED 1469598103934665603ULL

#endif /* IAS15_LIMITS_H */
