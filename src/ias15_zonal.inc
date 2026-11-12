/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * ZONAL HARMONICS: J2 and J4 of one body of every system, about the z
 * axis, in the run's own format (ROUND4.md, P1).
 *
 * src/ias15_cft.c includes this file just before gravity_body(), which
 * calls zonal_accel() after its fold whenever a problem's
 * `zonal k R J2 J4` line was read; energy_all() calls zonal_potential()
 * and the record's header zonal_record(). read_problem() parses the line
 * into zonal_k, ZR, ZJ2 and ZJ4, exact in the run's format, and keys the
 * corrector schedule with them (sched_key()).
 *
 * Round 4's P0 lands this seam EMPTY: read_problem() then refuses the
 * line by name, so nothing below can run. P1 fills these three bodies and
 * deletes that one refusal. Its brief says what they must compute and the
 * gates that hold them (ROUND4.md, P1).
 */

/* Each body's acceleration from body zonal_k's J2 and J4, added to a, and
 * the reaction on body zonal_k: P1's. */
static void zonal_accel(void){
    (void)zonal_k; (void)ZR; (void)ZJ2; (void)ZJ4;
}

#ifndef IAS15_CFT_LIBRARY
/* The J2 and J4 potential, into each system's epot: P1's. */
static void zonal_potential(V epot){
    (void)epot;
}

/* The record's line naming the zonal term: P1's. */
static void zonal_record(void){
}
#endif /* !IAS15_CFT_LIBRARY */
