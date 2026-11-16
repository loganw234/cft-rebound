/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * ONE SYSTEM'S IAS15 STEP IN THE HARDWARE'S BINARY64 (ROUND5.md, P2).
 *
 * This is src/ias15_cft.c's step under `--arith fma --format fp64`, its
 * default compensated sums (`--cs kahan`), PRS23 (adaptive_mode 2), the
 * gravity of a problem's `active n` line with testparticle_type 0, J2 and
 * J4 from src/ias15_zonal.inc, the exact-time clock and energy_all(),
 * OPERATION FOR OPERATION on the same operands in the same order. The port
 * issues each operation as a vector over every lane; every lane's result
 * depends on its own lane alone, so a body's or a coordinate's chain may
 * run on its own thread, in any order against the others, and give the
 * same bits. What couples lanes is reproduced in the port's order:
 *   - gravity: each body's partners in ascending index, each pair's
 *     arithmetic as gravity_reference() and the gravity-* programs do it,
 *     from the pair's own (i > j) side, the sum from +0;
 *   - J2 and J4's reaction on body k: each J2 reaction in ascending body
 *     order, then each J4 one;
 *   - the corrector's maxima and PRS23's minimum: selections over the
 *     normal values, which no order changes;
 *   - energy_all(): its sums in its own order.
 * Each comment names the port's line it mirrors.
 *
 * ONE SOURCE, TWO EXECUTORS. Each function below is a template on X, the
 * executor: X::each(f) runs f for every body (the host twin a loop, the
 * kernel one thread a body, then a barrier) and X::one(f) runs f once for
 * the system (thread 0, then a barrier). Reductions go through the shared
 * block, sh, written by each() and read by one(); in the kernel that block
 * is __shared__ memory, so a reduction never uses a wave intrinsic and the
 * kernel is the same on 32-wide warps and on RDNA3's wave32 and wave64.
 *
 * WHAT THE PORT COMPUTES AND THIS DOES NOT, a stated limit. The port's
 * gravity computes every pair of the triangle and reads only those its
 * rows name (src/ias15_cft.c:1218): under `active 7` a wave's 1,596 pairs
 * of two test particles are computed and never read. They change no bit;
 * they could raise a flag, and the only one they can raise in a ring is
 * DIVBYZERO from two test particles at exactly one position (a
 * coordinate difference is a multiple of an ulp of the coordinates, so it
 * never underflows when not zero). That flag stops the port and not this.
 * Everything else the port computes and discards is computed here too,
 * for its flags: a rejected step's end-of-step update, the reaction of a
 * test particle, the energy's test-particle terms and every pair of its
 * potential, PRS23's fifth sum, and the prediction that a ratio over 20
 * then zeroes.
 */
#ifndef WAVE_STEP_H
#define WAVE_STEP_H

#include "wave_arith.h"

namespace wave {

#define WN WAVE_MAXN

/* Phases, for a flag's name (wave_main.cpp prints them). */
enum : unsigned { PH_NONE = 0, PH_GRAVITY, PH_ZONAL, PH_STEPSTART, PH_PREDICT, PH_CORRECT, PH_PCERROR,
                  PH_CONTROL, PH_UPDATE, PH_CLOCK, PH_NEXTSTEP, PH_ENERGY, PH_SQRT7, PH_SETUP, PH_N };

/* The zonal coefficients' indices: src/ias15_constants.h's IAS15_ZONAL_*,
 * restated as an enum so device code needs no host header; wave_main.cpp
 * asserts each against the header's own #define. */
enum { ZJ2_NUM = 0, ZJ2_DEN, ZJ2_S1, ZJ2_S0, ZJ2_Z0, ZJ4_NUM, ZJ4_DEN, ZJ4_S2, ZJ4_S1, ZJ4_S0, ZJ4_Z0, ZJ4_Z1,
       ZP2_S1, ZP2_S0, ZP2_DEN, ZP4_S2, ZP4_S1, ZP4_S0, ZP4_DEN, ZN };

/* Everything a run shares: the Gauss-Radau constants of
 * src/ias15_constants.h at binary64, the FMA form's reciprocals, the
 * literals, and the problem's own constants. Filled once by the host. */
struct Consts {
    double h[8], rr[28], c[21], d[21], rinv[28], hf[49], zk[ZN];
    double kposr[8], kvelr[8], khalf, invb;     /* 1/((m+2)(m+3)), 1/(m+2), 1/2, 1/0.25: the port's vdiv */
    double tol;                                 /* 1e-16 * 2^-pc_tol_shift */
    double eps, eps5040;                        /* epsilon, epsilon*5040 */
    double G, soft2, mindt;                     /* softening squared and the step floor: +0, as the port's default */
    double zr, zj2, zj4;
    int N, na, zonal_on, zonal_k, max_iter, adaptive, flag_abort;
    unsigned abort_mask;
};

/* One system's state: the port's vectors for this system's lanes, laid
 * out [field][coordinate][body] so the threads of a workgroup read along a
 * row; then its scalars. */
struct Sys {
    double x[3][WN], v[3][WN], a[3][WN];
    double x0[3][WN], v0[3][WN], a0[3][WN], csx[3][WN], csv[3][WN], at[3][WN];
    double g[7][3][WN], b[7][3][WN], e[7][3][WN], csb[7][3][WN], er[7][3][WN], br[7][3][WN];
    double m[WN];
    double dt, dtlast, tplain, thi, tlo;        /* SDT, SDTLAST, STPLAIN, STHI, STLO */
    double s7;                                  /* sqrt7(epsilon*5040), the port's S7B */
    int s7_ready;                               /* computed at the first step control, as the port's s7_stale */
    long long done, rejected, max_exceeded, pc_total;
    int pc_max;
    unsigned flag_bits, flag_first;             /* every flag raised, and the first by phase */
    int stopped;                                /* an aborting flag: the system took no further step */
    long long attempts;                         /* step attempts made, accepted or not */
    int last_iters, last_capped, last_accepted; /* the last attempt's, for a schedule */
    double last_dt;                             /* the dt the last attempt tried */
};

/* The workgroup's shared block: positions for gravity, the reductions, J2
 * and J4's reactions, and the system's scalars of the attempt. */
struct Shared {
    double xs[3][WN];
    double red[WN], red2[WN];
    double rz2[3][WN], rz4[3][WN];
    double pot[WN * (WN - 1) / 2];              /* the energy's pair terms, (i > j) at i(i-1)/2 + j */
    double D1[8], D2[8];
    double spce, spcelast, dtdone, dtnew, dtnB, ratio, q[8], energy;
    int iters, stop, accepted, skip, use_r, capped;
    unsigned fbits[WN], ffirst[WN];
};

/* binomial(n, k) as an exact double, the port's KBIN (kint(binom(...))) */
__host__ __device__ WAVE_INLINE double binom_d(int n, int k){
    long long r = 1;
    for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
    return (double)r;
}

/* The host twin's executor: each() a loop over the bodies, one() a call. */
struct HostExec {
    int N;
    template <class F> __host__ __device__ void each(F f){ for (int p = 0; p < N; p++) f(p); }
    template <class F> __host__ __device__ void one(F f){ f(); }
};

template <class X>
__host__ __device__ inline void zonal_accel(X &ex, const Consts &K, struct Sys &S, struct Shared &sh, Fl &fl);

/* ---- compensated summation: add_cs(), Kahan's (src/ias15_cft.c:1668) --- */
__host__ __device__ WAVE_INLINE void add_cs(double &p, double &cs, double inp, Fl &fl){
    const double y = w_sub(inp, cs, fl);        /* T3 = inp - cs */
    const double t = w_add(p, y, fl);           /* T4 = p + T3 */
    const double u = w_sub(t, p, fl);           /* T2 = T4 - p */
    cs = w_sub(u, y, fl);                       /* cs = T2 - T3 */
    p = t;
}

/* ---- sqrt7 (src/ias15_cft.c:2251), REBOUND's Newton iteration ---------- */
__host__ __device__ inline double sqrt7(double ain, Fl &fl){
    double aa = ain, scale = 1.0;
    while (aa < 1e-7 && is_normal(aa)){ scale = w_mul(scale, 0.1, fl); aa = w_mul(aa, 1e7, fl); }
    while (1e2 < aa && is_normal(aa)){ scale = w_mul(scale, 10.0, fl); aa = w_mul(aa, 1e-7, fl); }
    double xx = 1.0;
    for (int k = 0; k < 20; k++){
        double x6 = w_mul(xx, xx, fl);
        x6 = w_mul(x6, xx, fl); x6 = w_mul(x6, xx, fl); x6 = w_mul(x6, xx, fl); x6 = w_mul(x6, xx, fl);
        const double q = w_div(aa, x6, fl);
        double d = w_sub(q, xx, fl);
        d = w_div(d, 7.0, fl);
        xx = w_add(xx, d, fl);
    }
    return w_mul(xx, scale, fl);
}

/* ---- gravity (src/ias15_cft.c:1559, gravity_reference() at :1601) ------
 * Body p receives from each partner q in ascending index: under `active
 * na`, testparticle_type 0, an active body from every other active body,
 * a test particle from every active body (build_scatter(), :1244). The
 * pair (i, j), i > j, gives i ((-prefact) m_j) d and j (prefact m_i) d,
 * d = x_i - x_j (the gravity-force program). Then J2 and J4. */
template <class X>
__host__ __device__ inline void gravity(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    ex.each([&](int p){ sh.xs[0][p] = S.x[0][p]; sh.xs[1][p] = S.x[1][p]; sh.xs[2][p] = S.x[2][p]; });
    ex.each([&](int p){
        fl.phase = PH_GRAVITY;
        double ax = 0.0, ay = 0.0, az = 0.0;                    /* ax = 0: the fold starts from +0 */
        const int na = K.na;
        for (int q = 0; q < na; q++){
            if (q == p) continue;
            const int i = p > q ? p : q, j = p > q ? q : p;
            const double dx = w_sub(sh.xs[0][i], sh.xs[0][j], fl);
            const double dy = w_sub(sh.xs[1][i], sh.xs[1][j], fl);
            const double dz = w_sub(sh.xs[2][i], sh.xs[2][j], fl);
            double s = w_mul(dx, dx, fl);                       /* dx*dx + dy*dy + dz*dz + softening2 */
            double t = w_mul(dy, dy, fl);
            s = w_add(s, t, fl);
            t = w_mul(dz, dz, fl);
            s = w_add(s, t, fl);
            s = w_add(s, K.soft2, fl);
            const double r = w_sqrt(s, fl);
            double r3 = w_mul(r, r, fl);                        /* _r*_r*_r */
            r3 = w_mul(r3, r, fl);
            const double pf = w_div(K.G, r3, fl);               /* prefact = G/(_r*_r*_r) */
            double c;
            if (p == i) c = w_mul(w_neg(pf), S.m[j], fl);       /* prefactj = -prefact*m_j */
            else        c = w_mul(pf, S.m[i], fl);              /* prefacti = prefact*m_i */
            ax = w_add(ax, w_mul(c, dx, fl), fl);
            ay = w_add(ay, w_mul(c, dy, fl), fl);
            az = w_add(az, w_mul(c, dz, fl), fl);
        }
        S.a[0][p] = ax; S.a[1][p] = ay; S.a[2][p] = az;
    });
    if (K.zonal_on) zonal_accel(ex, K, S, sh, fl);
}

/* ---- J2 and J4 (src/ias15_zonal.inc) ----------------------------------- */
struct ZGeom { double dx, dy, dz, r2, r, c2, gm; };

/* zonal_geometry(): dx = x_j - x_k ..., r2 = (dx*dx + dy*dy) + dz*dz,
 * r = sqrt(r2), costheta2 = dz*dz/r2, G*M */
__host__ __device__ WAVE_INLINE ZGeom zonal_geometry(const Consts &K, const double *const x[3], const double *m,
                                                     int j, Fl &fl){
    const int k = K.zonal_k;
    ZGeom z;
    z.dx = w_sub(x[0][j], x[0][k], fl);
    z.dy = w_sub(x[1][j], x[1][k], fl);
    z.dz = w_sub(x[2][j], x[2][k], fl);
    const double dxx = w_mul(z.dx, z.dx, fl), dyy = w_mul(z.dy, z.dy, fl), dzz = w_mul(z.dz, z.dz, fl);
    z.r2 = w_add(dxx, dyy, fl);
    z.r2 = w_add(z.r2, dzz, fl);
    z.r = w_sqrt(z.r2, fl);
    z.c2 = w_div(dzz, z.r2, fl);
    z.gm = w_mul(K.G, m[k], fl);
    return z;
}

/* zonal_constant(): (n+1)Jn or Jn, times R nr times, left to right */
__host__ __device__ WAVE_INLINE double zonal_constant(const Consts &K, const double *num, double jv, int nr, Fl &fl){
    double cst = num ? w_mul(*num, jv, fl) : jv;
    for (int i = 0; i < nr; i++) cst = w_mul(cst, K.zr, fl);
    return cst;
}

/* zonal_apply(): the term on body j into (ax, ay, az), and its reaction */
__host__ __device__ WAVE_INLINE void zonal_apply(const ZGeom &z, double gmj, double pf, double f, double fz,
                                                 double &ax, double &ay, double &az, double r[3], Fl &fl){
    const double t = w_mul(z.gm, pf, fl);
    const double tf = w_mul(t, f, fl);
    const double tx = w_mul(tf, z.dx, fl), ty = w_mul(tf, z.dy, fl);
    const double tfz = w_mul(t, fz, fl);
    const double tz = w_mul(tfz, z.dz, fl);
    ax = w_add(ax, tx, fl); ay = w_add(ay, ty, fl); az = w_add(az, tz, fl);
    const double u = w_mul(gmj, pf, fl);
    const double uf = w_mul(u, f, fl);
    r[0] = w_mul(uf, z.dx, fl); r[1] = w_mul(uf, z.dy, fl);
    const double ufz = w_mul(u, fz, fl);
    r[2] = w_mul(ufz, z.dz, fl);
}

/* Body j's J2 and J4 (j != k), from its acceleration so far; its two
 * reactions on k into r2 and r4 */
__host__ __device__ inline void zonal_body(const Consts &K, const double *const x[3], const double *m, int j,
                                           double &ax, double &ay, double &az, double r2[3], double r4[3], Fl &fl){
    const ZGeom z = zonal_geometry(K, x, m, j, fl);
    const double gmj = w_mul(K.G, m[j], fl);
    /* J2: prefac = 3.*J2*R*R/r2/r2/r/2.; fac = 5.*costheta2 - 1.; z's fac - 2. */
    double cst = zonal_constant(K, &K.zk[ZJ2_NUM], K.zj2, 2, fl);
    double pf = w_div(cst, z.r2, fl);
    pf = w_div(pf, z.r2, fl);
    pf = w_div(pf, z.r, fl);
    pf = w_div(pf, K.zk[ZJ2_DEN], fl);
    double f = w_mul(K.zk[ZJ2_S1], z.c2, fl);
    f = w_add(f, K.zk[ZJ2_S0], fl);
    double fz = w_add(f, K.zk[ZJ2_Z0], fl);
    zonal_apply(z, gmj, pf, f, fz, ax, ay, az, r2, fl);
    /* J4: prefac = 5.*J4*R*R*R*R/r2/r2/r2/r/8.; fac = 63.*c2*c2 - 42.*c2 + 3.; z's fac + 12. - 28.*c2 */
    cst = zonal_constant(K, &K.zk[ZJ4_NUM], K.zj4, 4, fl);
    pf = w_div(cst, z.r2, fl);
    pf = w_div(pf, z.r2, fl);
    pf = w_div(pf, z.r2, fl);
    pf = w_div(pf, z.r, fl);
    pf = w_div(pf, K.zk[ZJ4_DEN], fl);
    f = w_mul(K.zk[ZJ4_S2], z.c2, fl);
    f = w_mul(f, z.c2, fl);
    double w = w_mul(K.zk[ZJ4_S1], z.c2, fl);
    f = w_add(f, w, fl);
    f = w_add(f, K.zk[ZJ4_S0], fl);
    fz = w_add(f, K.zk[ZJ4_Z0], fl);
    w = w_mul(K.zk[ZJ4_Z1], z.c2, fl);
    fz = w_add(fz, w, fl);
    zonal_apply(z, gmj, pf, f, fz, ax, ay, az, r4, fl);
}

/* zonal_gives(): body j gives k a reaction iff it is not k and it acts on
 * the active bodies (testparticle_type 0: j < N_active) */
__host__ __device__ WAVE_INLINE int zonal_gives(const Consts &K, int j){ return j != K.zonal_k && j < K.na; }

/* zonal_accel(): each body but k from k's J2 and J4, then k's reaction */
template <class X>
__host__ __device__ inline void zonal_accel(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    if (K.N < 2) return;
    ex.each([&](int p){
        if (p == K.zonal_k) return;
        fl.phase = PH_ZONAL;
        const double *const xs[3] = { sh.xs[0], sh.xs[1], sh.xs[2] };
        double r2[3], r4[3];
        double ax = S.a[0][p], ay = S.a[1][p], az = S.a[2][p];
        zonal_body(K, xs, S.m, p, ax, ay, az, r2, r4, fl);
        S.a[0][p] = ax; S.a[1][p] = ay; S.a[2][p] = az;
        for (int c = 0; c < 3; c++){ sh.rz2[c][p] = r2[c]; sh.rz4[c][p] = r4[c]; }
    });
    ex.one([&]{
        fl.phase = PH_ZONAL;
        const int k = K.zonal_k;
        double ak[3] = { S.a[0][k], S.a[1][k], S.a[2][k] };
        for (int j = 0; j < K.N; j++)
            if (zonal_gives(K, j)) for (int c = 0; c < 3; c++) ak[c] = w_sub(ak[c], sh.rz2[c][j], fl);
        for (int j = 0; j < K.N; j++)
            if (zonal_gives(K, j)) for (int c = 0; c < 3; c++) ak[c] = w_sub(ak[c], sh.rz4[c][j], fl);
        for (int c = 0; c < 3; c++) S.a[c][k] = ak[c];
    });
}

/* zonal_potential() for body j != k: U2 and U4, each to be added to epot
 * in ascending body order (all J2 terms, then all J4 terms) */
__host__ __device__ inline void zonal_pot_body(const Consts &K, const double *const x[3], const double *m, int j,
                                               double &u2, double &u4, Fl &fl){
    const ZGeom z = zonal_geometry(K, x, m, j, fl);
    const double gmm = w_mul(z.gm, m[j], fl);                   /* (G*M)*m_j */
    double cst = zonal_constant(K, nullptr, K.zj2, 2, fl);      /* J2*R*R */
    double pp = w_div(cst, z.r2, fl);
    pp = w_div(pp, z.r, fl);
    pp = w_div(pp, K.zk[ZP2_DEN], fl);
    double pl = w_mul(K.zk[ZP2_S1], z.c2, fl);
    pl = w_add(pl, K.zk[ZP2_S0], fl);
    double pu = w_mul(gmm, pp, fl);
    u2 = w_mul(pu, pl, fl);
    cst = zonal_constant(K, nullptr, K.zj4, 4, fl);             /* J4*R*R*R*R */
    pp = w_div(cst, z.r2, fl);
    pp = w_div(pp, z.r2, fl);
    pp = w_div(pp, z.r, fl);
    pp = w_div(pp, K.zk[ZP4_DEN], fl);
    pl = w_mul(K.zk[ZP4_S2], z.c2, fl);
    pl = w_mul(pl, z.c2, fl);
    const double w = w_mul(K.zk[ZP4_S1], z.c2, fl);
    pl = w_add(pl, w, fl);
    pl = w_add(pl, K.zk[ZP4_S0], fl);
    pu = w_mul(gmm, pp, fl);
    u4 = w_mul(pu, pl, fl);
}

/* --zonal-selftest-from: the term alone at given points, as the port's
 * --zonal-selftest evaluates it (src/ias15_zonal.inc, zonal_selftest()):
 * every acceleration from +0, then the potential from +0 */
template <class X>
__host__ __device__ inline void zonal_selftest_one(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl, double *pot){
    ex.each([&](int p){
        for (int c = 0; c < 3; c++){ sh.xs[c][p] = S.x[c][p]; S.a[c][p] = 0.0; }
    });
    zonal_accel(ex, K, S, sh, fl);
    ex.each([&](int p){
        if (p == K.zonal_k) return;
        const double *const xs[3] = { sh.xs[0], sh.xs[1], sh.xs[2] };
        double u2, u4;
        zonal_pot_body(K, xs, S.m, p, u2, u4, fl);
        sh.rz2[0][p] = u2; sh.rz4[0][p] = u4;
    });
    ex.one([&]{
        double ep = 0.0;
        for (int j = 0; j < K.N; j++) if (j != K.zonal_k && j < K.na) ep = w_add(ep, sh.rz2[0][j], fl);
        for (int j = 0; j < K.N; j++) if (j != K.zonal_k && j < K.na) ep = w_add(ep, sh.rz4[0][j], fl);
        *pot = ep;
    });
}

/* ---- energy_all() (src/ias15_cft.c:3115) --------------------------------
 * 0.5*m*(vx*vx + vy*vy + vz*vz) for every body; -G m_i m_j / r for every
 * pair, (G*m_i)*m_j with i the larger; then the sums in REBOUND's order
 * over the N_interact (= N_active) bodies; the zonal potential; + 0. Every
 * body's and pair's terms are formed, as the port forms them, for their
 * flags; only the sums read the active ones. Into sh.energy. */
template <class X>
__host__ __device__ inline void energy(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    ex.each([&](int p){ sh.xs[0][p] = S.x[0][p]; sh.xs[1][p] = S.x[1][p]; sh.xs[2][p] = S.x[2][p]; });
    ex.each([&](int p){
        fl.phase = PH_ENERGY;
        const double half = w_div(1.0, 2.0, fl);                /* vdiv(half, K1, K2) */
        double bt = w_mul(S.v[0][p], S.v[0][p], fl);
        double bu = w_mul(S.v[1][p], S.v[1][p], fl);
        bt = w_add(bt, bu, fl);
        bu = w_mul(S.v[2][p], S.v[2][p], fl);
        bt = w_add(bt, bu, fl);
        bu = w_mul(half, S.m[p], fl);                           /* vkmul(bu, halfb, mass) */
        sh.red[p] = w_mul(bu, bt, fl);
        /* the pairs whose larger index is p, as the pair list holds them */
        for (int q = 0; q < p; q++){
            const double dx = w_sub(sh.xs[0][p], sh.xs[0][q], fl);
            const double dy = w_sub(sh.xs[1][p], sh.xs[1][q], fl);
            const double dz = w_sub(sh.xs[2][p], sh.xs[2][q], fl);
            double tt = w_mul(dx, dx, fl);
            double uu = w_mul(dy, dy, fl);
            tt = w_add(tt, uu, fl);
            uu = w_mul(dz, dz, fl);
            tt = w_add(tt, uu, fl);
            tt = w_sqrt(tt, fl);
            uu = w_mul(K.G, S.m[p], fl);                        /* G*m_j*m_i, j the larger index */
            uu = w_mul(uu, S.m[q], fl);
            sh.pot[p * (p - 1) / 2 + q] = w_div(uu, tt, fl);
        }
        if (K.zonal_on && p != K.zonal_k){
            const double *const xs[3] = { sh.xs[0], sh.xs[1], sh.xs[2] };
            double u2, u4;
            zonal_pot_body(K, xs, S.m, p, u2, u4, fl);
            sh.rz2[0][p] = u2; sh.rz4[0][p] = u4;
        }
    });
    ex.one([&]{
        fl.phase = PH_ENERGY;
        const int na = K.na;                                    /* N_interact under testparticle_type 0 */
        double ekin = 0.0, epot = 0.0;
        for (int i = 0; i < na; i++) ekin = w_add(ekin, sh.red[i], fl);
        for (int i = 0; i < na; i++)
            for (int j = i + 1; j < na; j++) epot = w_sub(epot, sh.pot[j * (j - 1) / 2 + i], fl);
        if (K.zonal_on && K.N >= 2){
            for (int j = 0; j < K.N; j++) if (j != K.zonal_k && j < na) epot = w_add(epot, sh.rz2[0][j], fl);
            for (int j = 0; j < K.N; j++) if (j != K.zonal_k && j < na) epot = w_add(epot, sh.rz4[0][j], fl);
        }
        double out = w_add(ekin, epot, fl);
        sh.energy = w_add(out, 0.0, fl);                        /* + energy_offset, which is +0 */
    });
}

/* ---- the step criterion: choose_timestep() and dtnew_prs23() ------------
 * (src/ias15_cft.c:2553 and :2470), PRS23 at adaptive_mode 2 */
template <class X>
__host__ __device__ inline void choose_timestep(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    ex.one([&]{ fl.phase = PH_CONTROL; sh.dtnB = w_div(sh.dtdone, 0.25, fl); });   /* dt_done / safety_factor */
    ex.each([&](int p){
        fl.phase = PH_CONTROL;
        double sum[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };            /* each from 0, the three components in order */
        for (int c = 0; c < 3; c++){
            double B[7];
            for (int m = 0; m < 7; m++) B[m] = S.b[m][c][p];
            const double a0 = S.a0[c][p];
            double tmp, sq;
            const double s0 = w_mul(a0, a0, fl);
            tmp = w_add(a0, B[0], fl);
            for (int m = 1; m < 7; m++) tmp = w_add(tmp, B[m], fl);
            const double s1 = w_mul(tmp, tmp, fl);
            tmp = B[0];
            for (int m = 1; m < 7; m++){ sq = w_mul(B[m], (double)(m + 1), fl); tmp = w_add(tmp, sq, fl); }
            const double s2 = w_mul(tmp, tmp, fl);
            tmp = w_mul(B[1], 2.0, fl);                         /* KY4[1] = 1*2 */
            for (int m = 2; m < 7; m++){ sq = w_mul(B[m], (double)(m * (m + 1)), fl); tmp = w_add(tmp, sq, fl); }
            const double s3 = w_mul(tmp, tmp, fl);
            tmp = w_mul(B[2], 6.0, fl);                         /* KY5[2] = 1*2*3 */
            for (int m = 3; m < 7; m++){ sq = w_mul(B[m], (double)((m - 1) * m * (m + 1)), fl); tmp = w_add(tmp, sq, fl); }
            const double s4 = w_mul(tmp, tmp, fl);
            sum[0] = w_add(sum[0], s0, fl); sum[1] = w_add(sum[1], s1, fl); sum[2] = w_add(sum[2], s2, fl);
            sum[3] = w_add(sum[3], s3, fl); sum[4] = w_add(sum[4], s4, fl);
        }
        /* a particle REBOUND skips (a0^2 not normal) is given 1s, and never read */
        const int use = is_normal(sum[0]);
        const double y1 = use ? sum[1] : 1.0, y2 = use ? sum[2] : 1.0, y3 = use ? sum[3] : 1.0;
        const double num = w_mul(2.0, y1, fl);                  /* PRS23: 2*y2 / (y3 + sqrt(y4*y2)) */
        double tmp = w_mul(y3, y1, fl);
        tmp = w_sqrt(tmp, fl);
        const double den = w_add(y2, tmp, fl);
        const double ts2 = w_div(num, den, fl);
        sh.red[p] = (use && is_normal(ts2)) ? ts2 : u2d(0x7ff0000000000000ull);
        (void)sum[4];
    });
    ex.one([&]{
        if (!S.s7_ready){                                       /* dtnew_prs23()'s s7_stale: once, at the first */
            fl.phase = PH_SQRT7;
            S.s7 = sqrt7(K.eps5040, fl);
            S.s7_ready = 1;
        }
        fl.phase = PH_CONTROL;
        double mints2 = u2d(0x7ff0000000000000ull);             /* from +inf */
        for (int p = 0; p < K.N; p++) if (sh.red[p] < mints2) mints2 = sh.red[p];
        const int ok = is_normal(mints2);
        const double mim = ok ? mints2 : 1.0;
        double r = w_sqrt(mim, fl);
        r = w_mul(r, sh.dtdone, fl);
        const double dtnA = w_mul(r, S.s7, fl);
        double dtnew = ok ? dtnA : sh.dtnB;
        /* REBOUND's floor, min_dt +0 by default: the comparison is false */
        const double ad = w_abs(dtnew);
        const double dtnF = w_copysign(K.mindt, dtnew);
        if (ad < K.mindt) dtnew = dtnF;
        /* fabs(dt_new/dt_done) < safety_factor -> reject; larger than
         * 1/safety_factor -> clamp to dt_done/safety_factor */
        const double rr = w_div(dtnew, sh.dtdone, fl);
        const double ar = w_abs(rr);
        const int reject = ar < 0.25, larger = 1.0 < ar, past = K.invb < rr;
        if (!reject && larger && past) dtnew = sh.dtnB;
        sh.dtnew = dtnew;
        sh.accepted = !reject;
    });
}

/* ---- one step attempt: step_attempt() (src/ias15_cft.c:2885) ------------ */
template <class X>
__host__ __device__ inline void step_attempt(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    gravity(ex, K, S, sh, fl);                                  /* a at the step's start */
    ex.each([&](int p){
        fl.phase = PH_STEPSTART;
        for (int c = 0; c < 3; c++){
            S.x0[c][p] = S.x[c][p]; S.v0[c][p] = S.v[c][p]; S.a0[c][p] = S.a[c][p];
            for (int m = 0; m < 7; m++) S.csb[m][c][p] = 0.0;
            /* g from b */
            for (int m = 0; m < 6; m++){
                double t1 = w_mul(S.b[6][c][p], K.d[15 + m], fl);
                for (int j = 5; j > m; j--) t1 = w_add(t1, w_mul(S.b[j][c][p], K.d[j * (j - 1) / 2 + m], fl), fl);
                S.g[m][c][p] = w_add(t1, S.b[m][c][p], fl);
            }
            S.g[6][c][p] = S.b[6][c][p];
        }
    });
    ex.one([&]{
        fl.phase = PH_STEPSTART;
        for (int n = 1; n < 8; n++){
            sh.D1[n] = w_mul(S.dt, K.h[n], fl);                 /* fl(dt h_n) */
            sh.D2[n] = w_mul(sh.D1[n], K.khalf, fl);            /* and its exact half */
        }
        sh.spce = 1e300; sh.spcelast = 2.0; sh.iters = 0; sh.capped = 0;
        S.attempts++;
        S.last_dt = S.dt;
    });
    for (;;){
        ex.one([&]{
            int stop = 0;
            if (sh.spce < K.tol) stop = 1;
            else if (sh.iters > 2 && sh.spcelast <= sh.spce) stop = 1;
            else if (sh.iters >= K.max_iter){ S.max_exceeded++; sh.capped = 1; stop = 1; }
            if (!stop){ sh.spcelast = sh.spce; sh.spce = 0.0; sh.iters++; }
            sh.stop = stop;
        });
        if (sh.stop) break;
        for (int n = 1; n < 8; n++){
            /* predict_positions(n), the FMA form */
            ex.each([&](int p){
                fl.phase = PH_PREDICT;
                const double *hf = &K.hf[7 * (n - 1)];
                for (int c = 0; c < 3; c++){
                    double t1 = S.b[6][c][p];
                    for (int lvl = 0; lvl < 6; lvl++) t1 = w_fma(t1, hf[lvl], S.b[5 - lvl][c][p], fl);
                    t1 = w_fma(t1, hf[6], S.a0[c][p], fl);
                    t1 = w_fma(t1, sh.D2[n], S.v0[c][p], fl);
                    t1 = w_mul(t1, sh.D1[n], fl);
                    t1 = w_sub(t1, S.csx[c][p], fl);
                    S.x[c][p] = w_add(t1, S.x0[c][p], fl);
                }
            });
            gravity(ex, K, S, sh, fl);
            /* at = a; correct(n) */
            ex.each([&](int p){
                fl.phase = PH_CORRECT;
                double mat = 0.0, mtmp = 0.0;
                for (int c = 0; c < 3; c++){
                    const double at = S.a[c][p];
                    S.at[c][p] = at;
                    double gk = at, gkcs = 0.0;                 /* gk = at; gk_cs = gravity_cs, +0 */
                    const double na0 = w_neg(S.a0[c][p]);
                    double y = w_sub(na0, gkcs, fl), t = w_add(gk, y, fl), u = w_sub(t, gk, fl);
                    gkcs = w_sub(u, y, fl); gk = t;
                    y = w_sub(0.0, gkcs, fl); t = w_add(gk, y, fl); u = w_sub(t, gk, fl);   /* csa0 is +0 */
                    gkcs = w_sub(u, y, fl); gk = t;
                    const int base = n * (n - 1) / 2;
                    const double told = S.g[n - 1][c][p];
                    double gnew = w_mul(gk, K.rinv[base], fl);
                    for (int m = 1; m < n; m++){
                        gnew = w_sub(gnew, S.g[m - 1][c][p], fl);
                        gnew = w_mul(gnew, K.rinv[base + m], fl);
                    }
                    S.g[n - 1][c][p] = gnew;
                    const double tmp = w_sub(gnew, told, fl);
                    const int cbase = (n - 1) * (n - 2) / 2;
                    for (int m = 0; m < n - 1; m++){
                        const double prod = w_mul(tmp, K.c[cbase + m], fl);
                        add_cs(S.b[m][c][p], S.csb[m][c][p], prod, fl);
                    }
                    add_cs(S.b[n - 1][c][p], S.csb[n - 1][c][p], tmp, fl);
                    if (n == 7){
                        /* pc_error(): the maxima over the normal values, from +0 */
                        const double aa = w_abs(at), ta = w_abs(tmp);
                        if (is_normal(aa) && aa > mat) mat = aa;
                        if (is_normal(ta) && ta > mtmp) mtmp = ta;
                    }
                }
                if (n == 7){ sh.red[p] = mat; sh.red2[p] = mtmp; }
            });
            if (n == 7) ex.one([&]{
                fl.phase = PH_PCERROR;
                double maxak = 0.0, maxb6 = 0.0;
                for (int p = 0; p < K.N; p++){
                    if (sh.red[p] > maxak) maxak = sh.red[p];
                    if (sh.red2[p] > maxb6) maxb6 = sh.red2[p];
                }
                sh.spce = w_div(maxb6, maxak, fl);
            });
        }
    }
    ex.one([&]{
        S.pc_total += sh.iters;
        if (sh.iters > S.pc_max) S.pc_max = sh.iters;
        S.last_iters = sh.iters; S.last_capped = sh.capped;
        sh.dtdone = S.dt;                                       /* SDTDONE = SDT */
    });
    if (K.adaptive){
        choose_timestep(ex, K, S, sh, fl);
        ex.one([&]{ if (!sh.accepted) S.rejected++; S.dt = sh.dtnew; });
    }else{
        ex.one([&]{ sh.accepted = 1; });
    }
    /* A rejected step restores the step's start; the end-of-step update is
     * still computed, as the port computes it on every lane before it
     * restores the system, and only an accepted step keeps it. */
    ex.each([&](int p){
        fl.phase = PH_UPDATE;
        const double dt = sh.dtdone;
        const int acc = sh.accepted;
        for (int c = 0; c < 3; c++){
            if (!acc){ S.x[c][p] = S.x0[c][p]; S.v[c][p] = S.v0[c][p]; S.a[c][p] = S.a0[c][p]; }
            double x0 = S.x0[c][p], csx = S.csx[c][p], v0 = S.v0[c][p], csv = S.csv[c][p];
            const double v0s = v0, a0 = S.a0[c][p];
            for (int m = 6; m >= 0; m--){
                double t1 = w_mul(S.b[m][c][p], K.kposr[m], fl);
                t1 = w_mul(t1, dt, fl); t1 = w_mul(t1, dt, fl);
                add_cs(x0, csx, t1, fl);
            }
            double t1 = w_mul(a0, K.kposr[7], fl);
            t1 = w_mul(t1, dt, fl); t1 = w_mul(t1, dt, fl);
            add_cs(x0, csx, t1, fl);
            t1 = w_mul(v0s, dt, fl);
            add_cs(x0, csx, t1, fl);
            for (int m = 6; m >= 0; m--){
                t1 = w_mul(S.b[m][c][p], K.kvelr[m], fl);
                t1 = w_mul(t1, dt, fl);
                add_cs(v0, csv, t1, fl);
            }
            t1 = w_mul(a0, dt, fl);
            add_cs(v0, csv, t1, fl);
            if (acc){ S.x0[c][p] = x0; S.csx[c][p] = csx; S.v0[c][p] = v0; S.csv[c][p] = csv; }
        }
    });
    /* t += dt_done: the plain sum and the exact pair. A rejected attempt
     * adds +0, which is exact and raises nothing. */
    ex.one([&]{
        fl.phase = PH_CLOCK;
        const double dtadd = sh.accepted ? sh.dtdone : 0.0;
        S.tplain = w_add(S.tplain, dtadd, fl);
        double tr, terr;
        w_augadd(S.thi, dtadd, tr, terr, fl);
        S.tlo = w_add(S.tlo, terr, fl);
        S.thi = tr;
        S.last_accepted = sh.accepted;
        if (sh.accepted){ S.dtlast = sh.dtdone; S.done++; }
        /* the prediction's ratio: dt/dt_done after an accepted step,
         * dt/dt_last_done after a rejected one, none before any step was
         * accepted (s_eq(SDTLAST, K0): -0 counts as zero) */
        fl.phase = PH_NEXTSTEP;
        sh.skip = 1; sh.use_r = 0;
        double rdiv = 1.0;
        if (sh.accepted){ rdiv = sh.dtdone; sh.skip = 0; }
        else if (!(S.dtlast == 0.0)){ rdiv = S.dtlast; sh.skip = 0; sh.use_r = 1; }
        if (!sh.skip){
            sh.ratio = w_div(S.dt, rdiv, fl);
            const double q1 = sh.ratio;
            sh.q[1] = q1;
            sh.q[2] = w_mul(q1, q1, fl);
            sh.q[3] = w_mul(q1, sh.q[2], fl);
            sh.q[4] = w_mul(sh.q[2], sh.q[2], fl);
            sh.q[5] = w_mul(sh.q[2], sh.q[3], fl);
            sh.q[6] = w_mul(sh.q[3], sh.q[3], fl);
            sh.q[7] = w_mul(sh.q[3], sh.q[4], fl);
        }
    });
    ex.each([&](int p){
        fl.phase = PH_NEXTSTEP;
        const int acc = sh.accepted, skip = sh.skip, use_r = sh.use_r;
        const int zero = !skip && 20.0 < sh.ratio;
        for (int c = 0; c < 3; c++){
            S.x[c][p] = S.x0[c][p]; S.v[c][p] = S.v0[c][p];
            if (acc) for (int m = 0; m < 7; m++){ S.er[m][c][p] = S.e[m][c][p]; S.br[m][c][p] = S.b[m][c][p]; }
            if (skip) continue;
            double E_[7], B_[7], be[7], eo[7];
            for (int m = 0; m < 7; m++){
                E_[m] = use_r ? S.er[m][c][p] : S.e[m][c][p];
                B_[m] = use_r ? S.br[m][c][p] : S.b[m][c][p];
            }
            for (int m = 0; m < 7; m++) be[m] = w_sub(B_[m], E_[m], fl);
            for (int m = 0; m < 7; m++){
                /* e_m = q_{m+1} (b6 C(7,m+1) + ... + b_{m+1} C(m+2,m+1) + b_m), left to right */
                double s = 0.0;
                int first = 1;
                for (int j = 6; j > m; j--){
                    const double u = w_mul(B_[j], binom_d(j + 1, m + 1), fl);
                    if (first){ s = u; first = 0; } else s = w_add(s, u, fl);
                }
                s = first ? B_[m] : w_add(s, B_[m], fl);
                eo[m] = w_mul(sh.q[m + 1], s, fl);
            }
            for (int m = 0; m < 7; m++){
                const double bo = w_add(eo[m], be[m], fl);
                S.e[m][c][p] = zero ? 0.0 : eo[m];
                S.b[m][c][p] = zero ? 0.0 : bo;
            }
        }
    });
}

/* Each thread's flags into the system's: the union, and the first raise
 * of the lowest thread that raised one. */
template <class X>
__host__ __device__ inline void gather_flags(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl){
    ex.each([&](int p){ sh.fbits[p] = fl.bits; sh.ffirst[p] = fl.first; });
    ex.one([&]{
        for (int p = 0; p < K.N; p++){
            if (sh.fbits[p] && !S.flag_first) S.flag_first = sh.ffirst[p];
            S.flag_bits |= sh.fbits[p];
        }
        if (K.flag_abort && (S.flag_bits & K.abort_mask)) S.stopped = 1;
    });
}

/* A sample block: step attempts until `target` steps are accepted, or the
 * system stops on a flag; then the energy of the sample line. */
template <class X>
__host__ __device__ inline void run_block(X &ex, const Consts &K, Sys &S, Shared &sh, Fl &fl, long long target){
    while (S.done < target && !S.stopped){
        step_attempt(ex, K, S, sh, fl);
        gather_flags(ex, K, S, sh, fl);
    }
    if (!S.stopped){
        energy(ex, K, S, sh, fl);
        gather_flags(ex, K, S, sh, fl);
    }
}

} /* namespace wave */

#endif /* WAVE_STEP_H */
