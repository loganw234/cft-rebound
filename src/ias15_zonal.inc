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
 * corrector schedule with them (sched_key()). main() calls zonal_ready()
 * once the problem and its active line are read, so a line this term
 * cannot honour is refused before anything is printed or issued.
 *
 * THE ORDER IS REBOUNDx'S. Body k, mass M, reference radius R, acts on
 * every other body j of its system, at (dx, dy, dz) = x_j - x_k, through
 * the expressions of REBOUNDx's gravitational_harmonics.c, J2's function
 * and then J4's, each evaluated as C evaluates it, left to right:
 *
 *   r2 = dx*dx + dy*dy + dz*dz;  r = sqrt(r2);  costheta2 = dz*dz/r2;
 *   J2: prefac = 3.*J2*R*R/r2/r2/r/2.;           fac = 5.*costheta2 - 1.;
 *       a_j += G*M*prefac*fac*dx, ...*dy, G*M*prefac*(fac - 2.)*dz;
 *       a_k -= G*m_j*prefac*fac*dx, ...*dy, G*m_j*prefac*(fac - 2.)*dz;
 *   J4: prefac = 5.*J4*R*R*R*R/r2/r2/r2/r/8.;    fac = 63.*costheta2*costheta2
 *                                                     - 42.*costheta2 + 3.;
 *       z's factor is (fac + 12. - 28.*costheta2), and the rest as J2's.
 *
 * so (((G*M)*prefac)*fac)*dx, and a_j = (a_j + J2's) + J4's, and body k's
 * is a_k minus each J2 reaction in ascending body order, then each J4
 * reaction. Why this order:
 *   - there is no FMA in it, as there is none in gravity, so one order
 *     serves both --arith forms;
 *   - it is the expression REBOUND's users run, so a binary64 comparison
 *     with REBOUNDx can be operation for operation (ROUND4.md, decision 4);
 *   - it is the same expression at every format, so raising the format is
 *     a controlled experiment, as it is for gravity.
 * BELIEVED, not read: REBOUNDx is not pinned here (decision 4) and was not
 * read for this file; the expressions above are this file's recollection of
 * its source. The first operation-for-operation comparison confirms them or
 * corrects this paragraph, and the gates below hold the physics either way.
 *
 * Bits, not spelling. Each coefficient is signed in the header and added,
 * so `fac - 1.` is fac + (-1): IEEE 754 defines subtraction as addition of
 * the negated operand, so the bits are the same, zeros included. The
 * divides by 2. and 8. stay divides. 3*J2*R*R and 5*J4*R^4 are the same for
 * every body, so each is formed once a call, in that association, which is
 * the value REBOUNDx forms for each body. G*M is per system: an ensemble's
 * body k has its own mass in each.
 *
 * Test particles. Under a problem's active line (testparticle_type 0) a
 * test particle feels the term and gives no reaction: the reaction is
 * added only for the bodies that act on active ones, the same N_interact
 * rule as gravity's rows and energy_all()'s sums. So a test particle's mass
 * moves no active body here either. (REBOUNDx subtracts G*m_j*... for every
 * j; for a massless one that is a signed zero, which moves no bit.) Body k
 * must be active, and a zonal line naming a test particle is refused.
 *
 * The potential's order is this file's own, as REBOUNDx's was not read:
 * after energy_all()'s pair terms, epot plus each interacting body's J2
 * potential in ascending body order, then each one's J4 potential, with
 *   U2 = ((G*M)*m_j) * (J2*R*R/r2/r/2.) * (3.*costheta2 - 1.)
 *   U4 = ((G*M)*m_j) * (J4*R*R*R*R/r2/r2/r/8.)
 *                    * (35.*costheta2*costheta2 - 30.*costheta2 + 3.)
 * the J terms of -(G M m_j / r)[1 - J2 (R/r)^2 P2 - J4 (R/r)^4 P4], each
 * with P_n's coefficients over its denominator.
 *
 * The coefficients are src/ias15_constants.h's zonal[], which
 * tools/gen_constants.py derives from P2 and P4 and checks against
 * mpmath's derivative of the potential; each is a small integer, exact in
 * every format. A body costs one square root and ten divides a force
 * evaluation here, two of them by 2 and 8.
 *
 * --zonal-selftest evaluates the term and the potential at a fixed set of
 * points and traces every library operation, its operands and result as
 * bits, and every element moved between them; tools/check_zonal.py holds
 * the trace to cft-fp256's golden model and the results to the potential's
 * gradient in mpmath. Nothing there restates the order: it reads it from
 * the trace.
 */

#ifdef IAS15_CFT_LIBRARY
/* The library half reads no problem file, so no zonal line: never reached. */
static void zonal_accel(void){
    (void)zonal_k; (void)ZR; (void)ZJ2; (void)ZJ4;
}
#else
static void puthex(const V a_);   /* the records', defined below the engine */

/* ------------------------------------------------------------------ */
/* The trace (--zonal-selftest)                                         */
/* ------------------------------------------------------------------ */
static int zt_on;                 /* trace every operation and every move */
static unsigned long zt_ops, zt_els;
static V zt_ka, zt_kb, zt_one;    /* a call's operands as they were, and one element */
static size_t zt_cap;
static V ZK[IAS15_N_ZONAL];        /* the header's zonal[], broadcast NMAX wide */
static char zk_name[IAS15_N_ZONAL][8];

/* An operand's name in the trace: its own, or a coefficient's k<i>. */
static const char *zt_name(const char *n, const V p){
    if (n) return n;
    for (int c = 0; c < IAS15_N_ZONAL; c++) if (p == ZK[c]) return zk_name[c];
    return "?";
}

/* An element's bits, most significant byte first: the binary interchange
 * encoding as one integer. */
static void zt_bits(const V v){
    for (size_t b = ESZ; b-- > 0; ) printf("%02x", v[b]);
}

static void zt_leaf(const char *name, const V v, size_t n){
    if (!zt_on) return;
    for (size_t i = 0; i < n; i++){ printf("zin %s %zu ", name, i); zt_bits(E(v, i)); printf("\n"); }
}

/* One element moved, never computed: a gather, a scatter or a broadcast. */
static void zmove(const char *dn, V d, size_t di, const char *sn, const V s, size_t si){
    memcpy(E(d, di), E(s, si), ESZ);
    if (zt_on) printf("zcp %s %zu %s %zu\n", dn, di, sn, si);
}

/* One library call of the term, d = a op b over n elements, through the
 * program's own correctly rounded wrappers. Traced, the call is made as it
 * always is, and then each element is issued again on its own to read its
 * own flags; the two must give the same bits. */
static void zcall(char op, const char *dn, V d, const char *an, const V a, const char *bn, const V b, size_t n){
    if (zt_on){
        if (zt_cap < n){
            free(zt_ka); free(zt_kb); free(zt_one);
            zt_ka = valloc(n); zt_kb = valloc(n); zt_one = valloc(1); zt_cap = n;
        }
        memcpy(zt_ka, a, n * ESZ);
        if (b) memcpy(zt_kb, b, n * ESZ);
    }
    switch (op){
    case '+': vadd(d, a, b, n); break;
    case '-': vsub(d, a, b, n); break;
    case '*': vmul(d, a, b, n); break;
    case '/': vdiv(d, a, b, n); break;
    case 'r': vsqrt(d, a, n); break;
    default: die("zonal: no operation '%c'", op);
    }
    if (!zt_on) return;
    static const char *const nm[] = { "add", "sub", "mul", "div", "sqrt" };
    const int k = op == '+' ? 0 : op == '-' ? 1 : op == '*' ? 2 : op == '/' ? 3 : 4;
    zt_ops++;
    printf("zop %lu %s %s %s %s %zu\n", zt_ops, nm[k], dn, zt_name(an, a), b ? zt_name(bn, b) : "-", n);
    for (size_t i = 0; i < n; i++){
        uint32_t fl = 0;
        cft_status st;
        const V x1 = E(zt_ka, i), y1 = E(zt_kb, i);
        switch (op){
        case '+': st = cft_run(dev, CFT_ADD, F, CFT_RNE, x1, NULL, y1, zt_one, 1, &fl, NULL); break;
        case '-': st = cft_run(dev, CFT_SUB, F, CFT_RNE, x1, NULL, y1, zt_one, 1, &fl, NULL); break;
        case '*': st = cft_run(dev, CFT_MUL, F, CFT_RNE, x1, y1, NULL, zt_one, 1, &fl, NULL); break;
        case '/': st = cft_div(dev, F, CFT_RNE, x1, y1, zt_one, 1, &fl, NULL); break;
        default:  st = cft_sqrt(dev, F, CFT_RNE, x1, zt_one, 1, &fl, NULL); break;
        }
        if (st != CFT_OK) die("zonal trace: %s, one element: %s (%s)", nm[k], cft_strerror(st), cft_last_error());
        if (memcmp(zt_one, E(d, i), ESZ))
            die("zonal trace: operation %lu (%s) element %zu gives other bits issued alone than in its call", zt_ops, nm[k], i);
        printf("zel %lu %zu ", zt_ops, i); zt_bits(x1); printf(" ");
        if (b){ zt_bits(y1); printf(" "); }
        zt_bits(E(d, i)); printf(" %x\n", (unsigned)fl);
        zt_els++;
    }
}

/* ------------------------------------------------------------------ */
/* The term                                                             */
/* ------------------------------------------------------------------ */
static int zready;
static size_t zna, zni, zreact;    /* N_active, the bodies that act on active ones, those of them not k */
static size_t zcap_m, zcap_e;
/* E (N-1) elements each: system s's bodies other than k, ascending, at
 * s (N-1) + jj */
static V zxj, zyj, zzj, zxk, zyk, zzk, zmj, zmk, zdx, zdy, zdz, zdxx, zdyy, zdzz, zr2, zr, zc2, zgm, zgmj,
         zcstb, zpf, zf, zfz, zw, zt, ztf, ztx, zty, ztfz, ztz, zu, zuf, zrx, zry, zufz, zrz, zax, zay, zaz,
         zgmm, zpp, zpl, zpu;
static V *const zbuf_m[] = { &zxj, &zyj, &zzj, &zxk, &zyk, &zzk, &zmj, &zmk, &zdx, &zdy, &zdz, &zdxx, &zdyy,
         &zdzz, &zr2, &zr, &zc2, &zgm, &zgmj, &zcstb, &zpf, &zf, &zfz, &zw, &zt, &ztf, &ztx, &zty, &ztfz, &ztz,
         &zu, &zuf, &zrx, &zry, &zufz, &zrz, &zax, &zay, &zaz, &zgmm, &zpp, &zpl, &zpu };
static V zak, zrs, zus, zcst;       /* 3E, 3E, E and 1 */

static void zonal_ready(void){
    int nchecked = 0;
    if (zready) return;
    zready = 1;
    zna = (eng_active == (size_t)-1 || eng_active > N) ? N : eng_active;
    zni = eng_tptype ? N : zna;
    if (zonal_k >= zna)
        die("the zonal line's body %zu is a test particle (active %zu): J2 and J4 belong to one of the first "
            "N_active bodies, which act on the others", zonal_k, zna);
    for (int c = 0; c < IAS15_N_ZONAL; c++){
        ZK[c] = kderived(ias15_zonal_hex[FI][c], ias15_zonal_dec[c], "zonal", c, &nchecked);
        snprintf(zk_name[c], sizeof zk_name[c], "k%d", c);
    }
    zreact = 0;
    for (size_t j = 0; j < zni; j++) zreact += j != zonal_k;
    const size_t m = E * (N - 1);
    if (zcap_m < m){
        for (size_t i = 0; i < sizeof zbuf_m / sizeof zbuf_m[0]; i++){ free(*zbuf_m[i]); *zbuf_m[i] = valloc(m); }
        zcap_m = m;
    }
    if (zcap_e < E){ free(zak); free(zrs); free(zus); zak = valloc(3 * E); zrs = valloc(3 * E); zus = valloc(E); zcap_e = E; }
    if (!zcst) zcst = valloc(1);
}

/* The geometry both halves of the term read: the bodies gathered, then
 * dx, dy and dz, r2 = (dx*dx + dy*dy) + dz*dz, r, costheta2 and G*M. With
 * with_a, each body's acceleration is gathered too. */
static void zonal_geometry(int with_a){
    const size_t n1 = N - 1, m = E * n1;
    for (size_t s = 0; s < E; s++){
        const size_t pk = s * N + zonal_k;
        for (size_t jj = 0; jj < n1; jj++){
            const size_t e = s * n1 + jj, pj = s * N + (jj < zonal_k ? jj : jj + 1);
            zmove("xj", zxj, e, "x", x, 3 * pj); zmove("yj", zyj, e, "x", x, 3 * pj + 1); zmove("zj", zzj, e, "x", x, 3 * pj + 2);
            zmove("xk", zxk, e, "x", x, 3 * pk); zmove("yk", zyk, e, "x", x, 3 * pk + 1); zmove("zk", zzk, e, "x", x, 3 * pk + 2);
            zmove("mj", zmj, e, "mass", mass, pj); zmove("mk", zmk, e, "mass", mass, pk);
            if (with_a){ zmove("ax", zax, e, "a", a, 3 * pj); zmove("ay", zay, e, "a", a, 3 * pj + 1); zmove("az", zaz, e, "a", a, 3 * pj + 2); }
        }
    }
    zcall('-', "dx", zdx, "xj", zxj, "xk", zxk, m);          /* dx = p.x - source.x */
    zcall('-', "dy", zdy, "yj", zyj, "yk", zyk, m);
    zcall('-', "dz", zdz, "zj", zzj, "zk", zzk, m);
    zcall('*', "dxx", zdxx, "dx", zdx, "dx", zdx, m);
    zcall('*', "dyy", zdyy, "dy", zdy, "dy", zdy, m);
    zcall('*', "dzz", zdzz, "dz", zdz, "dz", zdz, m);
    zcall('+', "r2", zr2, "dxx", zdxx, "dyy", zdyy, m);       /* r2 = dx*dx + dy*dy + dz*dz */
    zcall('+', "r2", zr2, "r2", zr2, "dzz", zdzz, m);
    zcall('r', "r", zr, "r2", zr2, NULL, NULL, m);            /* r = sqrt(r2) */
    zcall('/', "c2", zc2, "dzz", zdzz, "r2", zr2, m);         /* costheta2 = dz*dz/r2 */
    zcall('*', "gm", zgm, "G", G, "mk", zmk, m);              /* G*source.m */
}

/* One pass's constant, (n+1)*Jn or Jn, times R n times, left to right, and
 * broadcast to every body: a call's one copy of what REBOUNDx forms for
 * each body. */
static void zonal_constant(const char *kn, const V k, const char *jn, const V jv, int nr){
    const size_t m = E * (N - 1);
    if (k) zcall('*', "cst", zcst, kn, k, jn, jv, 1);
    else   zmove("cst", zcst, 0, jn, jv, 0);
    for (int i = 0; i < nr; i++) zcall('*', "cst", zcst, "cst", zcst, "R", ZR, 1);
    for (size_t e = 0; e < m; e++) zmove("cstb", zcstb, e, "cst", zcst, 0);
}

/* A pass's term on every body and its reaction on body k, from the pass's
 * prefac (zpf), fac (zf) and z's factor (zfz):
 *   a_j += ((G*M)*prefac)*fac*dx, ...*dy, ((G*M)*prefac)*z's*dz
 *   a_k -= ((G*m_j)*prefac)*fac*dx, ...*dy, ((G*m_j)*prefac)*z's*dz,
 * the last in ascending body order, for the bodies that act on k only. */
static void zonal_apply(void){
    const size_t n1 = N - 1, m = E * n1;
    zcall('*', "t", zt, "gm", zgm, "pf", zpf, m);
    zcall('*', "tf", ztf, "t", zt, "f", zf, m);
    zcall('*', "tx", ztx, "tf", ztf, "dx", zdx, m);
    zcall('*', "ty", zty, "tf", ztf, "dy", zdy, m);
    zcall('*', "tfz", ztfz, "t", zt, "fz", zfz, m);
    zcall('*', "tz", ztz, "tfz", ztfz, "dz", zdz, m);
    zcall('+', "ax", zax, "ax", zax, "tx", ztx, m);
    zcall('+', "ay", zay, "ay", zay, "ty", zty, m);
    zcall('+', "az", zaz, "az", zaz, "tz", ztz, m);
    zcall('*', "u", zu, "gmj", zgmj, "pf", zpf, m);
    zcall('*', "uf", zuf, "u", zu, "f", zf, m);
    zcall('*', "rx", zrx, "uf", zuf, "dx", zdx, m);
    zcall('*', "ry", zry, "uf", zuf, "dy", zdy, m);
    zcall('*', "ufz", zufz, "u", zu, "fz", zfz, m);
    zcall('*', "rz", zrz, "ufz", zufz, "dz", zdz, m);
    for (size_t jj = 0; jj < n1; jj++){
        if ((jj < zonal_k ? jj : jj + 1) >= zni) continue;    /* a test particle gives no reaction */
        for (size_t s = 0; s < E; s++){
            zmove("rs", zrs, 3 * s, "rx", zrx, s * n1 + jj);
            zmove("rs", zrs, 3 * s + 1, "ry", zry, s * n1 + jj);
            zmove("rs", zrs, 3 * s + 2, "rz", zrz, s * n1 + jj);
        }
        zcall('-', "ak", zak, "ak", zak, "rs", zrs, 3 * E);
    }
}

/* Each body's acceleration from body zonal_k's J2 and J4, added to a, and
 * the reaction on body zonal_k. */
static void zonal_accel(void){
    zonal_ready();
    if (N < 2) return;
    const size_t n1 = N - 1, m = E * n1;
    for (size_t s = 0; s < E; s++)
        for (int c = 0; c < 3; c++) zmove("ak", zak, 3 * s + (size_t)c, "a", a, 3 * (s * N + zonal_k) + (size_t)c);
    zonal_geometry(1);
    zcall('*', "gmj", zgmj, "G", G, "mj", zmj, m);           /* G*p.m */
    /* J2: prefac = 3.*J2*R*R/r2/r2/r/2.; fac = 5.*costheta2 - 1.; z's fac - 2. */
    zonal_constant(NULL, ZK[IAS15_ZONAL_J2_NUM], "J2", ZJ2, 2);
    zcall('/', "pf", zpf, "cstb", zcstb, "r2", zr2, m);
    zcall('/', "pf", zpf, "pf", zpf, "r2", zr2, m);
    zcall('/', "pf", zpf, "pf", zpf, "r", zr, m);
    zcall('/', "pf", zpf, "pf", zpf, NULL, ZK[IAS15_ZONAL_J2_DEN], m);
    zcall('*', "f", zf, NULL, ZK[IAS15_ZONAL_J2_S1], "c2", zc2, m);
    zcall('+', "f", zf, "f", zf, NULL, ZK[IAS15_ZONAL_J2_S0], m);
    zcall('+', "fz", zfz, "f", zf, NULL, ZK[IAS15_ZONAL_J2_Z0], m);
    zonal_apply();
    /* J4: prefac = 5.*J4*R*R*R*R/r2/r2/r2/r/8.; fac = 63.*costheta2*costheta2
     * - 42.*costheta2 + 3.; z's fac + 12. - 28.*costheta2 */
    zonal_constant(NULL, ZK[IAS15_ZONAL_J4_NUM], "J4", ZJ4, 4);
    zcall('/', "pf", zpf, "cstb", zcstb, "r2", zr2, m);
    zcall('/', "pf", zpf, "pf", zpf, "r2", zr2, m);
    zcall('/', "pf", zpf, "pf", zpf, "r2", zr2, m);
    zcall('/', "pf", zpf, "pf", zpf, "r", zr, m);
    zcall('/', "pf", zpf, "pf", zpf, NULL, ZK[IAS15_ZONAL_J4_DEN], m);
    zcall('*', "f", zf, NULL, ZK[IAS15_ZONAL_J4_S2], "c2", zc2, m);
    zcall('*', "f", zf, "f", zf, "c2", zc2, m);
    zcall('*', "w", zw, NULL, ZK[IAS15_ZONAL_J4_S1], "c2", zc2, m);
    zcall('+', "f", zf, "f", zf, "w", zw, m);
    zcall('+', "f", zf, "f", zf, NULL, ZK[IAS15_ZONAL_J4_S0], m);
    zcall('+', "fz", zfz, "f", zf, NULL, ZK[IAS15_ZONAL_J4_Z0], m);
    zcall('*', "w", zw, NULL, ZK[IAS15_ZONAL_J4_Z1], "c2", zc2, m);
    zcall('+', "fz", zfz, "fz", zfz, "w", zw, m);
    zonal_apply();
    for (size_t s = 0; s < E; s++){
        for (size_t jj = 0; jj < n1; jj++){
            const size_t e = s * n1 + jj, pj = s * N + (jj < zonal_k ? jj : jj + 1);
            zmove("a", a, 3 * pj, "ax", zax, e); zmove("a", a, 3 * pj + 1, "ay", zay, e); zmove("a", a, 3 * pj + 2, "az", zaz, e);
        }
        for (int c = 0; c < 3; c++) zmove("a", a, 3 * (s * N + zonal_k) + (size_t)c, "ak", zak, 3 * s + (size_t)c);
    }
}

/* epot plus each interacting body's zpu, in ascending body order, every
 * system at once: element s of each call is the add system s would do. */
static void zonal_sum(V epot){
    const size_t n1 = N - 1;
    for (size_t jj = 0; jj < n1; jj++){
        if ((jj < zonal_k ? jj : jj + 1) >= zni) continue;    /* REBOUND's energy counts the N_interact bodies */
        for (size_t s = 0; s < E; s++) zmove("us", zus, s, "pu", zpu, s * n1 + jj);
        zcall('+', "epot", epot, "epot", epot, "us", zus, E);
    }
}

/* The J2 and J4 potential, into each system's epot, after the pair terms. */
static void zonal_potential(V epot){
    zonal_ready();
    if (N < 2) return;
    const size_t m = E * (N - 1);
    zonal_geometry(0);
    zcall('*', "gmm", zgmm, "gm", zgm, "mj", zmj, m);         /* (G*M)*m_j */
    /* U2 = ((G*M)*m_j) * (J2*R*R/r2/r/2.) * (3.*costheta2 - 1.) */
    zonal_constant(NULL, NULL, "J2", ZJ2, 2);
    zcall('/', "pp", zpp, "cstb", zcstb, "r2", zr2, m);
    zcall('/', "pp", zpp, "pp", zpp, "r", zr, m);
    zcall('/', "pp", zpp, "pp", zpp, NULL, ZK[IAS15_ZONAL_P2_DEN], m);
    zcall('*', "pl", zpl, NULL, ZK[IAS15_ZONAL_P2_S1], "c2", zc2, m);
    zcall('+', "pl", zpl, "pl", zpl, NULL, ZK[IAS15_ZONAL_P2_S0], m);
    zcall('*', "pu", zpu, "gmm", zgmm, "pp", zpp, m);
    zcall('*', "pu", zpu, "pu", zpu, "pl", zpl, m);
    zonal_sum(epot);
    /* U4 = ((G*M)*m_j) * (J4*R*R*R*R/r2/r2/r/8.)
     *      * (35.*costheta2*costheta2 - 30.*costheta2 + 3.) */
    zonal_constant(NULL, NULL, "J4", ZJ4, 4);
    zcall('/', "pp", zpp, "cstb", zcstb, "r2", zr2, m);
    zcall('/', "pp", zpp, "pp", zpp, "r2", zr2, m);
    zcall('/', "pp", zpp, "pp", zpp, "r", zr, m);
    zcall('/', "pp", zpp, "pp", zpp, NULL, ZK[IAS15_ZONAL_P4_DEN], m);
    zcall('*', "pl", zpl, NULL, ZK[IAS15_ZONAL_P4_S2], "c2", zc2, m);
    zcall('*', "pl", zpl, "pl", zpl, "c2", zc2, m);
    zcall('*', "w", zw, NULL, ZK[IAS15_ZONAL_P4_S1], "c2", zc2, m);
    zcall('+', "pl", zpl, "pl", zpl, "w", zw, m);
    zcall('+', "pl", zpl, "pl", zpl, NULL, ZK[IAS15_ZONAL_P4_S0], m);
    zcall('*', "pu", zpu, "gmm", zgmm, "pp", zpp, m);
    zcall('*', "pu", zpu, "pu", zpu, "pl", zpl, m);
    zonal_sum(epot);
}

/* The record's line naming the zonal term, with what it acts on: the
 * bodies that feel it and those whose reaction body k takes, counted from
 * the same rule the term runs. */
static void zonal_record(void){
    zonal_ready();
    printf("# zonal body %zu of every system, about z: R=", zonal_k); puthex(ZR);
    printf(" J2="); puthex(ZJ2); printf(" J4="); puthex(ZJ4);
    printf("; it acts on %zu bodies a system, and body %zu takes the reaction of %zu\n", N - 1, zonal_k, zreact);
}

/* ------------------------------------------------------------------ */
/* --zonal-selftest                                                     */
/* ------------------------------------------------------------------ */
/* The fixed points: two systems of eight, body 2 oblate in each, with a
 * different mass in each; the first six bodies active and the last two
 * test particles, one massless and one not. Every value is an exact
 * binary64, so the points are the same at every format. Off the plane and
 * on it (dz exactly +0: z_j is z_k's bits), near and far, both hemispheres,
 * two nearly over a pole, a massless active body, and in each system a
 * test particle with a mass, whose reaction must not be taken. */
static const char *const zst_g = "0x1.3000000000000p-1";
static const char *const zst_rj[3] = { "0x1.8000000000000p+0", "0x1.0ae79d91288f7p-6", "-0x1.ea5fb854f9344p-11" };
static const char *const zst_body[2][8][4] = {   /* m, x, y, z */
    {
    { "0x1.0000000000000p-10", "0x1.03bcba9082217p+1", "0x1.159330aa2ca54p-1", "0x1.ffbf45cd2a1a0p-6" },   /* 0: r = 1.32 R, on the plane */
    { "0x1.8000000000000p-6", "0x1.65be1798b882bp+0", "-0x1.1911ba16234efp+0", "0x1.924e23a428094p-1" },  /* 1: r = 1.18 R, u = +0.427 */
    { "0x1.8000000000000p+3", "0x1.401060393f809p-3", "-0x1.80102e85d6419p-4", "0x1.ffbf45cd2a1a0p-6" },  /* 2: the oblate body, M = 12 */
    { "0x1.4000000000000p-2", "-0x1.1d32105571f09p+3", "0x1.7702082b33fcdp+2", "-0x1.1f75c43b5f29ep+2" }, /* 3: r = 7.84 R, u = -0.385 */
    { "0x0.0p+0", "0x1.1e5ab8f60e75cp-2", "-0x1.088b2cfc11e5dp-5", "0x1.406ecb2830e08p+1" },               /* 4: r = 1.65 R, u = +0.998, massless */
    { "0x1.0000000000000p-4", "-0x1.5ccce974cc5e6p+4", "-0x1.a166db16cbb12p+3", "0x1.ffbf45cd2a1a0p-6" },  /* 5: r = 16.99 R, on the plane */
    { "0x0.0p+0", "-0x1.2b6ef94dacf8cp+1", "0x1.d4876caada31dp+1", "0x1.4402c4d248954p+0" },               /* 6: a test particle, u = +0.264 */
    { "0x1.0000000000000p-12", "0x1.2844f08f846b8p+2", "-0x1.e9b6e3da0aab3p-1", "-0x1.5b27f667dedffp+1" }, /* 7: a test particle with a mass, u = -0.516 */
    }, {
    { "0x1.0000000000000p-8", "-0x1.68c6b5cc4c178p+0", "-0x1.27f938f6b9722p+0", "-0x1.70b923f6875e4p-1" }, /* 0: r = 1.26 R, u = -0.336 */
    { "0x1.0000000000000p-5", "0x1.e75ac8739c906p-2", "-0x1.889a6ed77c3eep+0", "-0x1.5ff96e1e1e56bp-4" },  /* 1: r = 1.27 R, on the plane */
    { "0x1.4000000000000p+2", "-0x1.1ffb14ab2a322p-2", "0x1.c005e3de7b0eep-3", "-0x1.5ff96e1e1e56bp-4" },  /* 2: the oblate body, M = 5 */
    { "0x1.c000000000000p-3", "0x1.4f3711584c052p+2", "0x1.db5e63f4576afp+2", "0x1.77c9a3f539023p+1" },   /* 3: r = 6.38 R, u = +0.316 */
    { "0x1.0000000000000p-9", "-0x1.60793a8cdfa84p-2", "0x1.a0d807f5bec28p-2", "-0x1.8abc66a47af01p+1" },  /* 4: r = 2.00 R, u = -0.998 */
    { "0x1.2000000000000p-4", "0x1.f3c62c12bb1bap+0", "0x1.1c76832102aeep+1", "0x1.67b798ef3ce40p+0" },   /* 5: r = 2.23 R, u = +0.445 */
    { "0x0.0p+0", "-0x1.a16f8db9f45cbp+1", "-0x1.19653b79b6807p-2", "-0x1.5ff96e1e1e56bp-4" },             /* 6: a test particle, on the plane */
    { "0x1.0000000000000p-14", "0x1.37ca6965ae252p+0", "-0x1.e7a1a63172343p+1", "0x1.eb9aa838d4df5p+0" },  /* 7: a test particle with a mass, u = +0.423 */
    } };

static int zonal_selftest(void){
    const size_t ns = 2, nb = 8, k = 2, act = 6;
    E = ns; N = nb; NB = E * N; N3 = 3 * NB; L = 3 * N; PS = N * (N - 1) / 2; P = E * PS;
    NMAX = N3 > P ? N3 : P; if (NMAX < 8) NMAX = 8;
    mass = valloc(NB); x = valloc(N3); a = valloc(N3); G = valloc(NMAX);
    { V g = from_text(zst_g, 1); vbcast(G, g, NMAX); free(g); }
    ZR = from_text(zst_rj[0], 1); ZJ2 = from_text(zst_rj[1], 1); ZJ4 = from_text(zst_rj[2], 1);
    zonal_k = k; zonal_on = 1;
    eng_active = act; eng_tptype = 0;
    for (size_t s = 0; s < ns; s++)
        for (size_t j = 0; j < nb; j++){
            V v = from_text(zst_body[s][j][0], 1); memcpy(E(mass, s * nb + j), v, ESZ); free(v);
            for (int c = 0; c < 3; c++){ v = from_text(zst_body[s][j][1 + c], 1); memcpy(E(x, 3 * (s * nb + j) + (size_t)c), v, ESZ); free(v); }
        }
    vzero(a, N3);
    zonal_ready();
    V ep = valloc(E); vzero(ep, E);
    printf("zonal selftest format=%s E=%zu N=%zu k=%zu active=%zu tptype=%d esz=%zu\n",
           cft_format_name(F), E, N, zonal_k, zna, eng_tptype, ESZ);
    printf("zst G "); zt_bits(G); printf("\nzst R "); zt_bits(ZR); printf("\nzst J2 "); zt_bits(ZJ2);
    printf("\nzst J4 "); zt_bits(ZJ4); printf("\n");
    for (size_t s = 0; s < E; s++)
        for (size_t j = 0; j < N; j++){
            const size_t p = s * N + j;
            printf("zst body %zu %zu ", s, j); zt_bits(E(mass, p));
            for (int c = 0; c < 3; c++){ printf(" "); zt_bits(E(x, 3 * p + (size_t)c)); }
            printf("\n");
        }
    zt_on = 1;
    zt_leaf("x", x, N3); zt_leaf("mass", mass, NB); zt_leaf("a", a, N3); zt_leaf("G", G, NMAX);
    zt_leaf("R", ZR, 1); zt_leaf("J2", ZJ2, 1); zt_leaf("J4", ZJ4, 1); zt_leaf("epot", ep, E);
    for (int c = 0; c < IAS15_N_ZONAL; c++) zt_leaf(zk_name[c], ZK[c], NMAX);
    printf("zst begin accel\n");
    zonal_accel();
    printf("zst begin potential\n");
    zonal_potential(ep);
    zt_on = 0;
    for (size_t s = 0; s < E; s++)
        for (size_t j = 0; j < N; j++){
            printf("zst acc %zu %zu", s, j);
            for (int c = 0; c < 3; c++){ printf(" "); zt_bits(E(a, 3 * (s * N + j) + (size_t)c)); }
            printf("\n");
        }
    for (size_t s = 0; s < E; s++){ printf("zst pot %zu ", s); zt_bits(E(ep, s)); printf("\n"); }
    printf("zonal selftest: %lu operations traced, %lu elements, at %s; %zu bodies a system feel the term and %zu "
           "react on body %zu; flags seen 0x%x\n", zt_ops, zt_els, cft_format_name(F), N - 1, zreact, zonal_k,
           (unsigned)flags_union);
    free(ep);
    return 0;
}
#endif /* !IAS15_CFT_LIBRARY */
