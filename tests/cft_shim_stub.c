/* SPDX-License-Identifier: GPL-3.0-or-later - see cft_shim_stub.h */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cft_shim_stub.h"

static cft_device *dev;

static void need_dev(void){
    if (!dev){
        if (cft_open(NULL, 0, &dev) != CFT_OK){
            fprintf(stderr, "cft_open(software) failed: %s\n", cft_last_error());
            exit(2);
        }
    }
}

/* ---- elementwise helpers, d may alias an input ------------------- */
static void vop(cft_op op, int F, const void *a, const void *b, const void *c,
                void *d, size_t n){
    cft_status st = cft_run(dev, op, (cft_format)F, CFT_RNE, a, b, c, d, n, NULL, NULL);
    if (st != CFT_OK){ fprintf(stderr, "cft_run: %s\n", cft_strerror(st)); exit(2); }
}
static void vadd(int F, const void *a, const void *c, void *d, size_t n){ vop(CFT_ADD, F, a, NULL, c, d, n); }
static void vmul(int F, const void *a, const void *b, void *d, size_t n){ vop(CFT_MUL, F, a, b, NULL, d, n); }
static void vneg(int F, const void *a, void *d, size_t n){ vop(CFT_NEG, F, a, NULL, NULL, d, n); }
static void vfma(int F, const void *a, const void *b, const void *c, void *d, size_t n){ vop(CFT_FMA, F, a, b, c, d, n); }

/* one wide element holding the value i/j, broadcast n times */
static unsigned char *kbcast(int F, int i, int j, size_t n){
    size_t w = cft_format_size((cft_format)F);
    unsigned char one[32], two[32], q[32];
    int32_t vi = i, vj = j;
    if (cft_cvt_from_i32(dev, (cft_format)F, CFT_RNE, &vi, one, 1, NULL) != CFT_OK) exit(2);
    if (cft_cvt_from_i32(dev, (cft_format)F, CFT_RNE, &vj, two, 1, NULL) != CFT_OK) exit(2);
    if (cft_div(dev, (cft_format)F, CFT_RNE, one, two, q, 1, NULL, NULL) != CFT_OK) exit(2);
    unsigned char *v = malloc(n * w);
    for (size_t k = 0; k < n; k++) memcpy(v + k * w, q, w);
    return v;
}

/* ---- the integrator --------------------------------------------- */

static void *shim_create(void){
    struct cft_ias15_state *s = calloc(1, sizeof *s);
    s->format = CFT_FP64;
    s->E = 1;
    s->epsilon = 1e-9;
    s->min_dt = 0.0;
    s->adaptive_mode = 2;
    s->max_iter = 12;
    s->arith_fma = 0;
    snprintf(s->cft_abi, sizeof s->cft_abi, "%d.%d",
             CFT_ABI_VERSION_MAJOR, CFT_ABI_VERSION_MINOR);
    s->constants_digest = 0;
    return s;
}

static void shim_free(void *p){
    if (!p) return;
    cft_archive_state_free((struct cft_ias15_state*)p);
    free(p);
}

static unsigned char **blob_of(struct cft_ias15_state *s, int i){
    switch (i){
        case 0: return &s->x0;  case 1: return &s->v0;  case 2: return &s->a0;
        case 3: return &s->csx; case 4: return &s->csv; case 5: return &s->csa0;
        default: break;
    }
    i -= 6;
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

static void shim_step(struct reb_simulation *r, void *p){
    struct cft_ias15_state *s = (struct cft_ias15_state*)p;
    need_dev();
    size_t n = s->n_elem, F = (size_t)s->format, w = cft_format_size((cft_format)s->format);
    if (n == 0 || n != 3 * r->N){
        fprintf(stderr, "shim_step: state not set up (n_elem %zu, 3N %zu)\n", n, 3 * r->N);
        exit(2);
    }

    unsigned char *dt   = malloc(n * w);
    unsigned char *tmp  = malloc(n * w);
    {   /* dt as a wide vector, exactly */
        double d = r->dt;
        double *row = malloc(n * sizeof(double));
        for (size_t i = 0; i < n; i++) row[i] = d;
        cft_archive_promote_doubles(dev, s->format, row, dt, n);
        free(row);
    }
    unsigned char *half = kbcast((int)F, 1, 2, n);

    /* positions and velocities, with carries that depend on both */
    vmul((int)F, s->v0, dt, tmp, n);              /* tmp = v0*dt      */
    vadd((int)F, s->x0, tmp, s->x0, n);           /* x0 += tmp        */
    vfma((int)F, tmp, half, s->csx, s->csx, n);   /* csx += tmp/2     */
    vfma((int)F, s->a0, dt, s->v0, s->v0, n);     /* v0 += a0*dt      */
    vfma((int)F, s->v0, half, s->csv, s->csv, n); /* csv += v0/2      */
    vneg((int)F, s->x0, s->a0, n);                /* a0 = -x0         */
    vfma((int)F, s->a0, half, s->csa0, s->csa0, n);

    /* the seven levels, each feeding the next array */
    for (int j = 0; j < 7; j++){
        unsigned char *k = kbcast((int)F, 1, j + 3, n);
        vfma((int)F, s->b  [j], k, s->g  [j], s->g  [j], n);
        vfma((int)F, s->e  [j], k, s->b  [j], s->b  [j], n);
        vfma((int)F, s->br [j], k, s->e  [j], s->e  [j], n);
        vfma((int)F, s->er [j], k, s->br [j], s->br [j], n);
        vfma((int)F, s->csb[j], k, s->er [j], s->er [j], n);
        vfma((int)F, s->g  [j], k, s->csb[j], s->csb[j], n);
        free(k);
    }

    /* the binary64 view REBOUND sees */
    double *row = malloc(n * sizeof(double));
    cft_archive_demote_doubles(dev, s->format, s->x0, row, n);
    for (size_t i = 0; i < r->N; i++){
        r->particles[i].x = row[3*i]; r->particles[i].y = row[3*i+1]; r->particles[i].z = row[3*i+2];
    }
    cft_archive_demote_doubles(dev, s->format, s->v0, row, n);
    for (size_t i = 0; i < r->N; i++){
        r->particles[i].vx = row[3*i]; r->particles[i].vy = row[3*i+1]; r->particles[i].vz = row[3*i+2];
    }
    free(row);

    r->t += r->dt;
    free(dt); free(tmp); free(half);
}

static const struct reb_integrator cft_shim_integrator = {
    .documentation = "parcel B archive gate stand-in; not IAS15",
    .step = shim_step,
    .synchronize = NULL,
    .create = shim_create,
    .free = shim_free,
    .did_add_particle = NULL,
    .will_remove_particle = NULL,
    /* The registered default is the binary64 list. cft_archive_bind()
     * swaps in the list whose element_size matches the run's format. */
    .field_descriptor_list = NULL,
};

void cft_shim_register(void){
    static int done = 0;
    if (done) return;
    struct reb_integrator it = cft_shim_integrator;
    it.field_descriptor_list = cft_archive_descriptor_list(CFT_FP64);
    reb_integrator_register(it, CFT_IAS15_INTEGRATOR_NAME);
    done = 1;
}

int cft_shim_setup(struct reb_simulation *r, int format){
    need_dev();
    void *p = reb_simulation_set_integrator(r, CFT_IAS15_INTEGRATOR_NAME);
    if (!p) return -1;
    struct cft_ias15_state *s = (struct cft_ias15_state*)p;
    s->format = format;
    s->E = 1;
    size_t n = 3 * r->N;
    if (cft_archive_state_alloc(s, n)) return -1;
    if (cft_archive_bind(r)) return -1;

    /* Seed every blob with distinct non-zero values: x0/v0 from the
     * particles, the rest from a bounded index formula. */
    double *row = malloc(n * sizeof(double));
    for (size_t i = 0; i < r->N; i++){
        row[3*i] = r->particles[i].x; row[3*i+1] = r->particles[i].y; row[3*i+2] = r->particles[i].z;
    }
    cft_archive_promote_doubles(dev, format, row, s->x0, n);
    for (size_t i = 0; i < r->N; i++){
        row[3*i] = r->particles[i].vx; row[3*i+1] = r->particles[i].vy; row[3*i+2] = r->particles[i].vz;
    }
    cft_archive_promote_doubles(dev, format, row, s->v0, n);
    for (int bi = 2; bi < 48; bi++){
        for (size_t j = 0; j < n; j++)
            row[j] = (double)(((size_t)bi * 7 + j * 13) % 97 + 1) / 256.0
                     * (((bi + j) & 1) ? -1.0 : 1.0);
        cft_archive_promote_doubles(dev, format, row, *blob_of(s, bi), n);
    }
    free(row);
    return 0;
}

int cft_shim_blobs_nonzero(struct reb_simulation *r){
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    size_t w = cft_format_size((cft_format)s->format);
    int count = 0;
    for (int i = 0; i < 48; i++){
        unsigned char *b = *blob_of(s, i);
        for (size_t k = 0; k < s->n_elem * w; k++) if (b[k]){ count++; break; }
    }
    return count;
}

unsigned char *cft_shim_snapshot(struct reb_simulation *r, size_t *len){
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    size_t w = cft_format_size((cft_format)s->format);
    size_t blob = s->n_elem * w;
    size_t total = 48 * blob + r->N * sizeof(struct reb_particle) + sizeof(double) * 2
                 + sizeof(size_t) * 2 + sizeof(int) * 4;
    unsigned char *buf = malloc(total);
    size_t o = 0;
    for (int i = 0; i < 48; i++){ memcpy(buf + o, *blob_of(s, i), blob); o += blob; }
    for (size_t i = 0; i < r->N; i++){
        struct reb_particle p = r->particles[i];
        p.ap = NULL; p.sim = NULL; p.name = NULL;   /* pointers are not state */
        memcpy(buf + o, &p, sizeof p); o += sizeof p;
    }
    memcpy(buf + o, &r->t, sizeof(double));  o += sizeof(double);
    memcpy(buf + o, &r->dt, sizeof(double)); o += sizeof(double);
    memcpy(buf + o, &s->n_elem, sizeof(size_t)); o += sizeof(size_t);
    memcpy(buf + o, &s->E, sizeof(size_t)); o += sizeof(size_t);
    memcpy(buf + o, &s->format, sizeof(int)); o += sizeof(int);
    memcpy(buf + o, &s->max_iter, sizeof(int)); o += sizeof(int);
    memcpy(buf + o, &s->arith_fma, sizeof(int)); o += sizeof(int);
    memcpy(buf + o, &s->adaptive_mode, sizeof(int)); o += sizeof(int);
    *len = o;
    return buf;
}
