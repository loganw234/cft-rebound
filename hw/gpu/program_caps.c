/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * What this port's sequencer programs ask of a device, as the device's own
 * loader reads them (round 5, P1). hw/gpu/where.py builds it against the
 * tree's libcft.a with ias15_cft's own link line and runs it; it is not part
 * of any build target.
 *
 *   program_caps ARTIFACT [--scratch-in N] [--consts N] IMAGE.cftp...
 *
 * Opens ARTIFACT (an empty string is the software backend) and prints the
 * sequencer capacities it publishes. Then each image is loaded through
 * cft_program_load() and either cft_program_get_info() is printed - the
 * constants it addresses, the scratch slots its instructions touch and its
 * block preloads and reads back, its deposits and instructions - or the
 * loader's refusal, by its status and its sentence. Nothing is run.
 *
 * The boundary controls rewrite one header word of each image before it is
 * loaded, in this process's copy only:
 *   --scratch-in N  scratch_io's low half, the slots a lane is preloaded
 *                   with (only for an image whose flags carry
 *                   CFT_PROG_FLAG_SCRATCH_IO; the others are skipped);
 *   --consts N      n_consts, the constants the program addresses (only for
 *                   a CFT_PROG_FLAG_BANK_EXT image, which stores none, so
 *                   its length does not change; the others are skipped).
 * The device's own depth must load and one past it must be refused, so the
 * numbers this prints are the loader's, not this file's.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "cft.h"

static uint32_t rd32(const unsigned char *p){
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}
static void wr32(unsigned char *p, uint32_t v){
    p[0] = (unsigned char)v; p[1] = (unsigned char)(v >> 8); p[2] = (unsigned char)(v >> 16); p[3] = (unsigned char)(v >> 24);
}

int main(int argc, char **argv){
    if (argc < 2){
        fprintf(stderr, "usage: program_caps ARTIFACT [--scratch-in N] [--consts N] IMAGE.cftp...\n");
        return 2;
    }
    const char *art = argv[1][0] ? argv[1] : NULL;
    long scratch_in = -1, consts = -1;
    cft_device *dev = NULL;
    cft_status st = cft_open(art, 0, &dev);
    if (st != CFT_OK){
        printf("open %s: refused: %s (%s)\n", art ? art : "software", cft_strerror(st), cft_last_error());
        return 3;
    }
    cft_caps c;
    memset(&c, 0, sizeof c);
    c.struct_size = sizeof c;
    if (cft_get_caps(dev, &c) != CFT_OK){
        printf("cft_get_caps: %s\n", cft_last_error());
        return 3;
    }
    printf("device %s backend=%s max_scratch=%u max_consts=%u max_deposits=%u max_insns=%u seq_features=0x%x\n",
           art ? art : "software", c.backend, (unsigned)c.max_scratch, (unsigned)c.max_consts,
           (unsigned)c.max_deposits, (unsigned)c.max_insns, (unsigned)c.seq_features);
    int bad = 0;
    for (int i = 2; i < argc; i++){
        if (!strcmp(argv[i], "--scratch-in") && i + 1 < argc){ scratch_in = atol(argv[++i]); continue; }
        if (!strcmp(argv[i], "--consts") && i + 1 < argc){ consts = atol(argv[++i]); continue; }
        FILE *f = fopen(argv[i], "rb");
        if (!f){ printf("image %s: cannot be read\n", argv[i]); bad = 1; continue; }
        fseek(f, 0, SEEK_END);
        long n = ftell(f);
        fseek(f, 0, SEEK_SET);
        unsigned char *buf = malloc((size_t)n > 0 ? (size_t)n : 1);
        if (!buf || fread(buf, 1, (size_t)n, f) != (size_t)n){ printf("image %s: short read\n", argv[i]); fclose(f); bad = 1; free(buf); continue; }
        fclose(f);
        const char *name = strrchr(argv[i], '/');
        name = name ? name + 1 : argv[i];
        if (n < 32){ printf("image %s: %ld bytes, no header\n", name, n); free(buf); bad = 1; continue; }
        uint32_t flags = rd32(buf + 24), scr = rd32(buf + 28);
        if (scratch_in >= 0){
            if (!(flags & CFT_PROG_FLAG_SCRATCH_IO)){ free(buf); continue; }
            wr32(buf + 28, (scr & 0xFFFF0000u) | ((uint32_t)scratch_in & 0xFFFFu));
        }
        if (consts >= 0){
            if (!(flags & CFT_PROG_FLAG_BANK_EXT)){ free(buf); continue; }
            wr32(buf + 12, (uint32_t)consts);
        }
        cft_program *p = NULL;
        st = cft_program_load(dev, buf, (size_t)n, &p);
        if (st != CFT_OK){
            printf("image %s: refused: %s (%s)\n", name, cft_strerror(st), cft_last_error());
            free(buf);
            continue;
        }
        cft_program_info info;
        memset(&info, 0, sizeof info);
        info.struct_size = sizeof info;
        if (cft_program_get_info(p, &info) != CFT_OK){
            printf("image %s: cft_program_get_info: %s\n", name, cft_last_error());
            bad = 1;
        }else
            printf("image %s: loaded n_insns=%u n_consts=%u max_deposits=%u flags=0x%x n_scratch_in=%u "
                   "n_scratch_out=%u scratch_used=%u\n", name, (unsigned)info.n_insns, (unsigned)info.n_consts,
                   (unsigned)info.max_deposits, (unsigned)info.flags, (unsigned)info.n_scratch_in,
                   (unsigned)info.n_scratch_out, (unsigned)info.scratch_used);
        cft_program_free(p);
        free(buf);
    }
    cft_close(dev);
    return bad;
}
