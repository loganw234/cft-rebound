/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * THE WAVE KERNEL'S ONE SOURCE, ON BOTH SIDES (ROUND5.md, P2).
 *
 * Every function the step runs is written once, as __host__ __device__
 * over double, in src/wave/. Under HIP (hipcc for an AMD GPU, or HIP's
 * NVIDIA layer, which is nvcc) the attributes are the compiler's own and
 * the same functions run in the kernel. Without HIP this header defines
 * both attributes away, so the same text compiles as plain C++ for the
 * host twin (`make wave-host`, MinGW g++ on the desktop and g++ on Linux).
 *
 * WHAT THE COMPILER MUST NOT DO, on either side. Each is a silent wrong
 * answer, and each has a flag in the Makefile and a check here:
 *   - contract a*b + c into an fma: -ffp-contract=off everywhere, and
 *     --fmad=false under nvcc. The only fma in the step is wave_fma(),
 *     written where the port's FMA form issues one (src/ias15_cft.c,
 *     predict_positions(), arith_fma);
 *   - relax IEEE semantics: no -ffast-math, no -fno-signed-zeros, no
 *     reassociation; under nvcc -prec-div=true -prec-sqrt=true;
 *   - flush subnormals: -fno-gpu-flush-denormals-to-zero under HIP-Clang,
 *     -ftz=false under nvcc;
 *   - call a library fma on the host: the host twin needs the CPU's own
 *     fused multiply-add (-mfma on x86-64), so that fma() is one
 *     correctly rounded operation and not whatever the C runtime ships.
 *     The #error below refuses a host build without it.
 * Gate 4 (src/wave/wave_corpus.cpp) MEASURES the arithmetic these flags
 * leave against libcft's software binary64, on the host and on each GPU;
 * nothing rests on the belief.
 */
#ifndef WAVE_HD_H
#define WAVE_HD_H

#if defined(__HIPCC__) || defined(__CUDACC__)
#  define WAVE_GPU_COMPILER 1
#  include <hip/hip_runtime.h>
#else
#  define WAVE_GPU_COMPILER 0
#  ifndef __host__
#    define __host__
#  endif
#  ifndef __device__
#    define __device__
#  endif
#endif

/* The device pass of a GPU compile (where there is no host CPU to ask). */
#if defined(__HIP_DEVICE_COMPILE__) || defined(__CUDA_ARCH__)
#  define WAVE_DEVICE_PASS 1
#else
#  define WAVE_DEVICE_PASS 0
#endif

#if !WAVE_DEVICE_PASS && (defined(__x86_64__) || defined(_M_X64)) && !defined(__FMA__)
#  error "src/wave on an x86-64 host needs -mfma: wave_fma() must be the CPU's single correctly rounded fused multiply-add, never a library routine"
#endif

#if defined(__FAST_MATH__)
#  error "src/wave must not be compiled with -ffast-math: every operation is the IEEE one, signed zeros and subnormals included"
#endif

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>

#if WAVE_GPU_COMPILER
#  define WAVE_INLINE __forceinline__
#else
#  define WAVE_INLINE inline
#endif

#ifndef WAVE_MAXN
#  define WAVE_MAXN 64          /* bodies a system: one 64-thread workgroup, GPU-E's shape a2 */
#endif

#endif /* WAVE_HD_H */
