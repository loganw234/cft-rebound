/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * THE KERNEL (ROUND5.md, P2): src/wave/wave_step.h's step in a HIP kernel,
 * one system to a workgroup, one thread a body. Included by wave_main.cpp
 * when it is built by hipcc (`make HIP=1 wave`, -DWAVE_HIP): hipcc for an
 * AMD GPU, or HIP's NVIDIA layer, which is nvcc.
 *
 * The executor: each() runs a body on its own thread and then waits at a
 * __syncthreads(); one() runs the system's scalar work on thread 0 and
 * waits. Every reduction - the corrector's maxima, PRS23's minimum, J2
 * and J4's reaction on the oblate body, the energy's sums - goes through
 * the workgroup's __shared__ block and is folded by thread 0 in the port's
 * order, never by a wave intrinsic, so the kernel is the same code on
 * NVIDIA's 32-wide warps and on RDNA3's wave32 and wave64.
 *
 * A launch is a sample block: every system steps until it has `target`
 * accepted steps (or stops on a flag), then forms its energy and writes
 * its sample. The state stays on the device between launches.
 */

struct DevExec {
    int N;
    template <class F> __device__ void each(F f){
        const int p = (int)threadIdx.x;
        if (p < N) f(p);
        __syncthreads();
    }
    template <class F> __device__ void one(F f){
        if (threadIdx.x == 0) f();
        __syncthreads();
    }
};

__constant__ Consts dK;

__global__ void wave_kernel(Sys *S, SampleOut *out, long long target, Attempt *log, int logcap, int *logn){
    __shared__ Shared sh;
    DevExec ex; ex.N = dK.N;
    Sys &y = S[blockIdx.x];
    Fl fl; fl.clear();
    while (y.done < target && !y.stopped){
        step_attempt(ex, dK, y, sh, fl);
        gather_flags(ex, dK, y, sh, fl);
        if (log && threadIdx.x == 0){
            const int i = logn[blockIdx.x]++;
            if (i < logcap){
                Attempt &a = log[(size_t)blockIdx.x * logcap + i];
                a.dt = y.last_dt; a.accepted = y.last_accepted; a.iters = y.last_iters; a.capped = y.last_capped;
            }
        }
        __syncthreads();
    }
    if (!y.stopped){
        energy(ex, dK, y, sh, fl);
        gather_flags(ex, dK, y, sh, fl);
    }
    if (threadIdx.x == 0) fill_sample(y, sh, dK.N, out[blockIdx.x]);
}

/* --zonal-selftest-from --device gpu: the term at the port's points, a
 * workgroup a system */
__global__ void zonal_selftest_kernel(Sys *S, double *pot, unsigned *flags){
    __shared__ Shared sh;
    DevExec ex; ex.N = dK.N;
    Fl fl; fl.clear();
    zonal_selftest_one(ex, dK, S[blockIdx.x], sh, fl, &pot[blockIdx.x]);
    gather_flags(ex, dK, S[blockIdx.x], sh, fl);
    if (threadIdx.x == 0) flags[blockIdx.x] = S[blockIdx.x].flag_bits;
}

#define HIPCHECK(x) do { hipError_t e_ = (x); if (e_ != hipSuccess) die(2, "%s: %s", #x, hipGetErrorString(e_)); } while (0)

/* The device, named: the one line every GPU record carries */
static std::string gpu_open(std::string &backend){
    int n = 0, dv = 0, rt = 0, drv = 0;
    HIPCHECK(hipGetDeviceCount(&n));
    if (n < 1) die(2, "--device gpu: no GPU is visible");
    HIPCHECK(hipSetDevice(0));
    HIPCHECK(hipGetDevice(&dv));
    hipDeviceProp_t pr;
    HIPCHECK(hipGetDeviceProperties(&pr, dv));
    HIPCHECK(hipRuntimeGetVersion(&rt));
    HIPCHECK(hipDriverGetVersion(&drv));
    char line[640];
#if defined(__HIP_PLATFORM_AMD__)
    const char *arch = pr.gcnArchName, *plat = "amd";
#else
    char archb[32]; snprintf(archb, sizeof archb, "sm_%d%d", pr.major, pr.minor);
    const char *arch = archb, *plat = "nvidia";
#endif
    snprintf(line, sizeof line, "hip:%d = %s (%s, %d multiprocessors, bus %04x:%02x:%02x.0, %zu MiB, warp %d), HIP platform %s, "
             "runtime %d, driver %d; of %d visible; a system a workgroup of N threads, the step in src/wave/",
             dv, pr.name, arch, pr.multiProcessorCount, pr.pciDomainID, pr.pciBusID, pr.pciDeviceID,
             (size_t)(pr.totalGlobalMem >> 20), pr.warpSize, plat, rt, drv, n);
    std::string b = "hip:" + std::to_string(dv) + ":";
    for (const char *c = arch; *c && *c != ':'; c++) b += *c;
    backend = b;
    return std::string(line);
}

struct GpuRunner {
    const Consts &K;
    std::vector<Sys> &S;
    std::vector<std::vector<Attempt>> *log;
    Sys *dS = nullptr;
    SampleOut *dOut = nullptr;
    Attempt *dLog = nullptr;
    int *dLogN = nullptr;
    int logcap = 0;
    GpuRunner(const Consts &K_, std::vector<Sys> &S_, std::vector<std::vector<Attempt>> *log_) : K(K_), S(S_), log(log_){
        const size_t E = S.size();
        HIPCHECK(hipMemcpyToSymbol(HIP_SYMBOL(dK), &K, sizeof K));
        HIPCHECK(hipMalloc((void **)&dS, E * sizeof(Sys)));
        HIPCHECK(hipMemcpy(dS, S.data(), E * sizeof(Sys), hipMemcpyHostToDevice));
        HIPCHECK(hipMalloc((void **)&dOut, E * sizeof(SampleOut)));
    }
    ~GpuRunner(){ hipFree(dS); hipFree(dOut); if (dLog) hipFree(dLog); if (dLogN) hipFree(dLogN); }
    void block(long long target, std::vector<SampleOut> &out){
        const size_t E = S.size();
        if (log){
            /* each block's attempts: its steps plus its rejections, bounded here */
            const long long have = out.empty() ? 0 : out[0].done;
            const int cap = (int)(4 * (target - have) + 64);
            if (cap > logcap){
                if (dLog) hipFree(dLog);
                if (dLogN) hipFree(dLogN);
                HIPCHECK(hipMalloc((void **)&dLog, E * (size_t)cap * sizeof(Attempt)));
                HIPCHECK(hipMalloc((void **)&dLogN, E * sizeof(int)));
                logcap = cap;
            }
            HIPCHECK(hipMemset(dLogN, 0, E * sizeof(int)));
        }
        const unsigned threads = (unsigned)K.N;
        hipLaunchKernelGGL(wave_kernel, dim3((unsigned)E), dim3(threads), 0, 0, dS, dOut, target, log ? dLog : nullptr,
                           logcap, log ? dLogN : nullptr);
        HIPCHECK(hipGetLastError());
        HIPCHECK(hipDeviceSynchronize());
        HIPCHECK(hipMemcpy(out.data(), dOut, E * sizeof(SampleOut), hipMemcpyDeviceToHost));
        if (log){
            std::vector<int> n(E);
            std::vector<Attempt> buf(E * (size_t)logcap);
            HIPCHECK(hipMemcpy(n.data(), dLogN, E * sizeof(int), hipMemcpyDeviceToHost));
            HIPCHECK(hipMemcpy(buf.data(), dLog, buf.size() * sizeof(Attempt), hipMemcpyDeviceToHost));
            for (size_t s = 0; s < E; s++){
                if (n[s] > logcap) die(3, "system %zu made %d step attempts in one sample block, past the log's %d", s, n[s], logcap);
                for (int i = 0; i < n[s]; i++) (*log)[s].push_back(buf[s * (size_t)logcap + i]);
            }
        }
    }
};

static unsigned gpu_zonal_selftest(const Consts &K, std::vector<Sys> &S, double *pot){
    const size_t E = S.size();
    Sys *dS; double *dp; unsigned *df;
    HIPCHECK(hipMemcpyToSymbol(HIP_SYMBOL(dK), &K, sizeof K));
    HIPCHECK(hipMalloc((void **)&dS, E * sizeof(Sys)));
    HIPCHECK(hipMalloc((void **)&dp, E * sizeof(double)));
    HIPCHECK(hipMalloc((void **)&df, E * sizeof(unsigned)));
    HIPCHECK(hipMemcpy(dS, S.data(), E * sizeof(Sys), hipMemcpyHostToDevice));
    hipLaunchKernelGGL(zonal_selftest_kernel, dim3((unsigned)E), dim3((unsigned)K.N), 0, 0, dS, dp, df);
    HIPCHECK(hipGetLastError());
    HIPCHECK(hipDeviceSynchronize());
    std::vector<unsigned> f(E);
    HIPCHECK(hipMemcpy(S.data(), dS, E * sizeof(Sys), hipMemcpyDeviceToHost));
    HIPCHECK(hipMemcpy(pot, dp, E * sizeof(double), hipMemcpyDeviceToHost));
    HIPCHECK(hipMemcpy(f.data(), df, E * sizeof(unsigned), hipMemcpyDeviceToHost));
    hipFree(dS); hipFree(dp); hipFree(df);
    unsigned u = 0;
    for (unsigned x : f) u |= x;
    return u;
}
