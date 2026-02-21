/*
 * GUTOE GPU — 3D Schrödinger solver: kernels + host driver
 * Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
 *
 * Compile (NVIDIA sm_86 = RTX 3070 Ti / Ampere):
 *   nvcc -O3 -arch=sm_86 -Xcompiler -fPIC -c cuda_main.cu -o cuda_main.o
 *   ar rcs libschrodinger.a cuda_main.o
 *
 * For AMD ROCm (HIP compatibility — same source):
 *   hipcc -O3 -Xcompiler -fPIC -c cuda_main.cu -o cuda_main.o
 *   ar rcs libschrodinger.a cuda_main.o
 */

#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Complex double arithmetic ───────────────────────────────────────────────*/

typedef struct { double re; double im; } C64;

static __device__ __host__ __forceinline__
C64 cadd(C64 a, C64 b)   { return {a.re+b.re, a.im+b.im}; }
static __device__ __host__ __forceinline__
C64 csub(C64 a, C64 b)   { return {a.re-b.re, a.im-b.im}; }
static __device__ __host__ __forceinline__
C64 cmulr(C64 a, double r){ return {a.re*r, a.im*r}; }
static __device__ __host__ __forceinline__
double cnorm2(C64 a)      { return a.re*a.re + a.im*a.im; }
/* Re(a* b) = a.re*b.re + a.im*b.im */
static __device__ __forceinline__
double re_conj_dot(C64 a, C64 b){ return a.re*b.re + a.im*b.im; }

/* ── Hex+z neighbour function ────────────────────────────────────────────────
 *
 * Layout: flat_idx = (z * hex_rows + r) * hex_cols + c
 * Returns 8 neighbors: 6 intra-layer hex + 2 inter-layer z±1
 */
static __device__ __forceinline__
void hex_nbrs_3d(int site, int hex_rows, int hex_cols, int layers,
                  int* out)
{
    int lsz = hex_rows * hex_cols;
    int z   = site / lsz;
    int rem = site % lsz;
    int r   = rem / hex_cols;
    int c   = rem % hex_cols;

    int dr[6], dc[6];
    if (r % 2 == 0) {
        dr[0]=-1;dc[0]=0;  dr[1]=-1;dc[1]=1;
        dr[2]= 0;dc[2]=-1; dr[3]= 0;dc[3]=1;
        dr[4]= 1;dc[4]=0;  dr[5]= 1;dc[5]=1;
    } else {
        dr[0]=-1;dc[0]=-1; dr[1]=-1;dc[1]=0;
        dr[2]= 0;dc[2]=-1; dr[3]= 0;dc[3]=1;
        dr[4]= 1;dc[4]=-1; dr[5]= 1;dc[5]=0;
    }
    for (int i = 0; i < 6; i++) {
        int nr = ((r+dr[i])%hex_rows+hex_rows)%hex_rows;
        int nc = ((c+dc[i])%hex_cols+hex_cols)%hex_cols;
        out[i] = (z*hex_rows+nr)*hex_cols+nc;
    }
    int zp = (z+layers-1)%layers, zn = (z+1)%layers;
    out[6] = (zp*hex_rows+r)*hex_cols+c;
    out[7] = (zn*hex_rows+r)*hex_cols+c;
}

/* ── Initialisation kernels ──────────────────────────────────────────────────*/

__global__
void init_rho_kernel(double* rho, int n, int center)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;
    rho[i] = -1.0 / (double)n;          /* neutralising background */
    if (i == center) rho[i] += 1.0;      /* +1 point charge         */
}

/* Hex Cartesian coordinates (device) */
static __device__ __forceinline__
void hex_cart(int r, int c, double* x, double* y)
{
    *x = (double)c - 0.5*(double)(r&1);
    *y = (double)r * 0.866025403784438;   /* sqrt(3)/2 */
}

__global__
void init_gaussian_kernel(C64* psi, int n, int hex_rows, int hex_cols, int layers,
                           double cx, double cy, double cz, double sigma)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int lsz = hex_rows*hex_cols;
    int z = i/lsz, rem = i%lsz;
    int r = rem/hex_cols, c = rem%hex_cols;

    double px, py;
    hex_cart(r, c, &px, &py);
    double dz = (double)z - cz;
    double d2 = (px-cx)*(px-cx) + (py-cy)*(py-cy) + dz*dz;

    double v = exp(-d2 / (2.0*sigma*sigma));
    psi[i] = {v, 0.0};
}

/* ── Jacobi step kernel ──────────────────────────────────────────────────────*/

__global__
void jacobi_step_kernel(const double* phi, const double* rho,
                          double* phi_new,
                          int n, int hex_rows, int hex_cols, int layers)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int nbrs[8];
    hex_nbrs_3d(i, hex_rows, hex_cols, layers, nbrs);

    double s = 0.0;
    for (int j = 0; j < 8; j++) s += phi[nbrs[j]];
    phi_new[i] = (s + 8.0*rho[i]) / 8.0;
}

/* ── Fused Hamiltonian + imaginary time step ─────────────────────────────────
 *
 * Reads psi once, writes psi_new once. Avoids second pass.
 * psi_new[i] = psi[i] - dtau * H*psi[i]
 */
__global__
void ham_and_step_kernel(const C64* psi, C64* psi_new,
                           const double* phi,
                           double alpha_em, double charge, double dtau,
                           int n, int hex_rows, int hex_cols, int layers)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int nbrs[8];
    hex_nbrs_3d(i, hex_rows, hex_cols, layers, nbrs);

    C64 ns = {0.0, 0.0};
    for (int j = 0; j < 8; j++) ns = cadd(ns, psi[nbrs[j]]);

    C64 kinetic   = csub(psi[i], cmulr(ns, 0.125));     /* /8 */
    C64 potential = cmulr(psi[i], alpha_em*charge*phi[i]);
    C64 h         = cadd(kinetic, potential);

    psi_new[i] = csub(psi[i], cmulr(h, dtau));
}

/* ── Energy expectation value: Re<psi|H|psi> ─────────────────────────────────
 *
 * Returns (e_total, e_kin, e_pot) by parallel reduction.
 */
__global__
void energy_kernel(const C64* psi, const double* phi,
                    double alpha_em, double charge,
                    double* out_total, double* out_kin, double* out_pot,
                    int n, int hex_rows, int hex_cols, int layers)
{
    extern __shared__ double sdata[];   /* 3 * blockDim.x doubles */
    int tid = threadIdx.x;
    int i   = blockIdx.x*blockDim.x + tid;

    double et = 0.0, ek = 0.0, ep = 0.0;
    if (i < n) {
        int nbrs[8];
        hex_nbrs_3d(i, hex_rows, hex_cols, layers, nbrs);
        C64 ns = {0.0,0.0};
        for (int j = 0; j < 8; j++) ns = cadd(ns, psi[nbrs[j]]);
        C64 kinetic   = csub(psi[i], cmulr(ns, 0.125));
        C64 potential = cmulr(psi[i], alpha_em*charge*phi[i]);
        ek = re_conj_dot(psi[i], kinetic);
        ep = re_conj_dot(psi[i], potential);
        et = ek + ep;
    }

    /* Three separate shared arrays */
    double* sk = sdata;
    double* se = sdata + blockDim.x;
    double* sp = sdata + 2*blockDim.x;
    sk[tid]=ek; se[tid]=et; sp[tid]=ep;
    __syncthreads();

    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            sk[tid] += sk[tid+s];
            se[tid] += se[tid+s];
            sp[tid] += sp[tid+s];
        }
        __syncthreads();
    }
    if (tid == 0) {
        out_kin  [blockIdx.x] = sk[0];
        out_total[blockIdx.x] = se[0];
        out_pot  [blockIdx.x] = sp[0];
    }
}

/* ── Norm reduction + scale ──────────────────────────────────────────────────*/

__global__
void norm_sq_kernel(const C64* psi, double* partial, int n)
{
    extern __shared__ double sm[];
    int tid = threadIdx.x;
    int i   = blockIdx.x*blockDim.x + tid;
    sm[tid] = (i < n) ? cnorm2(psi[i]) : 0.0;
    __syncthreads();
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) sm[tid] += sm[tid+s];
        __syncthreads();
    }
    if (tid == 0) partial[blockIdx.x] = sm[0];
}

__global__
void scale_kernel(C64* psi, double scale, int n)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) psi[i] = cmulr(psi[i], scale);
}

/* ── Host driver ─────────────────────────────────────────────────────────────*/

#define BLOCK 256
#define CUDA_CHECK(x) do { \
    cudaError_t e = (x); \
    if (e != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

extern "C"
void run_hydrogen_cuda(
    double alpha_em,
    int hex_rows, int hex_cols, int layers,
    int n_jacobi, int n_iter, int renorm_every,
    double dtau,
    double* out_e_total, double* out_e_kin, double* out_e_pot)
{
    int n       = hex_rows * hex_cols * layers;
    int n_blk   = (n + BLOCK - 1) / BLOCK;
    int lsz     = hex_rows * hex_cols;
    int center  = (layers/2)*lsz + (hex_rows/2)*hex_cols + (hex_cols/2);

    /* Allocate device memory */
    double *d_rho, *d_phi, *d_phi_new, *d_partial_e, *d_partial_k, *d_partial_p;
    C64    *d_psi, *d_psi2;
    CUDA_CHECK(cudaMalloc(&d_rho,       n       * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_phi,       n       * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_phi_new,   n       * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_psi,       n       * sizeof(C64)));
    CUDA_CHECK(cudaMalloc(&d_psi2,      n       * sizeof(C64)));
    CUDA_CHECK(cudaMalloc(&d_partial_e, n_blk   * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_partial_k, n_blk   * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_partial_p, n_blk   * sizeof(double)));

    /* Helper: host buffer for reductions */
    double* h_partial = (double*)malloc(n_blk * sizeof(double));

    /* ── Phase 1: Jacobi-Poisson for Coulomb field ───────────────────────── */
    CUDA_CHECK(cudaMemset(d_phi, 0, n * sizeof(double)));
    init_rho_kernel<<<n_blk, BLOCK>>>(d_rho, n, center);
    CUDA_CHECK(cudaDeviceSynchronize());

    for (int iter = 0; iter < n_jacobi; iter++) {
        jacobi_step_kernel<<<n_blk, BLOCK>>>(d_phi, d_rho, d_phi_new,
                                               n, hex_rows, hex_cols, layers);
        /* swap phi and phi_new */
        double* tmp = d_phi; d_phi = d_phi_new; d_phi_new = tmp;
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* ── Phase 2: Gaussian initialisation ───────────────────────────────── */
    double cx, cy, cz;
    {
        int r = hex_rows/2, c = hex_cols/2;
        cx = (double)c - 0.5*(double)(r&1);
        cy = (double)r * 0.866025403784438;
        cz = (double)(layers/2);
    }
    double sigma = 1.0 / alpha_em;
    if (sigma < 1.0) sigma = 1.0;   /* at least 1 lattice spacing */

    init_gaussian_kernel<<<n_blk, BLOCK>>>(d_psi, n, hex_rows, hex_cols, layers,
                                             cx, cy, cz, sigma);
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Initial normalisation */
    norm_sq_kernel<<<n_blk, BLOCK, BLOCK*sizeof(double)>>>(d_psi, d_partial_e, n);
    CUDA_CHECK(cudaMemcpy(h_partial, d_partial_e, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
    { double s=0.0; for(int i=0;i<n_blk;i++) s+=h_partial[i];
      scale_kernel<<<n_blk,BLOCK>>>(d_psi, 1.0/sqrt(s), n); }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* ── Phase 3: Imaginary time evolution ───────────────────────────────── */
    for (int iter = 0; iter < n_iter; iter++) {
        /* Fused: psi2 = psi - dtau * H*psi */
        ham_and_step_kernel<<<n_blk, BLOCK>>>(d_psi, d_psi2, d_phi,
                                               alpha_em, -1.0, dtau,
                                               n, hex_rows, hex_cols, layers);
        /* swap */
        C64* tmp = d_psi; d_psi = d_psi2; d_psi2 = tmp;

        /* Renormalise periodically */
        if ((iter+1) % renorm_every == 0) {
            norm_sq_kernel<<<n_blk, BLOCK, BLOCK*sizeof(double)>>>(d_psi, d_partial_e, n);
            CUDA_CHECK(cudaMemcpy(h_partial, d_partial_e, n_blk*sizeof(double),
                                  cudaMemcpyDeviceToHost));
            double s=0.0; for(int i=0;i<n_blk;i++) s+=h_partial[i];
            scale_kernel<<<n_blk,BLOCK>>>(d_psi, 1.0/sqrt(s), n);
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* ── Phase 4: Compute E = <psi|H|psi> ───────────────────────────────── */
    energy_kernel<<<n_blk, BLOCK, 3*BLOCK*sizeof(double)>>>(
        d_psi, d_phi, alpha_em, -1.0,
        d_partial_e, d_partial_k, d_partial_p,
        n, hex_rows, hex_cols, layers);
    CUDA_CHECK(cudaDeviceSynchronize());

    double e_total=0.0, e_kin=0.0, e_pot=0.0;
    CUDA_CHECK(cudaMemcpy(h_partial, d_partial_e, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) e_total+=h_partial[i];
    CUDA_CHECK(cudaMemcpy(h_partial, d_partial_k, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) e_kin+=h_partial[i];
    CUDA_CHECK(cudaMemcpy(h_partial, d_partial_p, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) e_pot+=h_partial[i];

    *out_e_total = e_total;
    *out_e_kin   = e_kin;
    *out_e_pot   = e_pot;

    /* ── Cleanup ─────────────────────────────────────────────────────────── */
    cudaFree(d_rho); cudaFree(d_phi); cudaFree(d_phi_new);
    cudaFree(d_psi); cudaFree(d_psi2);
    cudaFree(d_partial_e); cudaFree(d_partial_k); cudaFree(d_partial_p);
    free(h_partial);
}
