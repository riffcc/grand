/*
 * GUTOE GPU — 3D Schrödinger solver kernels
 * Copyright (C) 2026 Riff Labs, AGPL-3.0-or-later
 *
 * Compiles with:
 *   NVIDIA: nvcc -O3 -arch=sm_80 -c schrodinger.cu -o schrodinger.o
 *   AMD:    hipcc -O3 -c schrodinger.cu -o schrodinger.o
 *
 * The same source compiles for both architectures via HIP compatibility.
 * Include as: #include <hip/hip_runtime.h> for ROCm, or <cuda_runtime.h> for CUDA.
 */

#ifdef __HIP_PLATFORM_HCC__
  #include <hip/hip_runtime.h>
#else
  #include <cuda_runtime.h>
#endif

#include <math.h>

/* ── Types ──────────────────────────────────────────────────────────────────── */

typedef struct { double re; double im; } Complex;

static __device__ __forceinline__ Complex cadd(Complex a, Complex b) {
    return {a.re + b.re, a.im + b.im};
}
static __device__ __forceinline__ Complex csub(Complex a, Complex b) {
    return {a.re - b.re, a.im - b.im};
}
static __device__ __forceinline__ Complex cmulr(Complex a, double r) {
    return {a.re * r, a.im * r};
}
static __device__ __forceinline__ double cnorm2(Complex a) {
    return a.re * a.re + a.im * a.im;
}

/* ── Neighbour index for the hex+z lattice ───────────────────────────────────
 *
 * Layout: flat_idx = (z * hex_rows + r) * hex_cols + c
 * Intra-layer (6 hex):  even row: (-1,0),(-1,1),(0,-1),(0,1),(1,0),(1,1)
 *                       odd  row: (-1,-1),(-1,0),(0,-1),(0,1),(1,-1),(1,0)
 * Inter-layer (2):      z-1, z+1 (same r,c)
 */

static __device__ int hex_nbr_count() { return 8; }

static __device__ void hex_nbrs(int site, int hex_rows, int hex_cols, int layers,
                                  int* out, int* n_out)
{
    int layer_sz = hex_rows * hex_cols;
    int z   = site / layer_sz;
    int rem = site % layer_sz;
    int r   = rem / hex_cols;
    int c   = rem % hex_cols;

    /* intra-layer hex offsets */
    int dr[6], dc[6];
    if (r % 2 == 0) {
        dr[0]=-1; dc[0]=0;  dr[1]=-1; dc[1]=1;
        dr[2]= 0; dc[2]=-1; dr[3]= 0; dc[3]=1;
        dr[4]= 1; dc[4]=0;  dr[5]= 1; dc[5]=1;
    } else {
        dr[0]=-1; dc[0]=-1; dr[1]=-1; dc[1]=0;
        dr[2]= 0; dc[2]=-1; dr[3]= 0; dc[3]=1;
        dr[4]= 1; dc[4]=-1; dr[5]= 1; dc[5]=0;
    }

    int k = 0;
    for (int i = 0; i < 6; i++) {
        int nr = ((r + dr[i]) % hex_rows + hex_rows) % hex_rows;
        int nc = ((c + dc[i]) % hex_cols + hex_cols) % hex_cols;
        out[k++] = (z * hex_rows + nr) * hex_cols + nc;
    }
    /* inter-layer z±1 */
    int z_prev = (z + layers - 1) % layers;
    int z_next = (z + 1)          % layers;
    out[k++] = (z_prev * hex_rows + r) * hex_cols + c;
    out[k++] = (z_next * hex_rows + r) * hex_cols + c;
    *n_out = k;  /* = 8 */
}

/* ── Jacobi-Poisson step kernel ─────────────────────────────────────────────
 *
 * phi_new[i] = (sum_{j in nbrs} phi[j]  +  k * rho[i]) / k
 */
extern "C" __global__
void jacobi_step_3d(
    const double* __restrict__ rho,
    const double* __restrict__ phi,
          double* __restrict__ phi_new,
    int n, int hex_rows, int hex_cols, int layers)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int nbrs[8]; int k;
    hex_nbrs(i, hex_rows, hex_cols, layers, nbrs, &k);

    double sum = 0.0;
    for (int j = 0; j < k; j++) sum += phi[nbrs[j]];
    phi_new[i] = (sum + (double)k * rho[i]) / (double)k;
}

/* ── Hamiltonian application kernel ─────────────────────────────────────────
 *
 * h_psi[i] = (psi[i] - mean_nbrs(psi[i]))   +   alpha * charge * phi[i] * psi[i]
 */
extern "C" __global__
void apply_hamiltonian_3d(
    const Complex* __restrict__ psi,
    const double*  __restrict__ phi,
          Complex* __restrict__ h_psi,
    double alpha_em, double charge,
    int n, int hex_rows, int hex_cols, int layers)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    int nbrs[8]; int k;
    hex_nbrs(i, hex_rows, hex_cols, layers, nbrs, &k);

    Complex nbr_sum = {0.0, 0.0};
    for (int j = 0; j < k; j++) nbr_sum = cadd(nbr_sum, psi[nbrs[j]]);

    Complex kinetic  = csub(psi[i], cmulr(nbr_sum, 1.0 / (double)k));
    Complex potential = cmulr(psi[i], alpha_em * charge * phi[i]);

    h_psi[i] = cadd(kinetic, potential);
}

/* ── Imaginary time step kernel ─────────────────────────────────────────────
 *
 * psi[i] -= dtau * h_psi[i]   (normalisation handled separately)
 */
extern "C" __global__
void imaginary_time_step(
          Complex* __restrict__ psi,
    const Complex* __restrict__ h_psi,
    double dtau, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    psi[i] = csub(psi[i], cmulr(h_psi[i], dtau));
}

/* ── Norm squared reduction ──────────────────────────────────────────────────
 *
 * Computes sum |psi[i]|^2 for normalisation.
 * Call with n/BLOCK_SIZE blocks, then sum the partial results on host.
 */
extern "C" __global__
void norm_sq_reduce(
    const Complex* __restrict__ psi,
          double*  __restrict__ partial_sums,
    int n)
{
    extern __shared__ double sdata[];
    int tid = threadIdx.x;
    int i   = blockIdx.x * blockDim.x + tid;

    sdata[tid] = (i < n) ? cnorm2(psi[i]) : 0.0;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) partial_sums[blockIdx.x] = sdata[0];
}

/* ── Scale kernel (for normalisation) ───────────────────────────────────────
 *
 * psi[i] /= norm
 */
extern "C" __global__
void scale_psi(Complex* __restrict__ psi, double inv_norm, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    psi[i] = cmulr(psi[i], inv_norm);
}
