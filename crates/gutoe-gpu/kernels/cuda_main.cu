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

/* ── Portability: compile with nvcc (CUDA) or hipcc (ROCm/HIP) ───────────── */
#ifdef __HIP_PLATFORM_AMD__
#  include <hip/hip_runtime.h>
   /* Map CUDA runtime API names to HIP equivalents */
#  define cudaError_t            hipError_t
#  define cudaSuccess            hipSuccess
#  define cudaGetErrorString     hipGetErrorString
#  define cudaMalloc             hipMalloc
#  define cudaFree               hipFree
#  define cudaMemcpy             hipMemcpy
#  define cudaMemset             hipMemset
#  define cudaMemcpyDeviceToHost hipMemcpyDeviceToHost
#  define cudaDeviceSynchronize  hipDeviceSynchronize
#  define cudaMemGetInfo         hipMemGetInfo
#else
#  include <cuda_runtime.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BLOCK 256   /* threads per block — used by all kernels */
#define CUDA_CHECK(x) do { \
    cudaError_t _e = (x); \
    if (_e != cudaSuccess) { \
        fprintf(stderr, "CUDA error %s at %s:%d\n", cudaGetErrorString(_e), __FILE__, __LINE__); \
        exit(1); \
    } \
} while(0)

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

/* Hex Cartesian coordinates (device + host) */
static __device__ __host__ __forceinline__
void hex_cart(int r, int c, double* x, double* y)
{
    *x = (double)c - 0.5*(double)(r&1);
    *y = (double)r * 0.866025403784438;   /* sqrt(3)/2 */
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

/* Coulomb initialisation: φ = α/r (exact continuum solution).
 * Starting from the analytic solution instead of φ=0 means SOR only needs
 * to correct the lattice discretization error — O(a²/r²) residual — which
 * converges in ~50-100 sweeps instead of ~1400 from cold start. */
__global__
void init_coulomb_kernel(double* phi, int n,
                          int hex_rows, int hex_cols, int layers,
                          double cx, double cy, double cz, double alpha_em)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int lsz = hex_rows * hex_cols;
    int z = i / lsz, rem = i % lsz;
    int r = rem / hex_cols, c = rem % hex_cols;

    double px, py;
    hex_cart(r, c, &px, &py);
    double dz = (double)z - cz;
    double d  = sqrt((px-cx)*(px-cx) + (py-cy)*(py-cy) + dz*dz);
    phi[i] = (d > 0.5) ? alpha_em / d : 2.0 * alpha_em;   /* same regularisation as obc_phi */
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

/* ── Open-boundary-condition (OBC) helpers ───────────────────────────────────
 *
 * For each of the 8 neighbours of site (z,r,c) we return the flat index if
 * in-bounds, or -1 if outside the box.  Out-of-bounds sites use:
 *   ψ  → 0          (Dirichlet: wavefunction vanishes at wall)
 *   φ  → α / r      (Coulomb potential from proton at box centre)
 */

static __device__ __forceinline__
double obc_phi(int nr, int nc, int nz,
               int hex_rows, int hex_cols, int layers,
               double cx, double cy, double cz, double alpha_em)
{
    double px, py;
    hex_cart(nr, nc, &px, &py);
    double dz = (double)nz - cz;
    double d  = sqrt((px-cx)*(px-cx) + (py-cy)*(py-cy) + dz*dz);
    return (d > 0.5) ? alpha_em / d : 2.0 * alpha_em;   /* regularise r→0 */
}

static __device__ __forceinline__
void hex_nbrs_obc(int site, int hex_rows, int hex_cols, int layers,
                   int out_idx[8], int out_nr[8], int out_nc[8], int out_nz[8])
{
    int lsz = hex_rows * hex_cols;
    int z   = site / lsz, rem = site % lsz;
    int r   = rem / hex_cols, c = rem % hex_cols;

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
        int nr = r + dr[i], nc = c + dc[i];
        out_nr[i] = nr; out_nc[i] = nc; out_nz[i] = z;
        if (nr < 0 || nr >= hex_rows || nc < 0 || nc >= hex_cols)
            out_idx[i] = -1;
        else
            out_idx[i] = (z*hex_rows + nr)*hex_cols + nc;
    }
    /* z ± 1 */
    int zp = z - 1, zn = z + 1;
    out_nr[6]=r; out_nc[6]=c; out_nz[6]=zp;
    out_idx[6] = (zp < 0)      ? -1 : (zp*hex_rows+r)*hex_cols+c;
    out_nr[7]=r; out_nc[7]=c; out_nz[7]=zn;
    out_idx[7] = (zn >= layers) ? -1 : (zn*hex_rows+r)*hex_cols+c;
}

/* OBC Jacobi — boundary sites contribute their Coulomb value */
__global__
void jacobi_obc_kernel(const double* phi, const double* rho, double* phi_new,
                        int n, int hex_rows, int hex_cols, int layers,
                        double cx, double cy, double cz, double alpha_em)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int idx[8], nr[8], nc[8], nz[8];
    hex_nbrs_obc(i, hex_rows, hex_cols, layers, idx, nr, nc, nz);

    double s = 0.0;
    for (int j = 0; j < 8; j++) {
        s += (idx[j] >= 0) ? phi[idx[j]]
                            : obc_phi(nr[j], nc[j], nz[j],
                                      hex_rows, hex_cols, layers,
                                      cx, cy, cz, alpha_em);
    }
    phi_new[i] = (s + 8.0*rho[i]) / 8.0;
}

/* ── 6-color SOR OBC kernel (in-place) ───────────────────────────────────────
 *
 * Proper multi-color Gauss-Seidel for the hex+z lattice.
 *
 * The hex (triangular) lattice is NOT bipartite — (z+r+c)%2 coloring has
 * same-color neighbours, which causes divergence with ω > 1.
 *
 * Valid 3-coloring for the triangular offset-hex lattice:
 *   color_2d = (c + (r & 1)) % 3
 * Proof: for even row r, the six hex neighbours have offsets
 *   (dr,dc) ∈ {(±1,0),(±1,+1),(0,±1)}.  Substituting into (c+dc+(r+dr)&1)%3
 *   gives all three colours ≠ color_2d(r,c).  Verified analogously for odd r.
 *
 * Combined 6-coloring (3 hex × 2 z):
 *   color = (c + (r & 1)) % 3  +  (z & 1) * 3
 * Each of the 8 neighbours has a strictly different color → in-place GS update
 * is race-free within a single kernel launch.
 *
 * SOR relaxation ω = 2/(1+sin(π/L)) drives convergence O(L²)→O(L):
 *   L=161: ~200 passes (was 25 000 Jacobi)
 *   L=481: ~350 passes (was 216 000 Jacobi) — 600× speedup
 *
 * Usage: 6 sequential kernel launches per iteration (phases 0..5).
 * Kernels in the same CUDA stream are ordered, so no explicit sync between
 * phases — the stream barrier ensures phase k writes are visible to phase k+1.
 */
__global__
void sor_obc_inplace_kernel(double* phi, const double* rho, int phase, double omega,
                              int n, int hex_rows, int hex_cols, int layers,
                              double cx, double cy, double cz, double alpha_em)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int lsz = hex_rows * hex_cols;
    int z   = i / lsz, rem = i % lsz;
    int r   = rem / hex_cols, c = rem % hex_cols;
    int color = (c + (r & 1)) % 3  +  (z & 1) * 3;

    if (color != phase) return;   /* this thread sits out this phase */

    int idx[8], nr[8], nc[8], nz[8];
    hex_nbrs_obc(i, hex_rows, hex_cols, layers, idx, nr, nc, nz);

    double s = 0.0;
    for (int j = 0; j < 8; j++) {
        s += (idx[j] >= 0) ? phi[idx[j]]
                            : obc_phi(nr[j], nc[j], nz[j],
                                      hex_rows, hex_cols, layers,
                                      cx, cy, cz, alpha_em);
    }
    double phi_gs = (s + 8.0*rho[i]) / 8.0;
    phi[i] = (1.0 - omega)*phi[i] + omega*phi_gs;   /* in-place SOR update */
}

/* ── Compact SOR: map thread→site for one color only ────────────────────────
 *
 * Each color contains exactly n_color ≈ N/6 sites.  We launch only n_color
 * threads and reconstruct (z,r,c) from the compact index.
 *
 * Color decomposition:  color = (c + (r&1)) % 3 + (z&1) * 3
 *
 * For a given phase ∈ {0..5}:
 *   z_parity  = phase / 3          (0 or 1)
 *   hex_color = phase % 3          (0, 1, or 2)
 *
 * Within layers of the correct z-parity, we enumerate sites whose
 * hex_color matches.  For each row r in the layer, the matching columns
 * satisfy (c + (r&1)) % 3 == hex_color, i.e. c ≡ (hex_color - (r&1)) mod 3.
 * There are cols_in_row = ceil((hex_cols - first_c) / 3) such columns.
 *
 * We precompute n_color on the host, launch ceil(n_color/BLOCK) blocks,
 * and each thread reconstructs its (z,r,c) from a flat compact index.
 */
__global__
void sor_obc_compact_kernel(double* phi, const double* rho,
                             int phase, double omega,
                             int hex_rows, int hex_cols, int layers,
                             int n_color, int cols_per_row,
                             double cx, double cy, double cz, double alpha_em)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= n_color) return;

    int z_par  = phase / 3;
    int hc     = phase % 3;

    /* Decompose compact index → (half_z, r, col_idx) */
    int sites_per_layer = hex_rows * cols_per_row;
    int half_z = tid / sites_per_layer;
    int rem    = tid % sites_per_layer;
    int r      = rem / cols_per_row;
    int ci     = rem % cols_per_row;

    /* Reconstruct actual z and c */
    int z = half_z * 2 + z_par;
    if (z >= layers || r >= hex_rows) return;
    int first_c = ((hc - (r & 1)) % 3 + 3) % 3;
    int c = first_c + ci * 3;
    if (c >= hex_cols) return;

    int lsz = hex_rows * hex_cols;
    int i   = (z * hex_rows + r) * hex_cols + c;

    int idx[8], nr_a[8], nc_a[8], nz_a[8];
    hex_nbrs_obc(i, hex_rows, hex_cols, layers, idx, nr_a, nc_a, nz_a);

    double s = 0.0;
    for (int j = 0; j < 8; j++) {
        s += (idx[j] >= 0) ? phi[idx[j]]
                            : obc_phi(nr_a[j], nc_a[j], nz_a[j],
                                      hex_rows, hex_cols, layers,
                                      cx, cy, cz, alpha_em);
    }
    double phi_gs = (s + 8.0*rho[i]) / 8.0;
    phi[i] = (1.0 - omega)*phi[i] + omega*phi_gs;
}

/* OBC Hamiltonian + imaginary-time step — boundary ψ = 0 */
__global__
void ham_obc_kernel(const C64* psi, C64* psi_new, const double* phi,
                     double alpha_em, double charge, double dtau,
                     int n, int hex_rows, int hex_cols, int layers)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= n) return;

    int idx[8], nr[8], nc[8], nz[8];
    hex_nbrs_obc(i, hex_rows, hex_cols, layers, idx, nr, nc, nz);

    C64 ns = {0.0, 0.0};
    for (int j = 0; j < 8; j++)
        if (idx[j] >= 0) ns = cadd(ns, psi[idx[j]]);  /* else ψ=0 */

    C64 kinetic   = csub(psi[i], cmulr(ns, 0.125));
    C64 potential = cmulr(psi[i], alpha_em*charge*phi[i]);
    C64 h         = cadd(kinetic, potential);
    psi_new[i]    = csub(psi[i], cmulr(h, dtau));
}

/* OBC energy kernel */
__global__
void energy_obc_kernel(const C64* psi, const double* phi,
                        double alpha_em, double charge,
                        double* out_total, double* out_kin, double* out_pot,
                        int n, int hex_rows, int hex_cols, int layers)
{
    extern __shared__ double sdata[];
    int tid = threadIdx.x, i = blockIdx.x*blockDim.x + tid;

    double et=0.0, ek=0.0, ep=0.0;
    if (i < n) {
        int idx[8], nr[8], nc[8], nz[8];
        hex_nbrs_obc(i, hex_rows, hex_cols, layers, idx, nr, nc, nz);
        C64 ns = {0.0,0.0};
        for (int j = 0; j < 8; j++)
            if (idx[j] >= 0) ns = cadd(ns, psi[idx[j]]);
        C64 kinetic   = csub(psi[i], cmulr(ns, 0.125));
        C64 potential = cmulr(psi[i], alpha_em*charge*phi[i]);
        ek = re_conj_dot(psi[i], kinetic);
        ep = re_conj_dot(psi[i], potential);
        et = ek + ep;
    }
    double* sk = sdata; double* se = sdata+blockDim.x; double* sp = sdata+2*blockDim.x;
    sk[tid]=ek; se[tid]=et; sp[tid]=ep;
    __syncthreads();
    for (int s=blockDim.x/2; s>0; s>>=1) {
        if (tid<s) { sk[tid]+=sk[tid+s]; se[tid]+=se[tid+s]; sp[tid]+=sp[tid+s]; }
        __syncthreads();
    }
    if (tid==0) { out_kin[blockIdx.x]=sk[0]; out_total[blockIdx.x]=se[0]; out_pot[blockIdx.x]=sp[0]; }
}

/* ── OBC host driver ─────────────────────────────────────────────────────────
 *
 * Uses a small box of size hex_rows×hex_cols×layers with OPEN boundary
 * conditions: ψ=0 and φ=α/r at walls.  Equivalent to L=∞ periodic for a
 * localised state, but uses only n = hex_rows*hex_cols*layers allocations
 * regardless of the physical system size.
 */
extern "C"
void run_hydrogen_obc(
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

    double cx, cy, cz;
    { int r=hex_rows/2, c=hex_cols/2;
      hex_cart(r, c, &cx, &cy); cz=(double)(layers/2); }

    setbuf(stdout, NULL);  /* unbuffered — survive SIGKILL */

    /* Diagnostic: what does the runtime think is available? */
    {
        size_t mem_free = 0, mem_total = 0;
        cudaMemGetInfo(&mem_free, &mem_total);
        printf("  GPU memory: %.2f GiB free / %.2f GiB total\n",
               (double)mem_free / (1024.0*1024.0*1024.0),
               (double)mem_total / (1024.0*1024.0*1024.0));
        size_t need = (size_t)n * (2*sizeof(double) + 2*sizeof(C64));
        printf("  Need: %.2f GiB (rho+phi+psi+psi2)\n",
               (double)need / (1024.0*1024.0*1024.0));
    }

    double *d_rho, *d_phi, *d_pe, *d_pk, *d_pp;
    C64    *d_psi, *d_psi2;
    printf("  Allocating d_rho  (%.2f GiB)...\n", (double)((size_t)n*sizeof(double))/(1024.0*1024.0*1024.0));
    CUDA_CHECK(cudaMalloc(&d_rho,     (size_t)n * sizeof(double)));
    printf("  Allocating d_phi  (%.2f GiB)...\n", (double)((size_t)n*sizeof(double))/(1024.0*1024.0*1024.0));
    CUDA_CHECK(cudaMalloc(&d_phi,     (size_t)n * sizeof(double)));
    printf("  Allocating d_psi  (%.2f GiB)...\n", (double)((size_t)n*sizeof(C64))/(1024.0*1024.0*1024.0));
    CUDA_CHECK(cudaMalloc(&d_psi,     (size_t)n * sizeof(C64)));
    printf("  Allocating d_psi2 (%.2f GiB)...\n", (double)((size_t)n*sizeof(C64))/(1024.0*1024.0*1024.0));
    CUDA_CHECK(cudaMalloc(&d_psi2,    (size_t)n * sizeof(C64)));
    printf("  Allocating reduction buffers...\n");
    CUDA_CHECK(cudaMalloc(&d_pe,      n_blk *sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_pk,      n_blk *sizeof(double)));
    CUDA_CHECK(cudaMalloc(&d_pp,      n_blk *sizeof(double)));
    double* h_partial = (double*)malloc(n_blk * sizeof(double));
    printf("  All allocations succeeded!\n"); fflush(stdout);
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Phase 1: Coulomb warm-start + SOR correction.
     * φ = α/r is the exact continuum solution. SOR then corrects only the
     * lattice discretization error — residual ~ a²/r² — converging in ~50-100
     * sweeps instead of ~L sweeps from cold start (φ=0). */
    init_coulomb_kernel<<<n_blk,BLOCK>>>(d_phi, n, hex_rows, hex_cols, layers,
                                          cx, cy, cz, alpha_em);
    init_rho_kernel<<<n_blk,BLOCK>>>(d_rho, n, center);
    CUDA_CHECK(cudaDeviceSynchronize());
    {
        double omega = 2.0 / (1.0 + sin(M_PI / (double)hex_rows));
        /* Compact SOR: launch only N/6 threads per phase instead of N.
         * cols_per_row = max columns matching one hex_color in a row. */
        int cols_per_row = (hex_cols + 2) / 3;
        int half_layers  = (layers + 1) / 2;
        int n_color      = half_layers * hex_rows * cols_per_row;
        int n_color_blk  = (n_color + BLOCK - 1) / BLOCK;
        printf("  Phase 1: SOR Poisson (%d sweeps, ω=%.4f, %d threads/phase vs %d total)...\n",
               n_jacobi, omega, n_color, n); fflush(stdout);
        for (int iter = 0; iter < n_jacobi; iter++) {
            for (int phase = 0; phase < 6; phase++) {
                sor_obc_compact_kernel<<<n_color_blk,BLOCK>>>(
                    d_phi, d_rho, phase, omega,
                    hex_rows, hex_cols, layers,
                    n_color, cols_per_row,
                    cx, cy, cz, alpha_em);
            }
            if ((iter+1) % 100 == 0 || iter == 0) {
                CUDA_CHECK(cudaDeviceSynchronize());
                printf("    SOR sweep %d/%d\n", iter+1, n_jacobi); fflush(stdout);
            }
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    printf("  Phase 1 done.\n"); fflush(stdout);

    /* d_rho is dead after Poisson — free it to save N×8B for large L.
     * L=961: saves 7.1 GB, bringing peak from 43 GB to 36 GB. */
    cudaFree(d_rho); d_rho = NULL;
    printf("  Freed d_rho (%.2f GiB reclaimed).\n", (double)((size_t)n*sizeof(double))/(1024.0*1024.0*1024.0)); fflush(stdout);

    /* Phase 2: Gaussian init centred on proton */
    printf("  Phase 2: Gaussian init (σ=%.1f)...\n", 1.0/alpha_em > 1.0 ? 1.0/alpha_em : 1.0);
    double sigma = 1.0/alpha_em; if (sigma<1.0) sigma=1.0;
    init_gaussian_kernel<<<n_blk,BLOCK>>>(d_psi, n, hex_rows, hex_cols, layers,
                                           cx, cy, cz, sigma);
    CUDA_CHECK(cudaDeviceSynchronize());

    norm_sq_kernel<<<n_blk,BLOCK,BLOCK*sizeof(double)>>>(d_psi, d_pe, n);
    CUDA_CHECK(cudaMemcpy(h_partial, d_pe, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
    { double s=0.0; for(int i=0;i<n_blk;i++) s+=h_partial[i];
      scale_kernel<<<n_blk,BLOCK>>>(d_psi, 1.0/sqrt(s), n); }
    CUDA_CHECK(cudaDeviceSynchronize());
    printf("  Phase 2 done.\n"); fflush(stdout);

    /* Phase 3: OBC imaginary-time evolution */
    printf("  Phase 3: Imaginary time (%d steps, dtau=%.3f, τ=%.1f)...\n",
           n_iter, dtau, n_iter*dtau); fflush(stdout);
    for (int iter=0; iter<n_iter; iter++) {
        ham_obc_kernel<<<n_blk,BLOCK>>>(d_psi, d_psi2, d_phi,
            alpha_em, -1.0, dtau, n, hex_rows, hex_cols, layers);
        C64* tmp=d_psi; d_psi=d_psi2; d_psi2=tmp;
        if ((iter+1)%renorm_every==0) {
            norm_sq_kernel<<<n_blk,BLOCK,BLOCK*sizeof(double)>>>(d_psi, d_pe, n);
            CUDA_CHECK(cudaMemcpy(h_partial, d_pe, n_blk*sizeof(double), cudaMemcpyDeviceToHost));
            double s=0.0; for(int i=0;i<n_blk;i++) s+=h_partial[i];
            scale_kernel<<<n_blk,BLOCK>>>(d_psi, 1.0/sqrt(s), n);
            if ((iter+1) % 1000 == 0) {
                printf("    Step %d/%d\n", iter+1, n_iter); fflush(stdout);
            }
        }
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    printf("  Phase 3 done.\n"); fflush(stdout);

    /* Phase 4: energy */
    energy_obc_kernel<<<n_blk,BLOCK,3*BLOCK*sizeof(double)>>>(
        d_psi, d_phi, alpha_em, -1.0,
        d_pe, d_pk, d_pp, n, hex_rows, hex_cols, layers);
    CUDA_CHECK(cudaDeviceSynchronize());

    double et=0,ek=0,ep=0;
    CUDA_CHECK(cudaMemcpy(h_partial,d_pe,n_blk*sizeof(double),cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) et+=h_partial[i];
    CUDA_CHECK(cudaMemcpy(h_partial,d_pk,n_blk*sizeof(double),cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) ek+=h_partial[i];
    CUDA_CHECK(cudaMemcpy(h_partial,d_pp,n_blk*sizeof(double),cudaMemcpyDeviceToHost));
    for(int i=0;i<n_blk;i++) ep+=h_partial[i];

    *out_e_total=et; *out_e_kin=ek; *out_e_pot=ep;

    /* d_rho already freed after Poisson */
    cudaFree(d_phi);
    cudaFree(d_psi); cudaFree(d_psi2);
    cudaFree(d_pe); cudaFree(d_pk); cudaFree(d_pp);
    free(h_partial);
}

/* ── Host driver ─────────────────────────────────────────────────────────────*/

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
    CUDA_CHECK(cudaMemset(d_phi, 0, (size_t)n * sizeof(double)));
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
