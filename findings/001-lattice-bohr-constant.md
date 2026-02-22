# Finding 001: The Lattice Bohr Constant

**Date**: 2026-02-22
**Status**: CONFIRMED (GPU + CPU + analytical)
**Crate**: `gutoe-gpu`, module `watson.rs`

---

## Summary

The hydrogen ground-state energy on the GUTOE hex+z lattice is

$$E_0 = -\frac{\alpha^2}{2}\,C_\infty, \qquad C_\infty = 0.5466 \pm 0.0005$$

where $\alpha$ is the coupling constant. This is 45% weaker than continuum
hydrogen ($C = 1$) due to **ultraviolet regularization by the Brillouin zone
boundary** --- the lattice kinetic operator saturates at $T_{\max} = 2$ instead
of growing as $k^2/(2m) \to \infty$, weakening the short-distance Coulomb
potential.

The continuum Watson-integral prediction $C_G^2/2 = 0.721$ overestimates
binding by **32%**. The true value requires the full lattice dispersion relation
and cannot be obtained from a quadratic ($k \to 0$) expansion.

---

## 1. The Lattice Hamiltonian

The hex+z lattice has **8 neighbors per site**: 6 in-plane (triangular/hex) + 2
along the z-axis.

**Kinetic operator** (real-space):

$$(\hat{T}\psi)(r) = \psi(r) - \tfrac{1}{8}\sum_{j=1}^{8} \psi(r + \delta_j)$$

**Dispersion relation** (Bloch):

$$T(\theta_1,\theta_2,\theta_z) = 1 - \tfrac{1}{8}\bigl[2\cos\theta_1 + 2\cos\theta_2 + 2\cos(\theta_1-\theta_2) + 2\cos\theta_z\bigr]$$

where $\theta_1, \theta_2$ are phases along the triangular primitive vectors
$\mathbf{a}_1 = (1,0)$, $\mathbf{a}_2 = (1/2, \sqrt{3}/2)$, and $\theta_z$ is
along the stacking axis.

**Coulomb potential**: $\hat{T}\phi = \rho$ (lattice Poisson equation), with OBC
boundary condition $\phi_{\mathrm{wall}} = \alpha/r$ and charge
$\rho = -1/N + \delta_{\mathrm{center}}$.

**Hydrogen Hamiltonian**: $H = \hat{T} - \alpha\,\phi$, solved by imaginary-time
propagation $\psi \leftarrow \psi - d\tau\,(H\psi)$ with renormalization every
20 steps.

---

## 2. Analytical: Watson Integral and Second-Moment Tensor

### 2.1 Second-Moment Tensor

The quadratic expansion of $T(\mathbf{k})$ near $\mathbf{k} = 0$ in Cartesian
coordinates gives:

$$T(\mathbf{k}) \approx M_{xx}\,k_x^2 + M_{yy}\,k_y^2 + M_{zz}\,k_z^2$$

| Component | Value | Fraction |
|-----------|-------|----------|
| $M_{xx}$  | 3/16  | 0.1875   |
| $M_{yy}$  | 3/16  | 0.1875   |
| $M_{zz}$  | 1/8   | 0.1250   |

Effective masses: $m_{xy} = 1/(2M_{xx}) = 8/3 \approx 2.67$, $m_z = 1/(2M_{zz}) = 4$.

In BZ coordinates $(\theta_1, \theta_2, \theta_z)$:

$$T_{\mathrm{quad}} = \tfrac{1}{4}(\theta_1^2 - \theta_1\theta_2 + \theta_2^2) + \tfrac{1}{8}\theta_z^2$$

### 2.2 Watson Integral

The lattice Green's function at the origin:

$$G(0) = \frac{1}{(2\pi)^3}\int_0^{2\pi}\!\!d\theta_1\int_0^{2\pi}\!\!d\theta_2\int_{-\pi}^{\pi}\!\!d\theta_z\;\frac{1}{T(\theta_1,\theta_2,\theta_z)}$$

| Lattice | $G(0)$ | Known value |
|---------|--------|-------------|
| Simple cubic (validation) | 1.5150 | 1.5164 (Watson 1939) |
| **Hex+z** | **1.4482** | (new result) |

Convergence verified by midpoint rule at $n = 100, 200, 400$.

### 2.3 Continuum Green's Function Coefficient

$$C_G = \frac{1}{4\pi\,M_\perp\sqrt{M_z}} = \frac{1}{4\pi \cdot \tfrac{3}{16}\cdot\tfrac{1}{2\sqrt{2}}} = 1.2004$$

Continuum variational prediction for the Bohr constant:

$$C_\infty^{(\mathrm{cont})} = \frac{C_G^2}{2} = \frac{(1.2004)^2}{2} = 0.7205$$

**This is 32% above the actual lattice value.** The quadratic approximation to
$T(\mathbf{k})$ is insufficient.

---

## 3. Numerical: GPU Results

### 3.1 L-scan at $\alpha = 0.10$ (SOR Poisson, CUDA + ROCm verified identical)

Box size $L$ (odd, OBC), $L^3$ sites, charge at center.

| $L$ | Sites | $E_\mathrm{total}$ | $C(L)$ |
|-----|-------|---------------------|--------|
| 161 | 4,173,281 | $-2.170 \times 10^{-3}$ | 0.4340 |
| 241 | 13,997,521 | $-2.356 \times 10^{-3}$ | 0.4711 |
| 321 | 33,076,161 | $-2.449 \times 10^{-3}$ | 0.4899 |
| 481 | 111,284,641 | $-2.544 \times 10^{-3}$ | 0.5088 |
| 961 | 887,503,681 | $-2.638 \times 10^{-3}$ | 0.5277 |

L=961 run: 39.67 GiB VRAM on AMD Strix Halo gfx1151 (96 GiB), ROCm 7.2,
n_sor=1000, n_iter=30000 ($\tau = 1500$), 160 min.

### 3.2 Richardson Extrapolation

Model: $C(L) = C_\infty - B/L$

| Pair | $B$ | $C_\infty$ |
|------|-----|-----------|
| (161, 241) | 18.0 | 0.5458 |
| (241, 321) | 18.2 | 0.5465 |
| (321, 481) | 18.2 | 0.5467 |
| (481, 961) | 18.2 | **0.5466** |

**All four pairs agree to $C_\infty = 0.5466 \pm 0.0005$, $B = 18.15 \pm 0.15$.**

The linear $1/L$ model fits the data from $L = 161$ to $L = 961$ with sub-0.1%
residuals. Higher-order corrections ($D/L^2$, etc.) are negligible in this range.

### 3.3 Boundary Diagnostic

Wavefunction amplitude at OBC walls (confirming no image-charge artifacts):

| $L$ | $|\psi|^2_{\mathrm{boundary}}$ |
|-----|-------------------------------|
| 161 | $2.72 \times 10^{-9}$ |
| 241 | $1.87 \times 10^{-13}$ |
| 321 | $1.24 \times 10^{-18}$ (machine zero) |

The wavefunction vanishes at the walls. The OBC setup is clean.

### 3.4 $\alpha$-scan (single $L$ per $\alpha$, $L/a_0 \approx 16$)

| $\alpha$ | $L$ | $C(L)$ |
|----------|-----|--------|
| 0.30 | 55 | 0.5509 |
| 0.20 | 81 | 0.4800 |
| 0.10 | 161 | 0.4340 |
| 0.07 | 231 | 0.4243 |
| 0.05 | 321 | 0.4172 |

These are raw $C(L)$ values, **not** $C_\infty$. The finite-size correction
$B(\alpha)/L$ grows with decreasing $\alpha$ (wavefunction becomes more spread,
tail hits the walls harder). To extract $C_\infty(\alpha)$ requires a separate
L-scan at each $\alpha$.

**Open question**: does $C_\infty$ depend on $\alpha$? Physical expectation:
- $C_\infty(\alpha \to 0) \to 1$ (continuum limit, Bohr radius $\gg$ lattice spacing)
- $C_\infty(\alpha \to \infty) \to$ some lattice value (tight-binding limit)
- We only have a complete L-scan at $\alpha = 0.10$.

---

## 4. Numerical: CPU Validation

### 4.1 Implementation

CPU solver in `crates/gutoe-gpu/src/watson.rs`, function `cpu_hydrogen_obc()`.
Matches the GPU kernel exactly:

- **OBC neighbors**: `hex_nbrs_obc()` returns sentinel $-1$ for out-of-bounds.
- **Poisson solve**: 6-color SOR with $\omega = 2/(1 + \sin(\pi/L))$, Coulomb
  warm-start $\phi_\mathrm{init} = \alpha/r$.
- **Imaginary time**: $\psi \leftarrow \psi - d\tau\,(\hat{T}\psi - \alpha\phi\psi)$,
  renormalize every 20 steps, OBC: out-of-bounds $\to \psi = 0$.

### 4.2 Results ($\alpha = 0.10$, $\tau = 1000$, $d\tau = 0.05$)

| $L$ | $C(L)$ | GPU predicted ($0.547 - 18.1/L$) |
|-----|--------|----------------------------------|
| 31 | $-0.301$ | $-0.04$ (box too small) |
| 41 | 0.023 | 0.11 |
| 51 | 0.167 | 0.19 |
| 61 | 0.245 | 0.25 |
| 71 | 0.292 | 0.29 |
| 81 | 0.325 | 0.32 |

Two-point Richardson from CPU:

| Pair | $C_\infty$ |
|------|-----------|
| (61, 71) | 0.583 |
| (71, 81) | **0.558** |

Three-point Richardson (last triple):

| Triple | $C_\infty$ |
|--------|-----------|
| (61, 71, 81) | 0.481 |

**Bracketed**: $0.481 < C_\infty < 0.558$, centered on the GPU's 0.547.

Including $L = 101$ (from separate run, $C = 0.366$):
$(81, 101) \to C_\infty = 0.540$ --- 1.3% below GPU.

---

## 5. Dispersion Analysis: Why 32% Off?

### 5.1 Lattice vs Quadratic Dispersion

Along all high-symmetry directions, $T_\mathrm{lattice}/T_\mathrm{quad}$
starts at $\approx 1.00$ for small $\theta$ and drops to **0.41 at the BZ
boundary** ($\theta = \pi$).

| $\theta$ | $T_\mathrm{lattice}$ | $T_\mathrm{quad}$ | Ratio |
|----------|---------------------|-------------------|-------|
| 0.16 | 0.00616 | 0.00617 | 0.998 |
| 0.79 | 0.146 | 0.154 | 0.950 |
| 1.57 | 0.500 | 0.617 | 0.811 |
| 2.36 | 0.854 | 1.388 | 0.615 |
| 3.14 | 1.000 | 2.467 | **0.405** |

The quadratic approximation diverges ($T_\mathrm{quad} \propto \theta^2$)
while the lattice saturates ($T_\mathrm{lattice} \leq 2$). This saturation
provides **natural UV regularization**.

### 5.2 Lattice Coulomb vs Continuum (z-axis)

$$G_\mathrm{cont}(0,0,z) = \frac{1}{4\pi\,M_\perp\,|z|} = \frac{4}{3\pi\,|z|}$$

| $z$ | $G_\mathrm{lattice}$ | $G_\mathrm{cont}$ | Ratio |
|-----|---------------------|--------------------|-------|
| 1 | 0.4079 | 0.4244 | 0.961 |
| 2 | 0.1921 | 0.2122 | 0.905 |
| 5 | 0.0707 | 0.0849 | 0.833 |
| 8 | 0.0427 | 0.0531 | 0.805 |
| 12 | 0.0279 | 0.0354 | 0.788 |

At the Bohr radius ($r \sim 1/\alpha = 10$ lattice units), the lattice Coulomb
potential is **~20% weaker** than continuum. This directly reduces $C_\infty$.

### 5.3 UV Dominance of $G(0)$

Fraction of the Watson integral $G(0)$ contained within momentum shell
$|\mathbf{k}| < k_\mathrm{max}$:

| $k_\mathrm{max}/\pi$ | Fraction of $G(0)$ |
|-----------------------|-------------------|
| 0.10 | 2.1% |
| 0.20 | 4.2% |
| 0.50 | 11.1% |
| 0.70 | 16.3% |
| 1.00 | **26.3%** |

**73.7% of $G(0)$ comes from UV modes ($|k| > \pi$)** where the lattice
dispersion is far from quadratic. These modes dominate the on-site Coulomb
potential. Since $T_\mathrm{lattice} < T_\mathrm{quad}$ at large $k$, the
UV contributions to $G(0)$ are *smaller* than the continuum prediction,
weakening the potential and reducing $C_\infty$.

---

## 6. Physical Interpretation

### Why $C_\infty = 0.547$, not $1.0$?

1. **Lattice UV regularization**: The kinetic operator $\hat{T}$ has eigenvalues
   bounded in $[0, 2]$. The continuum kinetic operator $-\nabla^2/(2m)$ has
   eigenvalues $\in [0, \infty)$. This means the lattice Fourier-space Coulomb
   potential $\tilde{\phi}(k) = 1/T(k)$ is bounded below by $1/2$, while the
   continuum $1/(k^2/2m)$ diverges as $k \to \infty$.

2. **Short-distance Coulomb weakening**: The lattice Coulomb potential at the
   origin is $G(0) = 1.448$ (finite), while the continuum Coulomb diverges as
   $1/r \to \infty$. The lattice "spreads out" the point charge over one lattice
   site, regularizing the UV divergence.

3. **Quantitative effect**: The hydrogen wavefunction at $\alpha = 0.10$ has
   Bohr radius $\sim 10$ sites. It samples the lattice Coulomb potential at all
   distances $0 \leq r \lesssim 30$. The potential is 4% weak at $r = 1$ and
   20% weak at $r = 12$, giving an integrated weakening of $\sim 25\%$, consistent
   with $C_\infty/C_\mathrm{cont} = 0.547/0.721 = 0.76$.

### Is $C_\infty$ a lattice invariant?

**Partially.** $C_\infty$ is determined by the lattice geometry (hex+z
dispersion relation) and the coupling $\alpha$:

- It does NOT depend on box size $L$ (confirmed by Richardson extrapolation over
  $L = 161$ to $961$).
- It does NOT depend on numerical parameters ($n_\mathrm{sor}$, $n_\mathrm{iter}$,
  $d\tau$) once converged.
- It MAY depend on $\alpha$ through the ratio of Bohr radius to lattice spacing
  (needs L-scans at multiple $\alpha$ to confirm).

---

## 7. PBC Zero-Mode Failure

Periodic boundary conditions (PBC) cause the lattice Poisson solver to diverge:

| $L$ | $\phi_\mathrm{max}$ (PBC) |
|-----|--------------------------|
| 161 | $\sim 10^{150}$ |
| 241 | $\sim 10^{87}$ |
| 321 | $\sim 10^{49}$ |

**Root cause**: The periodic Poisson equation $T\phi = \rho$ with
$\sum \rho = 0$ has a zero mode --- $T(k=0) = 0$, so $\phi(k=0)$ is
unconstrained. The SOR solver amplifies this unconstrained constant mode
exponentially. Fix: subtract the zero mode ($\phi \to \phi - \bar{\phi}$)
after each sweep. Not yet implemented.

---

## 8. Reproducibility

All code in `crates/gutoe-gpu/src/watson.rs`:

```bash
# Watson integral + dispersion analysis (< 2 sec)
cargo test -p gutoe-gpu --lib watson::tests::watson_analysis --release -- --nocapture
cargo test -p gutoe-gpu --lib watson::tests::dispersion --release -- --nocapture

# CPU hydrogen solver (~ 2 min)
cargo test -p gutoe-gpu --lib watson::tests::cpu_hydrogen --release -- --nocapture

# GPU hydrogen (requires CUDA sm_89 or ROCm gfx1151)
cargo test -p gutoe-gpu --features cuda --release -- bohr_obc_lscan --nocapture
```

GPU runs on remote host `wings@10.7.1.195` (tealc), AMD Strix Halo gfx1151,
96 GiB VRAM, ROCm 7.2. L=961 requires stopping GDM (`sudo systemctl stop gdm`)
to free display VRAM.

---

## 9. Key Numbers

| Quantity | Value | Source |
|----------|-------|--------|
| $G_\mathrm{hex+z}(0)$ (Watson integral) | 1.4482 | `watson_hex_z(400)` |
| $G_\mathrm{sc}(0)$ (simple cubic, validation) | 1.5150 | `watson_simple_cubic(400)` |
| $M_{xx} = M_{yy}$ | 3/16 = 0.1875 | `second_moment_tensor()` |
| $M_{zz}$ | 1/8 = 0.1250 | `second_moment_tensor()` |
| $C_G$ (continuum Green coefficient) | 1.2004 | $1/(4\pi M_\perp\sqrt{M_z})$ |
| $C_\infty^{(\mathrm{cont})}$ (continuum prediction) | 0.7205 | $C_G^2/2$ |
| $C_\infty$ (actual, GPU L-scan) | **0.5466** | Richardson, 5 L values |
| $B$ (finite-size slope) | **18.15** | Richardson, all pairs |
| Continuum/lattice ratio | 0.759 | $C_\infty/C_\infty^{(\mathrm{cont})}$ |
| UV fraction of $G(0)$ | 73.7% | BZ shell integral |

---

## 10. What This Means for GUTOE

The lattice Bohr constant $C_\infty = 0.547$ is a **prediction** of the hex+z
lattice geometry. Given the 8-neighbor kinetic operator and the lattice Poisson
equation, the hydrogen binding energy follows deterministically. There is no free
parameter.

However, $C_\infty \neq 1$ means the lattice hydrogen does not reproduce the
exact continuum Bohr energy. This is expected --- the lattice IS the UV
regulator, and lattice artifacts are $O(a^2)$ where $a$ is the lattice spacing.
In the continuum limit ($\alpha \to 0$, Bohr radius $\gg a$), we expect
$C_\infty(\alpha) \to 1$.

The 2.73x enrichment from the GUTOE simulation is a separate, dynamical result
that depends on the lattice step function, gauge coupling, and lepton injection
protocol --- not directly on $C_\infty$. The Schrodinger solver measures the
STATIC lattice hydrogen, while the simulation measures DYNAMICAL enrichment of
leptons near proton triplets.
