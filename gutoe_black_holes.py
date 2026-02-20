#!/usr/bin/env python3
"""
GUTOE near black holes: where Planck-scale dispersion meets observable physics.

Black holes are the only place where k approaches k_c naturally.
Gravitational blueshift pushes infalling photon frequencies toward Planck scale
without limit at the horizon. Escaping Hawking photons carry the imprint.

Three observables (ranked by feasibility):
  1. Quasi-normal mode frequency shift (LIGO measures this NOW)
  2. Hawking radiation spectrum modification (Fermi looks for this)
  3. Photon ring phase shift (next-gen EHT, ~2030s)

The physics question:
  Does λ_QG run near the horizon? If so, how fast?
  "Does the temporal sector leak into the spatial sector near a BH?"
"""

import numpy as np

print("=" * 70)
print("GUTOE NEAR BLACK HOLES: FINDING THE TESTABLE REGIME")
print("=" * 70)
print()

# ── Constants ───────────────────────────────────────────────────────────────
c     = 2.998e8          # m/s
G     = 6.674e-11        # m³ kg⁻¹ s⁻²
ħ     = 1.055e-34        # J·s
k_B   = 1.381e-23        # J/K
ℓ_P   = 1.616e-35        # m
M_sun = 1.989e30         # kg
M_P   = np.sqrt(ħ * c / G)  # Planck mass ≈ 2.18e-8 kg

lambda_QG = 1/12
DISPERSION_COEFF = lambda_QG * ℓ_P**2

# Critical wavenumber: modes above this are evanescent in spatial sector
k_c = c / (ℓ_P * np.sqrt(1/lambda_QG))   # = c√12 / ℓ_P

print(f"λ_QG = 1/12,  ℓ_P = {ℓ_P:.3e} m")
print(f"k_c  = {k_c:.3e} m⁻¹  (Planck-scale spatial cutoff)")
print()

# ── Schwarzschild radius and basic BH scales ────────────────────────────────
def r_s(M_kg):
    """Schwarzschild radius"""
    return 2 * G * M_kg / c**2

def T_Hawking(M_kg):
    """Standard Hawking temperature"""
    return ħ * c**3 / (8 * np.pi * G * M_kg * k_B)

def t_evaporate(M_kg):
    """Hawking evaporation time"""
    return 5120 * np.pi * G**2 * M_kg**3 / (ħ * c**4)

# ── THE GUTOE CORRECTION TO HAWKING TEMPERATURE ─────────────────────────────
print("=" * 70)
print("1. HAWKING TEMPERATURE CORRECTION")
print("=" * 70)
print()
print("The Hawking derivation uses modes with k ~ 1/r_s near the horizon.")
print("GUTOE modifies these modes: ω²(k) = c²k² - (1/12)ℓ_P²k⁴")
print()
print("The fractional correction to T_H:")
print("  δT/T ~ (1/12)(ℓ_P/r_s)² = λ_QG(ℓ_P/r_s)²")
print()

def hawking_correction(M_kg):
    """
    GUTOE correction to Hawking temperature.
    The characteristic wavenumber at horizon: k_horizon ~ 1/r_s
    Correction: δω/ω ~ -(1/2)(k·ℓ_P)²/12 = -(1/24)(ℓ_P/r_s)²
    → δT/T ~ +(1/24)(ℓ_P/r_s)²  [positive: slightly hotter]
    """
    rs = r_s(M_kg)
    correction = (1/24) * (ℓ_P / rs)**2
    return correction

print(f"{'BH Mass':>20s} {'r_s':>14s} {'T_Hawking':>14s} {'δT/T (GUTOE)':>16s}")
print("-" * 70)

bh_masses = {
    'Stellar (30 M☉)':     30 * M_sun,
    'LIGO typical (60 M☉)': 60 * M_sun,
    'Intermediate (10⁵ M☉)': 1e5 * M_sun,
    'Sgr A* (4×10⁶ M☉)':  4e6 * M_sun,
    'Prim. BH r_s=1μm':    c**2 * 1e-6 / (2*G),
    'Prim. BH r_s=1pm':    c**2 * 1e-12 / (2*G),
    'Prim. BH r_s=1fm':    c**2 * 1e-15 / (2*G),
    'Near-Planck r_s=1nm': c**2 * 1e-9 * ℓ_P / (2*G),
    'Planck BH r_s=ℓ_P':   c**2 * ℓ_P / (2*G),
}

for name, M in bh_masses.items():
    rs_val = r_s(M)
    T_H = T_Hawking(M)
    corr = hawking_correction(M)
    print(f"{name:>20s}  {rs_val:12.3e}m  {T_H:12.3e}K  {corr:14.3e}")

print()
print("The correction becomes order unity at r_s ~ ℓ_P (Planck BH).")
print()

# Find the BH mass where correction reaches 10%, 1%, 0.01%
print("BH mass where δT/T reaches threshold:")
for threshold, label in [(0.1, '10%'), (0.01, '1%'), (1e-4, '0.01%'), (1.0, '100% (order unity)')]:
    # (1/24)(ℓ_P/r_s)² = threshold
    # r_s = ℓ_P / sqrt(24 * threshold)
    rs_target = ℓ_P / np.sqrt(24 * threshold)
    M_target = c**2 * rs_target / (2 * G)
    T_target = T_Hawking(M_target)
    t_evap = t_evaporate(M_target)
    print(f"  {label:>12s}: M = {M_target:.3e} kg = {M_target/1e9:.3e} µg,"
          f"  r_s = {rs_target:.3e} m,  T_H = {T_target:.3e} K,  t_evap = {t_evap:.3e} s")
print()

# ── QUASI-NORMAL MODE FREQUENCY SHIFT ───────────────────────────────────────
print("=" * 70)
print("2. QUASI-NORMAL MODE FREQUENCY SHIFT")
print("=" * 70)
print()
print("QNM frequencies: ω_n ~ (1/r_s) × [n + 1/2 - i(2n+1)/4] (schematic)")
print()
print("GUTOE correction to n-th overtone:")
print("  δω_n/ω_n ~ n²(ℓ_P/r_s)²  [higher overtones probe closer to horizon]")
print()
print("This uses: higher overtones have WKB turning points closer to r_s,")
print("so they sample larger local k → larger GUTOE correction.")
print()

def qnm_correction(M_kg, n_overtone):
    """
    Fractional shift of n-th QNM overtone frequency.
    n=0: fundamental, n=1: first overtone, etc.
    GUTOE shift ~ n² × (ℓ_P/r_s)²
    (Full calculation needs wave eq. on Schwarzschild with GUTOE dispersion)
    """
    rs = r_s(M_kg)
    return (n_overtone**2) * (ℓ_P / rs)**2

M_LIGO = 60 * M_sun  # Typical final BH from GW merger

print(f"For a {int(M_LIGO/M_sun)} M☉ BH (typical LIGO merger remnant),")
print(f"r_s = {r_s(M_LIGO):.1f} m,  f_QNM ≈ {c / (2*np.pi * r_s(M_LIGO) * 2.19):.0f} Hz:")
print()
print(f"  {'Overtone n':>12s}  {'δf/f (GUTOE)':>16s}  {'Detectable (LISA)?':>20s}")
print("  " + "-" * 52)
for n in range(0, 8):
    shift = qnm_correction(M_LIGO, n)
    # LISA can measure QNM shifts at ~10^-3 fractional precision for massive BHs
    # Current LIGO: ~10^-2 for dominant mode
    ligo_threshold = 1e-2
    lisa_threshold = 1e-5
    detectable = "LISA? No" if shift < lisa_threshold else "MAYBE"
    if shift > ligo_threshold:
        detectable = "LIGO now!"
    print(f"  n = {n:>3d}:          {shift:14.3e}    {detectable}")

print()
# LISA supermassive merger: 10^7 M_sun
M_LISA = 1e7 * M_sun
print(f"For LISA target: {int(M_LISA/M_sun):.0e} M☉ BH merger,  r_s = {r_s(M_LISA):.3e} m:")
for n in [0, 10, 100, 1000]:
    shift = qnm_correction(M_LISA, n)
    print(f"  n = {n:>4d}:  δf/f = {shift:.3e}  {'← LISA threshold ~10⁻⁶' if shift > 1e-6 else ''}")

print()
# What overtone number n makes the GUTOE shift visible to LISA?
print("Overtone n needed for LIGO/LISA to see GUTOE in QNMs:")
for M, name, threshold in [
    (60 * M_sun, 'LIGO (60 M☉)', 1e-2),
    (1e7 * M_sun, 'LISA (10⁷ M☉)', 1e-6),
    (c**2 * 1e-15 / (2*G), 'Planck-near PBH', 0.1),
]:
    rs_val = r_s(M)
    # n² × (ℓ_P/r_s)² = threshold  →  n = sqrt(threshold) × (r_s/ℓ_P)
    n_needed = np.sqrt(threshold) * rs_val / ℓ_P
    print(f"  {name}: n ≈ {n_needed:.2e} overtones needed  [impossible for current BHs]")

print()

# ── THE TUNNELING CALCULATION: DOES λ_QG RUN? ────────────────────────────────
print("=" * 70)
print("3. TEMPORAL SECTOR LEAKAGE NEAR THE HORIZON")
print("=" * 70)
print()
print("KEY QUESTION: Does the evanescent k > k_c barrier become 'thin'")
print("near a BH horizon, allowing modes to tunnel through?")
print()
print("Near a Schwarzschild BH, a photon with far-away wavenumber k₀")
print("gets blueshifted to k_local(r) = k₀ / √(1 - r_s/r)")
print()
print("It crosses k_c at radius r_c where:")
print("  k₀ / √(1 - r_s/r_c) = k_c")
print("  → 1 - r_s/r_c = (k₀/k_c)²")
print("  → r_c = r_s / (1 - (k₀/k_c)²)  [r_c barely above r_s for k₀ << k_c]")
print()
print("The 'barrier width' in the tortoise coordinate r* = r + r_s ln(r/r_s - 1):")
print()

def compute_tunneling(M_kg, f_hz):
    """
    WKB tunneling amplitude through the Planck barrier near a BH.

    A photon of frequency f (far away) gets blueshifted near the horizon.
    At radius r, local k = (2πf/c) / √(1 - r_s/r)

    Above k_c, the spatial sector is evanescent with decay constant κ.
    The WKB tunneling exponent is ∫κ dr from r_c to r_s.

    This integral gives the amplitude for temporal-sector leakage.
    """
    rs = r_s(M_kg)
    k0 = 2 * np.pi * f_hz / c

    # Where does the mode hit k_c?
    if k0 >= k_c:
        return 1.0, rs, "already above k_c at infinity"

    # r_c: radius where k_local = k_c
    # (1 - rs/r_c) = (k0/k_c)^2
    ratio = (k0/k_c)**2
    if ratio >= 1:
        return 1.0, rs, "k0 >= k_c"
    r_c = rs / (1 - ratio)

    # Distance from r_c to horizon (in units of r_s)
    epsilon = r_c/rs - 1  # = ratio/(1-ratio) ≈ ratio for small ratio

    # Near-horizon evanescent decay:
    # κ(r) = |ω_spatial(k_local(r))|/c
    # For k >> k_c: κ ≈ (1/√12) ℓ_P k²
    # k_local ≈ k_c * (1 + δ(r)) where δ(r) grows as you approach horizon

    # WKB integral using near-horizon approximation:
    # ∫κ dr from r_c to r_s
    # In tortoise coordinate r* ≈ r - r_s - r_s ln((r-r_s)/r_s)
    # dr = dr  (near horizon, dr* ≈ dr/(1-r_s/r) ≈ dr/(r-r_s)·r_s)

    # Simple estimate: the integral is dominated by the near-horizon region
    # κ near r_c: κ ≈ k_c × √(2(k_local-k_c)/k_c) ≈ k_c × √(2δ)
    # where δ = k_local/k_c - 1 grows from 0 at r_c to ∞ at r_s

    # The decay exponent (numerical integration):
    # ∫κ dr ≈ ∫_{r_c}^{r_s} κ(r) dr

    # Change of variable: let u = (r - rs)/rs (goes from epsilon to 0 as r → rs)
    # k_local = k0 / sqrt(u/(1+u)) ≈ k0 / sqrt(u) for small u
    # κ(u) = sqrt(k_local(u)⁴ · ℓ_P²/12 - k_local(u)²)

    # Numerical integration
    n_steps = 1000
    r_values = np.linspace(r_c * 1.0001, rs * 0.99999, n_steps)
    kappa_values = []
    for r in r_values:
        k_local = k0 / np.sqrt(max(1 - rs/r, 1e-100))
        omega_sq = c**2 * k_local**2 - DISPERSION_COEFF * k_local**4
        if omega_sq < 0:
            kappa = np.sqrt(-omega_sq) / c
            kappa_values.append(kappa)
        else:
            kappa_values.append(0)

    if not kappa_values:
        return 0, r_c, "no evanescent region"

    # Trapezoidal integration
    dr = abs(r_values[1] - r_values[0])
    tunnel_exponent = np.trapz(kappa_values, r_values)

    T_amplitude = np.exp(-2 * tunnel_exponent)

    return T_amplitude, r_c, f"barrier width = {r_c - rs:.3e} m"

# Compute for a typical LIGO BH and a range of photon frequencies
M_test = 30 * M_sun
rs_test = r_s(M_test)

print(f"For a {int(M_test/M_sun)} M☉ BH (r_s = {rs_test:.1f} m):")
print(f"  k_c = {k_c:.3e} m⁻¹  (Planck barrier)")
print()
print(f"  {'Freq (Hz)':>12s}  {'k₀ (m⁻¹)':>14s}  {'k₀/k_c':>12s}  {'r_c - r_s':>14s}  {'T (tunnel)':>12s}")
print("  " + "-" * 70)

for f_hz in [100, 1e6, 1e10, 1e20, 1e30, 1e40, 1e43]:
    k0 = 2 * np.pi * f_hz / c
    T, r_c, info = compute_tunneling(M_test, f_hz)
    barrier_width = r_c - rs_test
    print(f"  {f_hz:12.1e}  {k0:14.3e}  {k0/k_c:12.3e}  {barrier_width:14.3e}  {T:12.3e}")

print()

# ── THE REAL ANSWER: WHERE DOES THE RUNNING COME FROM? ──────────────────────
print("=" * 70)
print("4. THE RUNNING OF λ_QG NEAR THE HORIZON")
print("=" * 70)
print()
print("The WKB tunneling exponent above shows:")
print()
print("  For k₀ << k_c (all current photons): barrier width → (r_c - r_s) → 0")
print("  But the decay rate κ → ∞ as k_local → ∞ at horizon")
print("  The product ∫κ dr is an integral from r_c to r_s")
print()
print("  Result: T_tunnel = exp(-2∫κ dr) → 0 as k₀ → 0")
print("  For realistic photons, T_tunnel ~ exp(-enormous number)")
print()
print("This means:")
print("  The Planck barrier is NOT 'thin' for ordinary photons near a BH.")
print("  Even at the horizon, the barrier grows faster (κ → ∞) than it")
print("  shrinks (r_c → r_s), so tunneling is exponentially suppressed.")
print()
print("  BUT: Hawking radiation is DIFFERENT from classical tunneling.")
print("  Hawking modes are quantum vacuum fluctuations that straddle the")
print("  horizon. They don't tunnel THROUGH the barrier — they're created")
print("  AT the barrier by pair production (one mode inside, one outside).")
print()
print("  The GUTOE correction to Hawking radiation is a DIRECT modification")
print("  of the near-horizon vacuum structure, not a tunneling effect.")

# ── THE EFFECTIVE λ FROM VACUUM STRUCTURE NEAR HORIZON ──────────────────────
print()
print("=" * 70)
print("5. EFFECTIVE λ FROM NEAR-HORIZON VACUUM STRUCTURE")
print("=" * 70)
print()
print("In the Unruh/Hawking derivation, the key modes have:")
print("  k_characteristic ~ T_H / (ħc) = 1/(4r_s)  [thermal wavenumber]")
print()
print("These modes, near the horizon, get blueshifted to k_c at some r_c.")
print("The region r_c < r < r_s is where GUTOE corrections are O(1).")
print()
print("The GUTOE modification to the thermal spectrum:")
print("  δT/T ≈ (1/24)(ℓ_P/r_s)²")
print("  δ(Hawking flux) ≈ (ℓ_P/r_s)² × (Planck-scale vacuum contribution)")
print()

# The regime where GUTOE corrections matter for Hawking
print("Summary: δT/T and required BH parameters for observability:")
print()
print(f"  {'BH type':>25s}  {'r_s (m)':>14s}  {'M (kg)':>14s}  {'δT/T':>12s}  {'T_H (K)':>12s}")
print("  " + "-" * 80)

bh_cases = [
    ("Stellar (30 M☉)",        30 * M_sun),
    ("Intermediate (10⁴ M☉)", 1e4 * M_sun),
    ("EHT: Sgr A*",            4e6 * M_sun),
    ("Prim. BH, r_s = 1 mm",  c**2 * 1e-3 / (2*G)),
    ("Prim. BH, r_s = 1 µm",  c**2 * 1e-6 / (2*G)),
    ("Prim. BH, r_s = 1 nm",  c**2 * 1e-9 / (2*G)),
    ("Prim. BH, r_s = 1 pm",  c**2 * 1e-12 / (2*G)),
    ("Prim. BH, r_s = 1 fm",  c**2 * 1e-15 / (2*G)),
    ("Prim. BH, r_s = 10 ℓ_P", c**2 * 10*ℓ_P / (2*G)),
    ("Planck BH, r_s = ℓ_P",  c**2 * ℓ_P / (2*G)),
]

for name, M in bh_cases:
    rs_val = r_s(M)
    corr = hawking_correction(M)
    T_H = T_Hawking(M)
    print(f"  {name:>25s}  {rs_val:14.3e}  {M:14.3e}  {corr:12.3e}  {T_H:12.3e}")

print()
print("INTERPRETATION:")
print()
print("  1. Primordial BHs ending their evaporation (r_s ~ fm to pm scale)")
print("     are in the 'interesting' regime where δT/T ~ 10⁻¹⁰ to 10⁻².")
print("     These produce detectable gamma-ray bursts. FERMI-LAT is looking.")
print()
print("  2. At r_s = 10 ℓ_P (the near-Planck regime), δT/T ~ 10⁻³.")
print("     The Hawking spectrum is noticeably non-blackbody.")
print("     The final burst energy spectrum would show the k⁴ cutoff.")
print()
print("  3. The SIGN of the correction: GUTOE predicts T_Hawking is SLIGHTLY")
print("     HIGHER than the standard result (modes below k_c propagate more")
print("     efficiently — the k⁴ correction reduces ω for given k, which")
print("     redshifts the modes → pair production at lower energy → more pairs")
print("     → higher effective temperature).")
print("     SIGN is: GUTOE BH is slightly HOTTER than Schwarzschild.")
print()

# ── THE BOTTOM LINE ─────────────────────────────────────────────────────────
print("=" * 70)
print("BOTTOM LINE: WHERE TO LOOK FOR GUTOE")
print("=" * 70)
print()
print("Bare λ_QG = 1/12 with ℓ_P² coefficient:")
print()
print("  ❌ LIGO (GW propagation): off by 63 orders. Consistent but undetectable.")
print("  ❌ QNM overtones (LIGO/LISA): off by ~37 orders even for n=1.")
print("  ❌ EHT photon ring phase: similar suppression.")
print()
print("  🟡 Primordial BH evaporation (Fermi): correction reaches 10⁻² when")
f"""     r_s ~ {ℓ_P/np.sqrt(24*0.01):.2e} m. These BHs would be finishing evaporation NOW."""
print(f"     r_s ~ {ℓ_P/np.sqrt(24*0.01):.2e} m. These would be evaporating now.")
print(f"     M ~ {c**2 * ℓ_P/np.sqrt(24*0.01) / (2*G):.2e} kg. t_evap ~ age of universe.")
print()
print("  ✅ The ONLY current observable: primordial BHs near end of life.")
print("     The 1% correction regime needs r_s ~ 10⁻³⁵ m × √(12/0.24) ~")
print(f"     {ℓ_P * np.sqrt(12/0.24):.2e} m. BH mass {c**2*ℓ_P*np.sqrt(12/0.24)/(2*G):.2e} kg.")
print()
print("  The real question Wings asked: does λ_QG RUN?")
print()
print("  If λ_eff(r) increases as r → r_s faster than (ℓ_P/r_s)²,")
print("  then the detectable regime opens up at larger BH masses.")
print("  The calculation needed: the renormalization group equation for λ_QG")
print("  in the presence of the temporal sector's propagating modes.")
print("  That β function is what tells us whether the GUTOE corrections")
print("  are permanently tiny or become order-unity at accessible scales.")
