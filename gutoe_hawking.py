#!/usr/bin/env python3
"""
GUTOE Hawking radiation: the focused calculation Wings asked for.

Three questions:
  1. At what radius does local k enter the GUTOE-significant regime?
  2. How much of the Hawking thermal spectrum originates from
     modes that passed through that regime?
  3. Does modified dispersion change the TEMPERATURE or just the tail?

Answer to (3) first: GUTOE has a maximum propagating frequency
    ω_max = c²√3 / ℓ_P  (in Planck units: ω_max = √3)
The Hawking spectrum is a Planck distribution cut off at ω_max.
For T_H << T_Planck: cutoff is far above spectrum peak → temperature unchanged.
For T_H ~ T_Planck: cutoff truncates the spectrum → total power reduced.

This is the 63 orders back from black holes as amplifiers.
"""

import numpy as np
from scipy import integrate

print("=" * 70)
print("GUTOE HAWKING RADIATION: THE FOCUSED CALCULATION")
print("=" * 70)
print()

# ── Constants (natural Planck units where ħ = c = G = 1 for key calcs) ──────
c     = 2.998e8          # m/s
G     = 6.674e-11        # m³ kg⁻¹ s⁻²
ħ     = 1.055e-34        # J·s
k_B   = 1.381e-23        # J/K
ℓ_P   = 1.616e-35        # m
M_P   = np.sqrt(ħ * c / G)   # Planck mass = 2.18e-8 kg
T_P   = np.sqrt(ħ * c**5 / (G * k_B**2))  # Planck temperature

lambda_QG = 1/12

# The maximum propagating frequency in GUTOE:
# ω²(k) = c²k² - (1/12)ℓ_P²k⁴, maximum at k = k_c/√2
# ω_max = c × k_c / √(4/3) = c²√3 / ℓ_P
k_c   = c * np.sqrt(12) / ℓ_P           # Planck cutoff wavenumber
omega_max = c**2 * np.sqrt(3) / ℓ_P     # Maximum propagating frequency
T_max = ħ * omega_max / k_B              # Temperature equivalent of ω_max

print(f"λ_QG = 1/12")
print(f"k_c  = {k_c:.3e} m⁻¹")
print(f"ω_max = {omega_max:.3e} rad/s  (max propagating frequency)")
print(f"k_BT_max = ħω_max = {ħ*omega_max:.3e} J = {T_max:.3e} K")
print(f"T_Planck  = {T_P:.3e} K")
print(f"T_max / T_Planck = {T_max/T_P:.3f}  (≈ √3)")
print()

# ── Question 1: At what radius does local k become significant? ───────────────
print("=" * 70)
print("QUESTION 1: At what radius does local k enter the GUTOE regime?")
print("=" * 70)
print()

def r_s(M_kg):
    return 2 * G * M_kg / c**2

def T_Hawking(M_kg):
    return ħ * c**3 / (8 * np.pi * G * M_kg * k_B)

# For a Hawking photon with frequency ω (at infinity), the local frequency
# at radius r (outside the horizon) is:
#   ω_local(r) = ω / √(1 - r_s/r)
# The local wavenumber: k_local = ω_local / c

# GUTOE regime: k_local > k_c/10  (10% of Planck cutoff, ~1% correction)
# The radius where this happens for a photon at the peak Hawking frequency:
#   ω_peak ≈ 2.82 k_B T_H / ħ  (Wien's law)

def radius_of_gutoe_regime(M_kg, fraction_of_kc=0.1):
    """
    Radius where locally blueshifted Hawking photons enter the regime
    where k_local > fraction × k_c (i.e., GUTOE corrections become noticeable).
    """
    rs = r_s(M_kg)
    T_H = T_Hawking(M_kg)
    omega_peak = 2.82 * k_B * T_H / ħ  # Wien's law peak
    k_peak = omega_peak / c             # far-away wavenumber of peak photon

    k_threshold = fraction_of_kc * k_c

    if k_peak >= k_threshold:
        return rs, "Peak photon already above threshold at infinity!"

    # k_local = k_peak / sqrt(1 - rs/r) = k_threshold
    # 1 - rs/r = (k_peak/k_threshold)^2
    ratio = (k_peak / k_threshold)**2
    if ratio >= 1:
        return rs, "Already above threshold everywhere"

    r_c = rs / (1 - ratio)
    return r_c, f"r_c = {r_c:.3e} m = {r_c/rs:.6f} × r_s"

print("Radius where Hawking peak photons enter GUTOE-significant regime")
print("(defined as k_local > k_c/10, i.e., correction ~ 1%):")
print()
print(f"  {'BH mass':>20s}  {'r_s':>12s}  {'r_c':>14s}  {'(r_c-r_s)/r_s':>16s}  {'T_H':>12s}")
print("  " + "-" * 80)

for M, name in [
    (30 * 1.989e30, "30 M☉ (LIGO)"),
    (1e5 * 1.989e30, "10⁵ M☉ (ISCO)"),
    (4e6 * 1.989e30, "Sgr A*"),
    (1e8, "PBH 10⁸ kg"),
    (1e6, "PBH 10⁶ kg"),
    (1e4, "PBH 10⁴ kg"),
    (1e2, "PBH 100 kg"),
    (1.0, "PBH 1 kg"),
    (1e-4, "PBH 0.1 g"),
    (M_P * 10, "10 M_Planck"),
    (M_P, "1 M_Planck"),
    (M_P / 10, "0.1 M_Planck"),
]:
    rs = r_s(M)
    T_H = T_Hawking(M)
    r_c, info = radius_of_gutoe_regime(M)
    if isinstance(r_c, float):
        dr = (r_c - rs) / rs
        print(f"  {name:>20s}  {rs:12.3e}  {r_c:14.3e}  {dr:16.3e}  {T_H:12.3e} K")
    else:
        print(f"  {name:>20s}  {rs:12.3e}  {'AT HORIZON':>14s}  {'---':>16s}  {T_H:12.3e} K")

print()
print("The 'GUTOE-significant radius' r_c is exponentially close to r_s for all")
print("BHs until M ~ M_Planck. The Planck barrier is crossed essentially AT the horizon.")
print()

# ── Question 2: How much of the spectrum comes from trans-Planckian modes? ───
print("=" * 70)
print("QUESTION 2: What fraction of the Hawking spectrum is affected?")
print("=" * 70)
print()
print("The Hawking spectrum: P(ω) ∝ ω² / (exp(ħω/k_BT_H) - 1)")
print()
print("GUTOE modifies modes with ω > ω_max = √3 × T_Planck/ħk_B.")
print("The fraction of thermal energy in modes ω > ω_max:")
print("  f_affected = ∫_{ω_max}^∞ P(ω)dω / ∫_0^∞ P(ω)dω")
print()

def planck_integrand(x):
    """x³/(exp(x)-1) — the Stefan-Boltzmann integrand"""
    if x > 500:
        return 0
    return x**3 / (np.exp(x) - 1)

stefan_boltzmann_integral = np.pi**4 / 15  # = ∫₀^∞ x³/(e^x-1) dx

def fraction_above_cutoff(T_H):
    """Fraction of Hawking energy in modes that exceed ω_max."""
    x_max = ħ * omega_max / (k_B * T_H)
    if x_max > 500:
        return 0.0  # negligible
    # ∫_{x_max}^∞ x³/(exp(x)-1) dx
    result, _ = integrate.quad(planck_integrand, x_max, 500)
    return result / stefan_boltzmann_integral

def fraction_below_cutoff(T_H):
    """Fraction of Hawking energy in propagating modes (below ω_max)."""
    return 1 - fraction_above_cutoff(T_H)

def gutoe_power_ratio(T_H):
    """
    Ratio of GUTOE Hawking power to standard Hawking power.
    In GUTOE, modes above ω_max are either:
      (a) Suppressed entirely (evanescent spatial) → ratio < 1
      (b) Redirected to temporal sector → different outgoing channel

    The conserved quantity is energy. If temporal-sector modes contribute
    to the BH entropy but NOT to the escaping radiation, the effective
    evaporation rate is reduced.

    Here we compute the fraction assuming ALL modes above ω_max are blocked.
    This is a lower bound — some temporal-sector energy may contribute.
    """
    return fraction_below_cutoff(T_H)

print(f"  {'T_H / T_Planck':>15s}  {'x_max = ħω_max/k_BT_H':>24s}  {'Power ratio':>14s}  {'Correction':>12s}")
print("  " + "-" * 70)

for t_ratio in [1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.1, 0.3, 1.0, 3.0, 10.0]:
    T_H = t_ratio * T_P
    x_max = ħ * omega_max / (k_B * T_H)
    ratio = gutoe_power_ratio(T_H)
    correction = 1 - ratio
    print(f"  {t_ratio:15.2e}  {x_max:24.3e}  {ratio:14.6f}  {correction:12.3e}")

print()

# At what T_H does the correction reach 10%?
# We need fraction_above_cutoff(T_H) = 0.10
from scipy.optimize import brentq

def find_threshold_temperature(target_correction):
    def f(log_T_ratio):
        T_H = 10**log_T_ratio * T_P
        return fraction_above_cutoff(T_H) - target_correction

    try:
        log_ratio = brentq(f, -20, 2)
        return 10**log_ratio
    except:
        return None

print("Temperature (as T/T_Planck) where GUTOE power reduction reaches:")
for threshold in [0.001, 0.01, 0.1, 0.5]:
    T_ratio = find_threshold_temperature(threshold)
    if T_ratio:
        M = ħ * c**3 / (8 * np.pi * G * k_B * T_ratio * T_P)
        rs = r_s(M)
        t_evap = 5120 * np.pi * G**2 * M**3 / (ħ * c**4)
        print(f"  {threshold*100:.1f}%: T/T_P = {T_ratio:.3e}, M = {M:.3e} kg = {M/M_P:.2f} M_P,"
              f"  t_evap = {t_evap:.2e} s")

print()

# ── Question 3: Temperature change or spectral shape change? ─────────────────
print("=" * 70)
print("QUESTION 3: Does it change the temperature or the high-freq tail?")
print("=" * 70)
print()
print("Answer: BOTH, but in different regimes.")
print()
print("For T_H << T_Planck (all macroscopic BHs):")
print("  x_max = ħω_max/k_BT_H >> 1 (cutoff is far above peak)")
print("  Only modes in the exponentially suppressed Wien tail are affected.")
print("  → Pure high-frequency tail modification. Temperature UNCHANGED.")
print()
print("For T_H ~ T_Planck (near-Planck BHs, M ~ M_Planck):")
print("  x_max = ħω_max/k_BT_H ~ 1 (cutoff is near the spectral peak)")
print("  Significant fraction of modes are suppressed.")
print("  → TOTAL POWER REDUCED, evaporation slows down.")
print("  → The effective temperature as SEEN FROM INFINITY is reduced.")
print("  → The BH appears to be at a lower temperature than T_H predicts.")
print()
print("The effective temperature T_eff (fitting a blackbody to GUTOE spectrum):")

def effective_temperature(T_H):
    """
    Fit an effective temperature to the GUTOE-truncated Planck spectrum.
    The GUTOE spectrum is Planck(T_H) truncated at ω_max.
    We find T_eff such that the total power matches Stefan-Boltzmann(T_eff).
    P_GUTOE = σ T_eff⁴  →  T_eff = T_H × (fraction_below_cutoff)^(1/4)
    """
    frac = fraction_below_cutoff(T_H)
    return T_H * frac**0.25

print()
print(f"  {'T_H/T_Planck':>15s}  {'T_eff/T_Planck':>17s}  {'T_eff/T_H':>12s}  {'Interpretation':>30s}")
print("  " + "-" * 80)

for t_ratio in [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]:
    T_H = t_ratio * T_P
    T_eff = effective_temperature(T_H)
    ratio = T_eff / T_H
    if ratio > 0.99:
        interp = "Standard Hawking"
    elif ratio > 0.9:
        interp = "Small correction"
    elif ratio > 0.5:
        interp = "Significant modification"
    else:
        interp = "STRONGLY modified"
    print(f"  {t_ratio:15.3f}  {T_eff/T_P:17.4f}  {ratio:12.4f}  {interp:>30s}")

print()
print("CONCLUSION on temperature:")
print("  For T_H/T_Planck < 0.1 (M > 10 M_Planck): T_eff ≈ T_H (unchanged)")
print("  For T_H/T_Planck ~ 1 (M ~ M_Planck): T_eff < T_H (BH appears cooler)")
print("  The 'cooling' is physical: modes above ω_max don't escape spatially.")
print()

# ── The bottom line: primordial BHs and Fermi ─────────────────────────────────
print("=" * 70)
print("THE FERMI CONNECTION: PRIMORDIAL BHs EVAPORATING NOW")
print("=" * 70)
print()

# BH that finishes evaporating in the age of universe
t_universe = 4.35e17  # seconds
M_now = (t_universe * ħ * c**4 / (5120 * np.pi * G**2))**(1/3)
rs_now = r_s(M_now)
T_H_now = T_Hawking(M_now)
correction_now = 1 - gutoe_power_ratio(T_H_now)

print(f"Primordial BH completing evaporation NOW:")
print(f"  M = {M_now:.3e} kg  =  {M_now/M_P:.3e} M_Planck")
print(f"  r_s = {rs_now:.3e} m  =  {rs_now/ℓ_P:.3e} ℓ_P")
print(f"  T_H = {T_H_now:.3e} K  =  {T_H_now/T_P:.3e} T_Planck")
print(f"  GUTOE power correction: {correction_now:.3e}")
print()
print(f"  → Correction is {correction_now:.3e} — still negligible for Fermi.")
print()
print("Where is the correction 1% for Fermi?")
T_ratio_1pct = find_threshold_temperature(0.01)
if T_ratio_1pct:
    M_1pct = ħ * c**3 / (8 * np.pi * G * k_B * T_ratio_1pct * T_P)
    rs_1pct = r_s(M_1pct)
    t_evap_1pct = 5120 * np.pi * G**2 * M_1pct**3 / (ħ * c**4)
    print(f"  1% correction at M = {M_1pct:.3e} kg = {M_1pct/M_P:.3f} M_Planck")
    print(f"  r_s = {rs_1pct:.3e} m = {rs_1pct/ℓ_P:.3f} ℓ_P")
    print(f"  Evaporation time from that mass: t_evap = {t_evap_1pct:.3e} s")
    print()
    print(f"  These BHs are near-Planck. We can't detect their final burst.")
    print(f"  The burst energy: E = M_1pct × c² = {M_1pct * c**2:.3e} J = {M_1pct*c**2/1.602e-13:.3e} MeV")

print()
print("=" * 70)
print("THE ACTUAL ANSWER WINGS' QUESTION FORCES:")
print("=" * 70)
print()
print("The GUTOE dispersion correction to Hawking radiation:")
print()
print("  Is real, is computable, is well-defined.")
print("  Is negligible for ALL currently observable BHs.")
print("  Becomes order-unity ONLY at M < M_Planck.")
print()
print("  The 63 orders are NOT recovered by black holes as amplifiers —")
print("  not with the bare λ_QG = 1/12 coefficient.")
print()
print("  THE ONLY ESCAPE ROUTE: running coupling λ_QG(r).")
print()
print("  If λ_eff(r) grows near the horizon faster than the naive")
print("  (ℓ_P/r_s)² suppression, the detectable regime opens up.")
print()
print("  Specifically: if λ_eff(r) = λ_QG × (r_s/r - r_s/r_c)^{-α}")
print("  for some α > 0, then there's a mass scale where corrections")
print("  become observable before you hit M_Planck.")
print()
print("  That running is NOT in the classical dispersion relation.")
print("  It would come from the RG flow in the presence of the")
print("  temporal sector — which is what we proved propagates above k_c.")
print()
print("  The β function for λ_QG in GUTOE = the next calculation.")
