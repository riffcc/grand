/- 
 * GUTOE — Lean Kerr Reference Frame (low-res oracle scaffold)
 * Copyright (C) 2026 Riff Labs
 * AGPL-3.0-or-later
 *
 * A tiny deterministic reference intensity field in Lean used as a sanity
 * oracle for renderer parity checks. This is intentionally low resolution and
 * algebraic-first (correctness over speed).
 *
 * No `sorry`.
-/

import Mathlib
import Gutoe.KerrTracerEquations

namespace Gutoe.KerrReferenceFrame

open Real
open Gutoe.KerrTracerEquations

/-- Observer-frame alpha coordinate in `[-fov, fov]` from integer pixel index. -/
noncomputable def alphaCoord (w x : ℕ) (fov : ℝ) : ℝ :=
  let wf : ℝ := w
  let xf : ℝ := x
  (2 * ((xf + 0.5) / wf) - 1) * fov

/-- Observer-frame beta coordinate in `[-fov, fov]` from integer pixel index. -/
noncomputable def betaCoord (h y : ℕ) (fov : ℝ) : ℝ :=
  let hf : ℝ := h
  let yf : ℝ := y
  (1 - 2 * ((yf + 0.5) / hf)) * fov

/-- Algebraic Kerr intensity proxy from radial+polar potential magnitudes. -/
noncomputable def referenceIntensity
    (rObs r_s aStar thetaObs alpha beta : ℝ) : ℝ :=
  let ξ := xi alpha thetaObs
  let a := Gutoe.KerrGeometry.spinLength r_s aStar
  let η := eta alpha beta a thetaObs
  let rad := max 0 (radialPotential rObs ξ η r_s aStar)
  let pol := max 0 (polarPotential thetaObs ξ η r_s aStar)
  let s := rad + pol
  s / (1 + s)

theorem referenceIntensity_nonneg
    (rObs r_s aStar thetaObs alpha beta : ℝ) :
    0 ≤ referenceIntensity rObs r_s aStar thetaObs alpha beta := by
  unfold referenceIntensity
  set ξ : ℝ := xi alpha thetaObs
  set a : ℝ := Gutoe.KerrGeometry.spinLength r_s aStar
  set η : ℝ := eta alpha beta a thetaObs
  set rad : ℝ := max 0 (radialPotential rObs ξ η r_s aStar)
  set pol : ℝ := max 0 (polarPotential thetaObs ξ η r_s aStar)
  have hrad : 0 ≤ rad := by
    simp [rad]
  have hpol : 0 ≤ pol := by
    simp [pol]
  have hs : 0 ≤ rad + pol := add_nonneg hrad hpol
  exact div_nonneg hs (by linarith)

theorem referenceIntensity_lt_one
    (rObs r_s aStar thetaObs alpha beta : ℝ) :
    referenceIntensity rObs r_s aStar thetaObs alpha beta < 1 := by
  unfold referenceIntensity
  set ξ : ℝ := xi alpha thetaObs
  set a : ℝ := Gutoe.KerrGeometry.spinLength r_s aStar
  set η : ℝ := eta alpha beta a thetaObs
  set rad : ℝ := max 0 (radialPotential rObs ξ η r_s aStar)
  set pol : ℝ := max 0 (polarPotential thetaObs ξ η r_s aStar)
  have hrad : 0 ≤ rad := by simp [rad]
  have hpol : 0 ≤ pol := by simp [pol]
  have hs : 0 ≤ rad + pol := add_nonneg hrad hpol
  have hden : 0 < 1 + (rad + pol) := by linarith
  have hnum_lt_den : rad + pol < 1 + (rad + pol) := by linarith
  exact (div_lt_one hden).2 hnum_lt_den

/-- Minimal Float reference pixel (`0..255`) for quick parity snapshots. -/
def referencePixelProxyFloat
    (w h x y : Nat)
    (fov rObs r_s aStar thetaObs : Float) : Float :=
  let wf := Float.ofNat w
  let hf := Float.ofNat h
  let α := (2.0 * ((Float.ofNat x + 0.5) / wf) - 1.0) * fov
  let β := (1.0 - 2.0 * ((Float.ofNat y + 0.5) / hf)) * fov
  let s := Float.sin thetaObs
  let c := Float.cos thetaObs
  let m := 0.5 * r_s
  let a := aStar * m
  let ξ := -α * s
  let η := β * β + (α * α - a * a) * (c * c)
  let Δ := rObs * rObs - r_s * rObs + a * a
  let t := (rObs * rObs + a * a) - a * ξ
  let rad := max 0.0 (t * t - Δ * ((ξ - a) * (ξ - a) + η))
  let sinT := max (Float.sin thetaObs) 1e-6
  let cot := (Float.cos thetaObs) / sinT
  let pol := max 0.0 (η + a * a * c * c - ξ * ξ * cot * cot)
  let sum := rad + pol
  let i := sum / (1.0 + sum)
  max 0.0 (min 255.0 (255.0 * i))

/-- Kerr horizon radius `r+` in our `r_s` normalization. -/
def horizonOuterFloat (r_s aStar : Float) : Float :=
  let m := 0.5 * r_s
  let a := aStar * m
  m + Float.sqrt (max 0.0 (m * m - a * a))

/-- One branch of Kerr geodesic integration; returns accumulated geodesic activity. -/
partial def tracePhiKerrBranch
    (rObs r_s aStar alpha beta thetaObs maxLambda dLambda sgnThInit : Float)
    (fuel : Nat) : Float :=
  let s := Float.sin thetaObs
  let c := Float.cos thetaObs
  let m := 0.5 * r_s
  let a := aStar * m
  let ξ := -alpha * s
  let η := beta * beta + (alpha * alpha - a * a) * (c * c)
  let rPlus := horizonOuterFloat r_s aStar
  let b := Float.sqrt (alpha * alpha + beta * beta)
  let rStart := max (40.0 * r_s) (max (12.0 * b) 20.0)
  let maxSteps := UInt64.toNat (Float.toUInt64 (Float.ceil (maxLambda / (max dLambda 1e-6)))) + 1
  let maxFuel := min fuel maxSteps

  let rec loop (i : Nat) (r th phi sgnR sgnTh lam acc : Float) : Float :=
    if i = 0 then
      0.0
    else
      let sinTh := Float.sin th
      let sin2 := max (sinTh * sinTh) 1e-9
      let sigma := max 1e-12 (r * r + a * a * (Float.cos th) * (Float.cos th))
      let delta := r * r - r_s * r + a * a
      let p := (r * r + a * a) - a * ξ
      let rpot := p * p - delta * ((ξ - a) * (ξ - a) + η)
      let tpot := η + a * a * (Float.cos th) * (Float.cos th) - ξ * ξ * (Float.cos th) * (Float.cos th) / sin2
      let sgnR' := if rpot <= 1e-12 then -sgnR else sgnR
      let sgnTh' := if tpot <= 1e-12 then -sgnTh else sgnTh
      let rr := Float.sqrt (max 0.0 rpot)
      let thh := Float.sqrt (max 0.0 tpot)
      let dr := sgnR' * rr / sigma
      let dth := sgnTh' * thh / sigma
      let dphi := if Float.abs delta < 1e-9 then
        (ξ / sin2 - a) / sigma
      else
        (ξ / sin2 - a + a * p / delta) / sigma

      let rNew := r + dLambda * dr
      let thNew := max 1e-4 (min (3.141592653589793 - 1e-4) (th + dLambda * dth))
      let phiNew := phi + dLambda * dphi
      let lamNew := lam + dLambda

      let accNew :=
        acc
          + Float.abs (dphi * dLambda)
          + 0.15 * Float.abs (dth * dLambda)
          + 0.02 * Float.abs (dr * dLambda / max r_s 1e-6)
      if !(rNew.isFinite && thNew.isFinite && phiNew.isFinite) then
        accNew
      else if rNew <= rPlus * 1.001 then
        accNew
      else if sgnR' > 0.0 && rNew >= rStart * 0.995 then
        accNew
      else if lamNew >= maxLambda then
        accNew
      else
        loop (i - 1) rNew thNew phiNew sgnR' sgnTh' lamNew accNew

  loop maxFuel rStart thetaObs 0.0 (-1.0) sgnThInit 0.0 0.0

/-- Kerr geodesic-based reference pixel (`0..255`), using the stronger escaped branch. -/
def referencePixelFloat
    (w h x y : Nat)
    (fov rObs r_s aStar thetaObs : Float) : Float :=
  let wf := Float.ofNat w
  let hf := Float.ofNat h
  let α := (2.0 * ((Float.ofNat x + 0.5) / wf) - 1.0) * fov
  let β := (1.0 - 2.0 * ((Float.ofNat y + 0.5) / hf)) * fov
  let maxLambda := 30.0 * 3.141592653589793
  let dLambda := 0.02
  let up := tracePhiKerrBranch rObs r_s aStar α β thetaObs maxLambda dLambda 1.0 200000
  let down := tracePhiKerrBranch rObs r_s aStar α β thetaObs maxLambda dLambda (-1.0) 200000
  -- Use both admissible polar branches for a smooth reference intensity.
  -- `max` introduces a hard branch boundary in screen space.
  let activity := up + down
  let p := max 0.0 (min 1.0 (activity / (2.0 * maxLambda)))
  255.0 * Float.sqrt p

end Gutoe.KerrReferenceFrame
