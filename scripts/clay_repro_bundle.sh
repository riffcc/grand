#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LEAN_DIR="$ROOT_DIR/lean"
OUT_DIR="$ROOT_DIR/findings/assets/clay"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_FILE="$OUT_DIR/repro_${STAMP}.log"
THEOREM_FILE="$OUT_DIR/theorem_presence_${STAMP}.txt"

mkdir -p "$OUT_DIR"

MODULES=(
  "Gutoe.YangMillsStructuralGap"
  "Gutoe.YangMillsMassGap"
  "Gutoe.YangMillsContinuumSurvival"
  "Gutoe.HaarBridgeScaffold"
  "Gutoe.HaarMeasureHooks"
  "Gutoe.HaarExpectationDecomposition"
  "Gutoe.HaarFiberCollapse"
  "Gutoe.YangMillsWilsonBridge"
  "Gutoe.YangMillsWilsonEquivalence"
  "Gutoe.YangMillsFullStateLift"
  "Gutoe.YangMillsConstructiveQFT"
  "Gutoe.YangMillsGaugeScope"
  "Gutoe"
)

THEOREMS=(
  "lean/Gutoe/YangMillsStructuralGap.lean:doeblin_decomposition"
  "lean/Gutoe/YangMillsMassGap.lean:mass_gap_positive_of_doeblin_mode"
  "lean/Gutoe/YangMillsContinuumSurvival.lean:continuum_survival_gap_nonvanishing_of_z3_nn_schedule"
  "lean/Gutoe/HaarFiberCollapse.lean:normalized_expectation_reduce_to_center_of_quotient_normalization"
  "lean/Gutoe/YangMillsWilsonBridge.lean:c1_z3_to_su3_structural_construction"
  "lean/Gutoe/YangMillsWilsonBridge.lean:c2_counts_center_action_bijective"
  "lean/Gutoe/YangMillsWilsonBridge.lean:c3_wilson_gap_nonvanishing_from_clifford_z3"
  "lean/Gutoe/YangMillsWilsonEquivalence.lean:theorem_c_wilson_equivalence_domain_limits"
  "lean/Gutoe/YangMillsFullStateLift.lean:full_gap_positive_all_steps_of_wilson_center_schedule"
  "lean/Gutoe/YangMillsConstructiveQFT.lean:constructive_lane_gap_closure_of_wilson_equivalence_domain"
  "lean/Gutoe/YangMillsGaugeScope.lean:compact_simple_scope_supports_full_path2"
)

{
  echo "Clay reproducibility bundle"
  echo "timestamp_utc: $STAMP"
  echo "root: $ROOT_DIR"
  echo
  echo "[1/2] Lean builds"
  for mod in "${MODULES[@]}"; do
    echo ">>> lake build $mod"
    (cd "$LEAN_DIR" && lake build "$mod")
  done
  echo
  echo "[2/2] Theorem presence checks"
  for entry in "${THEOREMS[@]}"; do
    file="${entry%%:*}"
    name="${entry##*:}"
    if rg -n "^theorem ${name}\\b|^def ${name}\\b|^structure ${name}\\b" "$ROOT_DIR/$file" >/dev/null; then
      match="$(rg -n "^theorem ${name}\\b|^def ${name}\\b|^structure ${name}\\b" "$ROOT_DIR/$file" | head -n 1)"
      echo "OK: $file :: $match"
    else
      echo "MISSING: $file :: $name"
      exit 1
    fi
  done
} | tee "$LOG_FILE"

{
  echo "timestamp_utc: $STAMP"
  for entry in "${THEOREMS[@]}"; do
    file="${entry%%:*}"
    name="${entry##*:}"
    rg -n "^theorem ${name}\\b|^def ${name}\\b|^structure ${name}\\b" "$ROOT_DIR/$file" | head -n 1
  done
} > "$THEOREM_FILE"

echo
echo "Wrote:"
echo "  $LOG_FILE"
echo "  $THEOREM_FILE"
