#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LEAN_DIR="$ROOT_DIR/lean"
OUT_DIR="$ROOT_DIR/findings/assets/clay/submission"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_FILE="$OUT_DIR/submission_repro_${STAMP}.log"
THEOREM_FILE="$OUT_DIR/submission_theorem_presence_${STAMP}.txt"

mkdir -p "$OUT_DIR"

MODULES=(
  "Gutoe.YangMillsContinuumLimit"
  "Gutoe.YangMillsOSEndToEnd"
  "Gutoe.YangMillsWilsonEquivalence"
  "Gutoe.YangMillsOSCompletion"
  "Gutoe"
)

THEOREMS=(
  "lean/Gutoe/YangMillsContinuumLimit.lean:constructive_schwinger_family_exists"
  "lean/Gutoe/YangMillsOSEndToEnd.lean:grand331_end_to_end_os_reconstruction_of_domain"
  "lean/Gutoe/YangMillsWilsonEquivalence.lean:theorem_c_wilson_equivalence_domain_limits"
  "lean/Gutoe/YangMillsOSCompletion.lean:osGenerator_uniform_gap_floor_of_domain"
)

FILES_WITHOUT_SORRY=(
  "lean/Gutoe/YangMillsContinuumLimit.lean"
  "lean/Gutoe/YangMillsOSEndToEnd.lean"
  "lean/Gutoe/YangMillsWilsonEquivalence.lean"
  "lean/Gutoe/YangMillsOSCompletion.lean"
)

{
  echo "Clay submission reproducibility bundle"
  echo "timestamp_utc: $STAMP"
  echo "root: $ROOT_DIR"
  echo
  echo "[1/3] Lean builds"
  for mod in "${MODULES[@]}"; do
    echo ">>> lake build $mod"
    (cd "$LEAN_DIR" && lake build "$mod")
  done
  echo
  echo "[2/3] Theorem presence checks"
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
  echo
  echo "[3/3] No-sorry checks (target files)"
  for file in "${FILES_WITHOUT_SORRY[@]}"; do
    if rg -n "\\bsorry\\b" "$ROOT_DIR/$file" \
      | rg -v "^\\s*(--|/-|\\*|\\*/)" \
      | rg -F -v '`sorry`' >/dev/null; then
      echo "FOUND_SORRY: $file"
      rg -n "\\bsorry\\b" "$ROOT_DIR/$file" \
        | rg -v "^\\s*(--|/-|\\*|\\*/)" \
        | rg -F -v '`sorry`'
      exit 1
    else
      echo "OK: $file (no sorry)"
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
