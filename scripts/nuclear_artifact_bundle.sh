#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_BASE="${1:-/tmp/nuclear_chart}"
BUNDLE="${2:-/tmp/nuclear_chart_bundle}"

mkdir -p "$OUT_BASE" "$BUNDLE"

cd "$ROOT"

echo "[1/4] mass_periodic_report"
GUTOE_MASS_PERIODIC_OUT="$OUT_BASE" \
cargo run -q -p gutoe-physics --bin mass_periodic_report

echo "[2/4] AME2020 benchmark"
GUTOE_MASS_PERIODIC_OUT="$OUT_BASE" \
cargo run -q -p gutoe-physics --bin ame2020_benchmark

echo "[3/4] residual structure plot"
python3 "$ROOT/scripts/ame2020_residual_plot.py" \
  --csv "$OUT_BASE/ame2020_residuals.csv" \
  --out "$OUT_BASE/ame2020_residual_structure.png"

echo "[4/4] assemble bundle"
cp "$OUT_BASE/mass_periodic_report.json" "$BUNDLE/"
cp "$OUT_BASE/periodic_table_scoreboard.csv" "$BUNDLE/"
cp "$OUT_BASE/shell_gap_attenuation.csv" "$BUNDLE/"
cp "$OUT_BASE/tin_isotope_diagnostics.csv" "$BUNDLE/"
cp "$OUT_BASE/superheavy_closure_derivation.csv" "$BUNDLE/"
cp "$OUT_BASE/ame2020_benchmark.json" "$BUNDLE/"
cp "$OUT_BASE/ame2020_residuals.csv" "$BUNDLE/"
cp "$OUT_BASE/ame2020_residuals_top50.csv" "$BUNDLE/"
cp "$OUT_BASE/ame2020_residual_structure.png" "$BUNDLE/"

cat > "$BUNDLE/README.txt" <<'EOF'
GUTOE Nuclear Artifact Bundle

Contents:
- mass_periodic_report.json: periodic table + shell summary
- periodic_table_scoreboard.csv: per-Z stable-isotope counts
- shell_gap_attenuation.csv: magic-gap ratio diagnostics
- tin_isotope_diagnostics.csv: Sn gate-by-gate diagnostics
- superheavy_closure_derivation.csv: GRAND-263 closure scoring
- ame2020_benchmark.json: AME2020 RMS/MAE/bias metrics
- ame2020_residuals.csv: full residual table (pred-obs binding MeV)
- ame2020_residuals_top50.csv: largest absolute residuals
- ame2020_residual_structure.png: residual trends vs N and A
EOF

echo "Bundle ready: $BUNDLE"
