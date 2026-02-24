#!/usr/bin/env bash
set -euo pipefail

# ROCm/HIP smoke for GRAND GPU stack.
# Run this on the AMD box (e.g. 10.7.1.195 / tealc).

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[rocm] root: $ROOT_DIR"
echo "[rocm] rust: $(rustc --version)"

export HIPCC="${HIPCC:-/opt/rocm/bin/hipcc}"
export ROCM_PATH="${ROCM_PATH:-/opt/rocm}"
export GFX_ARCH="${GFX_ARCH:-gfx1100}"
if [[ -z "${HIP_WAVEFRONT_SIZE:-}" ]]; then
  if [[ "${GFX_ARCH}" == gfx11* ]]; then
    export HIP_WAVEFRONT_SIZE=32
  else
    export HIP_WAVEFRONT_SIZE=64
  fi
else
  export HIP_WAVEFRONT_SIZE
fi
export HIP_LIB_DIR="${HIP_LIB_DIR:-}"

echo "[rocm] HIPCC=$HIPCC"
echo "[rocm] ROCM_PATH=$ROCM_PATH"
echo "[rocm] GFX_ARCH=$GFX_ARCH"
echo "[rocm] HIP_WAVEFRONT_SIZE=$HIP_WAVEFRONT_SIZE"

if [[ ! -x "$HIPCC" ]]; then
  echo "[rocm] ERROR: hipcc not found at $HIPCC"
  exit 1
fi

if [[ -z "$HIP_LIB_DIR" ]]; then
  for d in "$ROCM_PATH/lib" "$ROCM_PATH/lib64" /opt/rocm/lib /opt/rocm/lib64; do
    if [[ -e "$d/libamdhip64.so" ]] || ls "$d"/libamdhip64.so.* >/dev/null 2>&1; then
      HIP_LIB_DIR="$d"
      break
    fi
  done
fi
if [[ -z "$HIP_LIB_DIR" ]]; then
  echo "[rocm] ERROR: could not locate libamdhip64.so (set HIP_LIB_DIR explicitly)"
  exit 1
fi
export LD_LIBRARY_PATH="$HIP_LIB_DIR:${LD_LIBRARY_PATH:-}"
echo "[rocm] HIP_LIB_DIR=$HIP_LIB_DIR"
echo "[rocm] LD_LIBRARY_PATH=$LD_LIBRARY_PATH"

echo "[rocm] 1/4 build gutoe-gpu with HIP"
cargo build -p gutoe-gpu --features rocm

echo "[rocm] 2/4 run GPU bohr smoke"
cargo test -p gutoe-gpu --features rocm -- --nocapture gpu_bohr_test

echo "[rocm] 3/4 run bohr convergence scan smoke"
cargo test -p gutoe-gpu --features rocm --release -- --nocapture bohr_convergence_scan

echo "[rocm] 4/4 run bh_render transfer parity (HIP feature build)"
BH_VALIDATE_GPU=1 cargo run -p gutoe-gpu --features rocm --bin bh_render -- m87star 1280x720

echo "[rocm] OK"
