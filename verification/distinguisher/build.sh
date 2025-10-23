#!/bin/bash
# Copyright (C) 2025 Hosein Hadipour
# Email: hsn.hadipour@gmail.com
# Date: September 24, 2025

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# build.sh — simple CUDA build helper with sanitizers
# Usage:
#   ./build.sh <mode> [-- <program args>]
#
# Build modes:
#   release    (default) Configure & build CUDA + CPU targets → build/
#              Use when you want the GPU `cuda_dlct` binary. Requires NVCC.
#   debug      Same as release but with Debug flags.
#   openmp     Configure with CUDA disabled and build OpenMP targets only → build/
#              Pass args after `--` to run the CPU/OpenMP tool immediately.
#
# Maintenance/analysis modes:
#   clean      Remove build directories (build, build-asan, build-csan)
#   asan       Build with AddressSanitizer and run self-test → build-asan/
#   csan       Build CUDA binary and run NVIDIA Compute Sanitizer (memcheck + leak) → build-csan/
#              Pass program args after `--`, e.g. `./build.sh csan -- --self-test`
#
# Notes:
#   • `openmp` works on macOS/Linux without a GPU; install libomp for threaded runs on Apple Clang.
#   • CUDA modes (release/debug/csan) require nvcc and the CUDA toolkit; they can be built on Linux servers.
#   • ASan mode exercises code with AddressSanitizer for memory error detection.

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Build dirs (separate to avoid sanitizer cross-contamination)
BUILD_DIR_DEFAULT="$SCRIPT_DIR/build"
BUILD_DIR_ASAN="$SCRIPT_DIR/build-asan"
BUILD_DIR_CSAN="$SCRIPT_DIR/build-csan"

TARGET_NAME="cuda_dlct"

usage() {
  sed -n '2,35p' "$0" | sed 's/^# \{0,1\}//'
}

nproc_portable() {
  if command -v nproc >/dev/null 2>&1; then nproc
  elif command -v sysctl >/dev/null 2>&1; then sysctl -n hw.logicalcpu
  elif command -v getconf >/dev/null 2>&1; then getconf _NPROCESSORS_ONLN
  else echo 4; fi
}

BUILD_TYPE="Release"
MODE="build"

# Parse mode and optional program args (after -- for csan)
while [[ $# -gt 0 ]]; do
  case "$1" in
    clean)    MODE="clean" ;;
    openmp)   MODE="openmp" ;;
    debug)    BUILD_TYPE="Debug" ;;
    release)  BUILD_TYPE="Release" ;;
    asan)     MODE="asan" ;;
    csan)     MODE="csan" ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;; # remaining args go to the program in csan mode
    *) break ;;
  esac
  shift
done
PROGRAM_ARGS=("$@")

# ---- Helpers ----
host_includes() {
  local incs=()
  for d in "$SCRIPT_DIR/src" "$SCRIPT_DIR/include"; do
    [[ -d "$d" ]] && incs+=("-I$d")
  done
  printf '%s\n' "${incs[@]}"
}

# ---- Modes ----
if [[ "$MODE" == "clean" ]]; then
  rm -rf "$BUILD_DIR_DEFAULT" "$BUILD_DIR_ASAN" "$BUILD_DIR_CSAN"
  echo "Removed $BUILD_DIR_DEFAULT $BUILD_DIR_ASAN $BUILD_DIR_CSAN"
  exit 0
fi

if [[ "$MODE" == "openmp" ]]; then
  mkdir -p "$BUILD_DIR_DEFAULT"
  cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR_DEFAULT" -DCMAKE_BUILD_TYPE="$BUILD_TYPE" -DENABLE_CUDA=OFF
  cmake --build "$BUILD_DIR_DEFAULT" --target openmp_dlct -j"$(nproc_portable)"
  dlct_binary="$BUILD_DIR_DEFAULT/openmp_dlct"
  echo "OpenMP DLCT binary: $dlct_binary"
  if [[ ${#PROGRAM_ARGS[@]} -gt 0 ]]; then
    echo "Running OpenMP DLCT binary with arguments: ${PROGRAM_ARGS[*]}"
    "$dlct_binary" "${PROGRAM_ARGS[@]}"
  else
    echo "Run with: $dlct_binary [options]"
  fi
  exit 0
fi

if [[ "$MODE" == "asan" ]]; then
  mkdir -p "$BUILD_DIR_ASAN"
  cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR_ASAN" -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer -g" \
    -DENABLE_CUDA=OFF
  cmake --build "$BUILD_DIR_ASAN" --target openmp_dlct -j"$(nproc_portable)"
  echo "Running OpenMP DLCT under AddressSanitizer with --self-test..."
  "$BUILD_DIR_ASAN/openmp_dlct" --self-test
  echo "ASan check complete. For CUDA device checks, use: ./build.sh csan -- --self-test"
  exit 0
fi

if [[ "$MODE" == "csan" ]]; then
  command -v compute-sanitizer >/dev/null 2>&1 || {
    echo "ERROR: compute-sanitizer not in PATH. Install NVIDIA Compute Sanitizer." >&2
    exit 1
  }
  # Ensure no ASan env interferes with Compute Sanitizer
  unset LD_PRELOAD || true
  unset ASAN_OPTIONS || true

  mkdir -p "$BUILD_DIR_CSAN"
  cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR_CSAN" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
  cmake --build "$BUILD_DIR_CSAN" -j"$(nproc_portable)"

  [[ -x "$BUILD_DIR_CSAN/$TARGET_NAME" ]] || {
    echo "ERROR: CUDA binary not found: $BUILD_DIR_CSAN/$TARGET_NAME" >&2
    exit 1
  }

  echo "Running $TARGET_NAME under NVIDIA Compute Sanitizer (memcheck + leak detection)..."
  compute-sanitizer --tool memcheck --leak-check full \
    "$BUILD_DIR_CSAN/$TARGET_NAME" "${PROGRAM_ARGS[@]}"
  exit 0
fi

# Default: regular CMake build (Release/Debug)
mkdir -p "$BUILD_DIR_DEFAULT"
cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR_DEFAULT" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
cmake --build "$BUILD_DIR_DEFAULT" -j"$(nproc_portable)"

echo
# Check what was actually built and report accordingly
if [[ -x "$BUILD_DIR_DEFAULT/cuda_dlct" ]]; then
  echo "CUDA DLCT binary: $BUILD_DIR_DEFAULT/cuda_dlct"
elif [[ -x "$BUILD_DIR_DEFAULT/openmp_dlct" ]]; then
  echo "OpenMP DLCT binary: $BUILD_DIR_DEFAULT/openmp_dlct"
fi

echo "Tips:"
echo "  • OpenMP (CPU-only):  ./build.sh openmp"
echo "  • Device sanitizers:  ./build.sh csan -- --self-test"
echo "  • Host sanitizers:    ./build.sh asan"
