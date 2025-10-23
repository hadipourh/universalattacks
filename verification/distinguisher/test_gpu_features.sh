#!/bin/bash
# GPU Testing Script for Orthros DLCT Tool
# Run this script on a GPU server with CUDA support
# Copyright (C) 2025 Hosein Hadipour

set -euo pipefail

echo "Orthros GPU Cryptanalysis Testing Suite"
echo "========================================"
echo

# First, build the CUDA targets (falls back to OpenMP binaries if CUDA is unavailable)
echo "1. Building CUDA targets..."
./build.sh release
echo

dlct_bin=""
dlct_label="CUDA"

if [[ -x "build/cuda_dlct" ]]; then
    dlct_bin="build/cuda_dlct"
elif [[ -x "build/openmp_dlct" ]]; then
    dlct_bin="build/openmp_dlct"
    dlct_label="OpenMP"
fi

if [[ -z "$dlct_bin" ]]; then
    echo "ERROR: Required binary not found under build/." >&2
    echo "Contents of build/:" >&2
    ls -la build/ >&2
    exit 1
fi

echo "✓ DLCT binary: $dlct_bin ($dlct_label)"
echo

if [[ "$dlct_label" != "CUDA" ]]; then
    echo "⚠ CUDA binary was not produced; running tests with $dlct_label target instead."
    echo "  Some performance metrics may not reflect GPU execution."
    echo
fi

echo "2. Running self-test..."
"$dlct_bin" --self-test
echo

echo "3. Testing 5-round PRF Distinguisher 0 (from paper, RB=1, RU=1, RM=4, RL=0, offset=2)"
"$dlct_bin" -r 5 -o 3 -m 2 \
    --diff-left 00002020000000000000000000000000 \
    --diff-right 00000000000044000000000000000000 \
    --mask 00000000000000000000000000002022 \
    -e 10 -s 10 -b 1024 -t 1024 -v 
echo

echo "✓ Test completed successfully!"