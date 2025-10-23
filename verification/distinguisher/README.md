# Distinguisher Verification Toolkit

A modular, GPU-accelerated toolkit for large-scale differential-linear cryptanalysis of Orthros PRF. It provides CUDA implementations for high-performance correlation computation, OpenMP fallbacks for CPU-only hosts, and visualization utilities for post-processing experimental data.

- [Distinguisher Verification Toolkit](#distinguisher-verification-toolkit)
  - [Overview](#overview)
  - [Core Features](#core-features)
  - [Architecture](#architecture)
  - [Prerequisites](#prerequisites)
    - [CUDA Path](#cuda-path)
    - [CPU/OpenMP Path](#cpuopenmp-path)
    - [Visualization](#visualization)
  - [Quick Start](#quick-start)
  - [Building From Source](#building-from-source)
  - [Usage Guide](#usage-guide)
    - [Self-Tests](#self-tests)
    - [Differential-Linear Analysis](#differential-linear-analysis)
    - [Visualization Workflow](#visualization-workflow)
  - [Example Scenarios](#example-scenarios)
  - [Testing and Quality Assurance](#testing-and-quality-assurance)
  - [Contributing](#contributing)
  - [License](#license)
  - [Contact](#contact)
  - [Acknowledgments](#acknowledgments)

## Overview

The toolkit provides GPU-accelerated and CPU-based tooling for differential-linear cryptanalysis.

Key components include:

- CUDA Differential-Linear Cryptanalysis Tool (DLCT) for Orthros
- OpenMP counterpart for CPU-only environments
- Automated GPU validation scripts
- Python visualization workflows for correlation distributions

## Core Features

- **GPU acceleration** for rapid experimentation on NVIDIA hardware (compute capability ≥ 8.0 recommended)
- **Support for Orthros branches** (left/right branch modes and PRF mode)
- **Configurable round offset** for reduced-round analysis
- **Flexible experiment configuration** with tunable rounds, sample powers, grid dimensions, and experiment counts
- **Robust entropy sourcing** via hardware RDRAND with secure fallbacks
- **Visualization utilities** for statistical summaries and publication-ready plots
- **Self-test coverage** for both CUDA kernels and host utilities

## Architecture

```
.
├── src/
│   ├── cuda_dlct.cu             # CUDA DLCT implementation for Orthros
│   ├── openmp_dlct.cpp          # OpenMP DLCT implementation
│   ├── correlation_analysis.cpp # Shared statistics and reporting
│   ├── visualizer.py            # Python plotting utilities
│   ├── orthros_core.h           # Orthros cipher primitives and round functions
│   ├── common_types.h           # Shared configuration/data structures
│   └── openmp/                  # OpenMP utilities
│       ├── dlct_runner.cpp      # OpenMP experiment runner
│       ├── random_utils.cpp     # Random number generation
│       └── random_utils.h       # Random utilities header
├── build.sh                     # Unified build helper
├── visualize.sh                 # Visualization wrapper
├── test_gpu_features.sh         # Automated GPU regression suite
├── requirements.txt             # Python dependencies
└── CMakeLists.txt               # CMake build definition
```

## Prerequisites

### CUDA Path

- NVIDIA GPU with compute capability 8.0 or newer
- CUDA Toolkit 11.0+
- CMake 3.18+
- GCC or Clang with C++17 support
- At least 8 GB GPU RAM recommended for high sample powers (≥ 2^27)

### CPU/OpenMP Path

- Modern C/C++ compiler with C++17 support
- OpenMP runtime (install `libomp` on macOS)

### Visualization

- Python 3.9+
- Install dependencies with `pip install -r requirements.txt`

## Quick Start

```bash
# Clone and enter the repository
git clone <repository-url>
cd code/orthros/keyrecoveryandverification/gpu

# Build with auto-detected backends (CUDA when available)
./build.sh

# Exercise GPU capabilities
./test_gpu_features.sh
```

## Building From Source

The `build.sh` helper centralizes build modes:

| Command                            | Description                                                                                  |
| ---------------------------------- | -------------------------------------------------------------------------------------------- |
| `./build.sh`                     | Release build. Produces CUDA binaries when `nvcc` is available, otherwise OpenMP binaries. |
| `./build.sh debug`               | Debug build with symbols and additional checks.                                              |
| `./build.sh openmp`              | Force OpenMP-only build (no CUDA).                                                           |
| `./build.sh clean`               | Remove build directories (`build/`, `build-csan/`).                     |
| `./build.sh csan -- --self-test` | Run the CUDA DLCT binary under NVIDIA Compute Sanitizer.                                     |

Binaries reside in the corresponding build directory (`build/` or `build-csan/`).

## Usage Guide

### Self-Tests

```bash
# GPU or OpenMP binaries (depending on build results)
./build/cuda_dlct --self-test   # or ./build/openmp_dlct
```

### Differential-Linear Analysis

```bash
./build/cuda_dlct \
  -r 5 \
  -o 3 \
  -m 2 \
  -d 00002020000000000000000000000000 \
  --diff-right 00000000000044000000000000000000 \
  --mask 00000000000000000000000000002022 \
  -e 2 \
  -s 10 \
  -b 256 \
  -t 512 \
  -v
```

Key arguments:

- `-r, --rounds` – number of cipher rounds (1 to 16, default: 4)
- `-o, --offset` – starting round offset (0 to 15, default: 4)
- `-m, --mode` – branch mode: 0=left, 1=right, 2=PRF (default: 0)
- `--diff-left` – left branch input difference (32 hex nibbles = 128 bits)
- `--diff-right` – right branch input difference (32 hex nibbles = 128 bits)
- `--mask` – output mask (32 hex nibbles = 128 bits)
- `--wildcard-mask` – mask with 'x' for variable positions
- `-e, --experiments` – number of experiment iterations (default: 16)
- `-s, --sample-power` – samples per thread = 2^s (default: 18, max: 30)
- `-b, --blocks` – number of CUDA blocks (default: 256)
- `-t, --threads` – threads per block (default: 512)
- `-g, --gendlct` – generate full DLCT table with single-bit patterns
- `-v, --verbose` – enable detailed output

### Visualization Workflow

```bash
pip install -r requirements.txt

# Run analysis and generate correlation data
./build/cuda_dlct -r 8 -o 4 -m 0 -e 32 -s 20 -v

# Visualize results
./visualize.sh
# or
python3 src/visualizer.py orthros_*.csv
```

Generated outputs include PDF, PNG, CSV, and textual summaries highlighting distribution fits and descriptive statistics.

## Example Scenarios

| Objective                                | Command Sequence                                                                                                                                                                           |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 5-round PRF Distinguisher 0 (from paper) | `./build/cuda_dlct -r 5 -o 3 -m 2 -d 00002020000000000000000000000000 --diff-right 00000000000044000000000000000000 --mask 00000000000000000000000000002022 -e 2 -s 10 -b 256 -t 512 -v` |
| Wildcard mask example                    | `./build/cuda_dlct -r 5 -o 3 -m 0 -d 00002020000000000000000000000000 --diff-right 00000000000000000000000000000000 --wildcard-mask 0000000000000000000000000000x0xx -e 10 -s 10 -v`     |
| Run self-test                            | `./build/cuda_dlct --self-test`                                                                                                                                                          |

## Testing and Quality Assurance

1. Execute CUDA/OpenMP self-tests: `./build/cuda_dlct --self-test`
2. Launch the GPU validation test (runs 5-round PRF Distinguisher 0): `./test_gpu_features.sh`
3. Optional CUDA sanitizer: `./build.sh csan -- --self-test`

Include hardware details, compiler versions, and full command lines when reporting issues.

## Contributing

1. Fork the repository and create a feature branch (`git checkout -b feature/my-change`).
2. Develop changes with clear commits and accompanying tests when applicable.
3. Run self-tests and sanitizer builds before opening a pull request.
4. Follow existing C++/CUDA style conventions and document complex logic.

When filing issues, provide:

- Build mode (release, openmp, csan, etc.)
- Compiler and/or CUDA toolkit versions
- Exact commands that reproduce the problem
- Logs, stack traces, or sanitizer reports

## License

Distributed under the GNU General Public License v3.0. Refer to the `LICENSE` file for the full text.

## Contact

For questions, issues, or contributions, please contact:

- **Hosein Hadipour**: [hsn.hadipour@gmail.com](mailto:hsn.hadipour@gmail.com)
- **GitHub Issues**: [https://github.com/hadipourh/universalattacks](https://github.com/hadipourh/universalattacks)

## Acknowledgments

- Orthros PRF specification and design
- NVIDIA CUDA platform and tooling
- OpenMP community and runtime maintainers
- Researchers advancing differential-linear cryptanalysis
