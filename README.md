# Universal Attacks: Turning Multiple Key-Dependent Attacks into Universal Attacks

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This repository contains the companion code for the paper **"Turning Multiple Key-Dependent Attacks into Universal Attacks"** by Hosein Hadipour, Yosuke Todo, Mostafizar Rahman, Maria Eichlseder, Ravi Anand, and Takanori Isobe.

## Overview

Key-dependent attacks are effective only for specific weak-key classes, limiting their practical impact. This work presents a statistical framework that combines multiple key-dependent distinguishers into universal attacks covering the full key space.

Applied to **Orthros-PRF**, a sum-of-permutations design where differential-based distinguishers hold only for a fraction of keys, this yields the **first universal 8-round differential-linear key-recovery attack** with median time complexity 2^119.58, whereas prior work reached at most 7 rounds in the weak-key setting.

## Key Contributions

### 1. Generic Framework for Universal Attacks

We propose a generic framework for combining multiple key-dependent attacks into a universal attack. The framework tests the secret key against multiple weak-key distinguishers using log-likelihood ratio statistics and aggregates the results to reduce effective key entropy across the full key space.

### 2. First Universal 8-Round Attack on Orthros-PRF

We apply our framework to Orthros-PRF and obtain the first universal 8-round differential-linear key-recovery attack with median time complexity 2^119.58, whereas prior work reached at most 7 rounds in the weak-key setting. The attack complexity varies across the key space: median complexity is 2^119.58 (for at least 50% of keys), while worst-case complexity is 2^126.92 (for approximately 48% of keys that are strong against all employed distinguishers).

### 3. Automated Attack Discovery with MILP

We develop a CP-MILP model for automated discovery of complete differential-linear attacks that integrates distinguisher search with key recovery. We extend the open-source [S-box Analyzer](https://github.com/hadipourh/sboxanalyzer) to support deterministic propagation of differences and linear masks, and discover multidimensional distinguishers covering up to 10 rounds in each Orthros branch, **improving the prior best by 4 rounds**.

### 4. Open-Source Implementation and Experimental Validation

We provide open-source implementations of our attack discovery tools and statistical framework, allowing independent verification and supporting future work on other primitives. Experimental observations that motivate and validate our statistical modeling are included in the paper. 

## Repository Organization

```
code/public/
├── README.md                          # This file
├── attack/                            # Attack discovery tools
│   ├── README.md                      # Detailed documentation
│   ├── attackprf.py                   # CP/MILP-based PRF attack search
│   ├── attackprf.mzn                  # MiniZinc model for PRF attacks
│   ├── attackbranch.py                # CP/MILP-based branch attack search
│   ├── attackbranch.mzn               # MiniZinc model for branch attacks
│   ├── difflin.c                      # Empirical correlation verification tool
│   ├── difflinprfmilp.py              # MILP model for PRF distinguishers
│   ├── draw.py                        # Visualization for branch attacks
│   ├── drawprf.py                     # Visualization for PRF attacks
│   ├── orthros.c                      # Orthros implementation
│   ├── input.yaml                     # Configuration file
│   └── tikzstyle/                     # LaTeX/TikZ style files
└── verification/                      # Statistical verification tools
    ├── distinguisher/                 # GPU-accelerated DLCT computation
    │   ├── README.md                  # Detailed documentation
    │   ├── CMakeLists.txt             # Build configuration
    │   ├── build.sh                   # Build script
    │   ├── test_gpu_features.sh       # Automated test suite
    │   └── src/                       # Source files
    │       ├── cuda_dlct.cu           # CUDA implementation
    │       ├── openmp_dlct.cpp        # OpenMP implementation
    │       ├── correlation_analysis.cpp  # Correlation analysis
    │       └── visualizer.py          # Result visualization
    └── keyrecovery/                   # Universal attack implementation
        ├── README.md                  # Detailed documentation
        ├── Makefile                   # Build configuration
        ├── keyrecovery.c              # Key recovery implementation (CPU)
        ├── keyrecovery.cu             # Key recovery implementation (GPU)
        ├── keyanalysis.cpp            # Key analysis tool
        ├── visualize_keyanalysis.py   # Visualization scripts
        └── llr/                       # LLR analysis and information gain
            ├── README.md              # Detailed documentation
            ├── Makefile               # Build configuration
            ├── main.cpp               # Information gain estimation
            ├── main.h                 # Header file
            ├── orthros.cpp            # Orthros cipher implementation
            └── plot_hist.py           # Visualization tools
```

## Repository Structure

This repository is organized into two main components:

### 1. Attack Discovery ([`attack/`](attack/))

Tools for automated discovery of differential-linear distinguishers and key-recovery attacks using CP/MILP optimization and constraint programming.

**Key Features:**
- MILP-based search for multidimensional differential-linear distinguishers
- Automated integration of distinguisher search with key recovery
- MiniZinc models for constraint-based attack discovery
- Visualization tools for attack paths

**Main Tools:**
- [`attackprf.py`](attack/attackprf.py) - CP/MILP-based search for PRF distinguishers
- [`attackbranch.py`](attack/attackbranch.py) - CP/MILP-based search for branch distinguishers
- [`difflin`](attack/difflin.c) - C tool for empirical correlation verification
- [`draw.py`](attack/draw.py), [`drawprf.py`](attack/drawprf.py) - LaTeX/TikZ visualization generators

**[→ See detailed documentation](attack/README.md)**

### 2. Statistical Verification ([`verification/`](verification/))

GPU-accelerated tools for statistical verification and key recovery using the universal attack framework.

#### 2.1 Distinguisher Verification ([`verification/distinguisher/`](verification/distinguisher/))

High-performance CUDA and OpenMP implementations for computing differential-linear correlation and verifying weak-key classifications.

**Key Features:**
- CUDA implementation optimized for NVIDIA GPUs
- OpenMP implementation for multi-core CPUs
- Support for single and multi-branch distinguishers
- Parallel correlation estimation with configurable precision

**[→ See detailed documentation](verification/distinguisher/README.md)**

#### 2.2 Key Recovery ([`verification/keyrecovery/`](verification/keyrecovery/))

Implementation of the universal key-recovery framework combining multiple weak-key distinguishers.

**Main Components:**
- **Distinguisher-based classification**: Test keys against weak-key distinguishers using LLR statistics
- **Information gain estimation**: Monte Carlo simulation to measure combined entropy reduction
- **Full attack implementation**: 6-round and 8-round universal attacks on Orthros-PRF

**Subdirectories:**
- [`llr/`](verification/keyrecovery/llr/) - Log-likelihood ratio analysis and information gain estimation

**[→ See detailed documentation](verification/keyrecovery/README.md)**

## Quick Start

### Prerequisites

Different components have different requirements. See individual README files for details.

**For Attack Discovery:**
- Python 3.8+
- MiniZinc 2.8+ with Gurobi solver
- LaTeX with TikZ (for visualization)

**For Statistical Verification:**
- CUDA Toolkit 11.0+ (for GPU acceleration)
- OpenMP support (for CPU parallelization)
- C++17 compiler (GCC 7+ or Clang 10+)

### Example Workflow

1. **Discover distinguishers** using CP/MILP-based search:
   ```bash
   cd attack
   # Example: 6-round PRF distinguisher (Distinguisher 0 from paper)
   python3 attackprf.py -offset 2 -RB 1 -RU 1 -RM 4 -RL 0 -d 10
   ```

2. **Verify distinguisher correlations** on GPU:
   ```bash
   cd verification/distinguisher
   ./build.sh
   # Example: 5-round PRF Distinguisher 0 from paper
   ./build/cuda_dlct -r 5 -o 3 -m 2 \
     -d 00002020000000000000000000000000 \
     --diff-right 00000000000044000000000000000000 \
     --mask 00000000000000000000000000002022 \
     -e 2 -s 10 -b 256 -t 512 -v
   ```

3. **Estimate information gain** from multiple distinguishers:
   ```bash
   cd verification/keyrecovery/llr
   make
   # Run Monte Carlo simulation over 10,000 keys for 5-round distinguishers
   ./main histogram-5r
   ```

## Citation

If you use this code in your research, please cite our paper:```bibtex
@article{hadipour2025universal,
  title={Turning Multiple Key-Dependent Attacks into Universal Attacks},
  author={Hadipour, Hosein and Todo, Yosuke and Rahman, Mostafizar and Eichlseder, Maria and Anand, Ravi and Isobe, Takanori},
  year={2025},
  note={To appear}
}

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details or visit [https://www.gnu.org/licenses/gpl-3.0](https://www.gnu.org/licenses/gpl-3.0).

## Contact

For questions, issues, or contributions, please contact:
- **Hosein Hadipour**: [hsn.hadipour@gmail.com](mailto:hsn.hadipour@gmail.com)
- **GitHub Issues**: [https://github.com/hadipourh/universalattacks](https://github.com/hadipourh/universalattacks)

## Acknowledgments and Short Story of the Work

I began this work during the final year of my (Hosein Hadipour) Ph.D., when I visited Takanori Isobe's group for about two weeks at the University of Hyogo, Japan, in early 2024. During that visit, I discovered several distinguishers for Orthros-PRF that significantly outperformed existing ones.

Shortly afterward, at FSE 2024 in Belgium, I shared these findings with Yosuke Todo, who was independently working on sum-of-permutation designs including Orthros-PRF. He identified a critical limitation in my distinguishers: the distinguishers were effective only for a fraction of the keys. This observation sparked a collaboration to develop a method for transforming the weak-key distinguishers into universal attacks, which culminated in the framework presented in this work.

The publication was delayed so long as I was busy with other commitments and we refined our methods and verified our results through extensive experimentation. Most of the research was completed after my Ph.D., during the first few months of my postdoctoral position at Ruhr University Bochum, Germany, under the supervision of Gregor Leander.

I would like to express my gratitude to Takanori Isobe and his lab members, Mostafizar Rahman and Ravi Anand, for their warm hospitality during my stay in Japan and their invaluable support throughout this project. I am also deeply grateful to Yosuke Todo and Maria Eichlseder for their collaboration, insightful discussions, and contributions that shaped this work. I also thank Gregor Leander for his support during my postdoctoral research.

---

**Note**: This repository contains the code for the published research. For the complete paper including detailed proofs and additional results, please refer to the publication.
