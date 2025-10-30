# Attack Discovery Tools

Tools for discovering and analyzing differential-linear distinguishers for the Orthros PRF using constraint programming and empirical verification.

- [Attack Discovery Tools](#attack-discovery-tools)
  - [Features](#features)
  - [Dependencies](#dependencies)
    - [Required Software](#required-software)
    - [Python Packages](#python-packages)
    - [LaTeX Requirements](#latex-requirements)
  - [Installation](#installation)
    - [1. Set Up Python Environment](#1-set-up-python-environment)
    - [2. Install MiniZinc](#2-install-minizinc)
    - [3. Install Gurobi (Optional)](#3-install-gurobi-optional)
  - [Building](#building)
    - [Compile the C Verification Tool](#compile-the-c-verification-tool)
    - [Clean Build Artifacts](#clean-build-artifacts)
    - [Verify Build](#verify-build)
  - [Usage](#usage)
    - [Single-Branch Distinguisher Search](#single-branch-distinguisher-search)
    - [Full PRF Distinguisher Search](#full-prf-distinguisher-search)
    - [Common Options](#common-options)
    - [Direct C Tool Usage](#direct-c-tool-usage)
  - [Output Files](#output-files)
    - [Generated Files](#generated-files)
    - [Compiling LaTeX Visualizations](#compiling-latex-visualizations)
  - [Troubleshooting](#troubleshooting)
    - [Common Issues](#common-issues)
  - [License](#license)
  - [Contact](#contact)


## Features

- **Automated Distinguisher Discovery**: Find optimal DL distinguishers using MiniZinc constraint solver
- **Single-Branch Analysis**: Test individual branches of the Orthros PRF independently
- **Full PRF Analysis**: Analyze the complete two-branch PRF construction
- **Correlation Verification**: Empirically verify distinguisher correlations using C implementation
- **Differential Effect Computation**: Calculate differential probabilities using Gurobi optimizer
- **Visualization**: Generate LaTeX-based graphical representations of discovered distinguishers

## Dependencies

### Required Software
- **Python 3.8+**: Core scripting environment
- **GCC compiler**: C99 support required for building verification tools
- **MiniZinc**: Constraint programming solver for distinguisher search
- **Gurobi Optimizer** (optional): For computing differential effects via MILP

### Python Packages
- `minizinc`: MiniZinc Python interface
- `pyyaml`: Configuration file parsing
- `gurobipy`: Gurobi Python API (only if using differential effect computation)

### LaTeX Requirements
The visualization tools generate LaTeX output files. The following are required for compiling these files:

**Required LaTeX packages:**
- `tikz`: For drawing distinguisher diagrams
- `xcolor`: For colored graphics
- Standard LaTeX distributions (TeX Live, MiKTeX, MacTeX)

**Required style files in `tikzstyle/` directory:**
- `tikzlibrarycipher.code.tex`: Custom TikZ library for cipher visualization
- `tugcolors.sty`: Color scheme definitions
- `orthros.sty`: Orthros-specific drawing styles
- `spn.sty`: SPN network visualization utilities

**Note:** These style files are included in the `tikzstyle/` directory and will be automatically found by LaTeX when using the provided `.latexmkrc` configuration.

## Installation

### 1. Set Up Python Environment

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install --upgrade pip
```

Install required Python packages:

```bash
pip install minizinc pyyaml
```

If you plan to use differential effect computation:

```bash
pip install gurobipy
```

### 2. Install MiniZinc

MiniZinc is required for constraint-based distinguisher search.

- **macOS:**: `brew install minizinc`
- **Linux:**: Follow the installation guide at: https://github.com/hadipourh/minizinc-installer-linux
- **Windows:**: Download and install from: https://www.minizinc.org/

### 3. Install Gurobi (Optional)

Gurobi is only required if you want to compute differential effects using `-de 1` option.

1. Download from: https://www.gurobi.com/
2. Obtain a free academic license or trial license
3. Install the software and activate your license

**For Linux users**, refer to this helper script: https://github.com/hadipourh/grabgurobi

## Building

### Compile the C Verification Tool

**Important:** The `difflin` C program must be compiled before running Python scripts.

Build the tool:

```bash
make
```

This compiles:
- `difflin`: Empirical correlation verification tool
- `orthros.o`: Orthros cipher implementation
- `utils.o`: Utility functions

### Clean Build Artifacts

To remove compiled binaries and object files:

```bash
make clean
```

### Verify Build

Test the compiled tool with the default 4-round example:

```bash
./difflin 6 4 70000000000000800000000000c00000 \
              70000000000000800000000000c00000 \
              000000000000ccc0000000000ccc0ccc 20 8 0
```

Or simply run without arguments to use built-in defaults:

```bash
./difflin
```

Expected output will show the correlation value, its base-2 logarithm, and sign.

## Usage

### Single-Branch Distinguisher Search

Search for differential-linear distinguishers on a single branch of Orthros.

**Syntax:**
```bash
python3 attackbranch.py -RU <rounds_upper> -RM <rounds_middle> -RL <rounds_lower> -offset <offset> -branchtype <0|1> [options]
```

**Example - 5-round distinguisher on branch 0:**
```bash
python3 attackbranch.py -RU 0 -RM 5 -RL 0 -offset 4 -branchtype 0
```

**Required Parameters:**
- `-RU <int>`: Number of rounds for upper part (differential phase, E_U)
- `-RM <int>`: Number of rounds in the middle part (E_M)
- `-RL <int>`: Number of rounds for lower part (linear phase, E_L)
- `-offset <int>`: Starting round number (0-indexed)
- `-branchtype <0|1>`: Branch selection (0 = left branch, 1 = right branch)

You can see the full list of options by running:
```bash
python3 attackbranch.py -h
```

### Full PRF Distinguisher Search

Search for distinguishers on the complete two-branch Orthros PRF construction.

**Syntax:**
```bash
python3 attackprf.py -RU <rounds_upper> -RM <rounds_middle> -RL <rounds_lower> -offset <offset> -RB <rb_offset> -d <degree> [options]
```

**Example - 6-round PRF distinguisher (Distinguisher 0 from paper):**
```bash
python3 attackprf.py -RU 1 -RM 4 -RL 0 -offset 2 -RB 1 -d 10
```

**Required Parameters:**
- `-RU <int>`: Number of rounds for upper part (differential phase, E_U)
- `-RM <int>`: Number of rounds in the middle part (E_M)
- `-RL <int>`: Number of rounds for lower part (linear phase, E_L)
- `-offset <int>`: Starting round number for the entire distinguisher (0-indexed)
- `-RB <int>`: Number of rounds included in the key recovery phase (default: 1)
- `-d <int>`, `--degree <int>`: Sample size exponent for middle part correlation estimation (N = 2^degree)

You can see the full list of options by running:
```bash
python3 attackprf.py -h
```

**Note:** The `--degree` parameter controls the number of samples used to empirically estimate correlations of the middle part (E_M). Typical values range from 8 to 16.

### Common Options

Both `attackbranch.py` and `attackprf.py` support the following options:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-np <threads>` | int | 8 | Number of parallel threads for constraint solver |
| `-tl <seconds>` | int | 36000 | Time limit for solver (10 hours) |
| `-sl <solver>` | str | cp-sat | Constraint programming solver to use |
| `-o <file>` | str | auto | Output file name for results |
| `-de <0\|1>` | int | 0 | Compute differential effect using Gurobi MILP (requires Gurobi license) |

### Direct C Tool Usage

The `difflin` tool provides low-level access for empirical correlation verification. This is useful for testing specific distinguisher configurations or integrating with external tools.

**Syntax:**
```bash
./difflin <rounds> <offset> <inputdiff0> <inputdiff1> <outputmask> <deg> <N> <mode>
```

**Parameters:**

- **`<rounds>`**: Number of rounds to evaluate
  - **Type:** Integer
  - **Description:** Specifies the number of Orthros rounds to analyze
  - **Example:** `4` evaluates a 4-round distinguisher

- **`<offset>`**: Round offset  
  - **Type:** Integer
  - **Description:** Starting round number (0-indexed). Enables testing distinguishers that begin at a specific round rather than round 0
  - **Example:** `4` starts analysis from round 4
  - **Use case:** Testing middle-round segments of the cipher

- **`<inputdiff0>`**: Input difference for branch 0
  - **Type:** 32-character hexadecimal string (128 bits)
  - **Description:** Input difference applied to the left branch of Orthros PRF
  - **Format:** Each character is a 4-bit nibble (0-F, case-insensitive)
  - **Example:** `00000000001000000000000000000000`

- **`<inputdiff1>`**: Input difference for branch 1
  - **Type:** 32-character hexadecimal string (128 bits)
  - **Description:** Input difference applied to the right branch of Orthros PRF
  - **Format:** Each character is a 4-bit nibble (0-F, case-insensitive)
  - **Example:** `00000000001000000000000000000000`

- **`<outputmask>`**: Output mask
  - **Type:** 32-character hexadecimal string (128 bits)
  - **Description:** Linear mask applied to the output for correlation computation
  - **Format:** Each character is a 4-bit nibble (0-F, case-insensitive)
  - **Example:** `00000000000000000000101100000000`
  - **Note:** Defines the linear mask at the output of the analyzed rounds

- **`<deg>`**: Degree parameter (sample size exponent)
  - **Type:** Integer (typically 18-24)
  - **Description:** Determines number of samples for correlation estimation: N₁ = 2^deg
  - **Example:** `20` → $2^{20}$ = 1048576 plaintext pairs per experiment
- **`<N>`**: Number of experiments
  - **Type:** Integer (typically 4-16)
  - **Description:** Number of independent experiments to average for statistical confidence
  - **Example:** `8` runs 8 independent trials with different keys and averages results
  - **Note:** Each experiment uses fresh random key and plaintext samples
  - **Recommendation:** Use 8+ for reliable correlation estimates

- **`<mode>`**: Branch mode
  - **Type:** Integer (0, 1, or 2)
  - **Values:**
    - `0` = Left branch only (single-branch distinguisher)
    - `1` = Right branch only (single-branch distinguisher)
    - `2` = Full PRF mode (both branches XORed)
  - **Example:** `0` tests the left branch in isolation

**Example Usage:**

Test with the built-in default parameters:
```bash
./difflin 4 4 00000000000000000000110111100000 \
             00000000000000008008000000000000 \
             00000000ff0f00000000000000000000 20 8 1
```

Or run without arguments to use defaults:
```bash
./difflin
```

**Output Format:**

```
Correlation: 1.000000
Logarithm2 of Correlation: 0.000000
Sign: 1
```

**Interpretation:**
- **Correlation:** Empirical correlation value (range: -1 to 1)
  - Values close to 0 indicate no bias (uniform distribution)
  - Values away from 0 indicate distinguishable bias
- **Logarithm2 of Correlation:** Base-2 log of absolute correlation (useful for comparison)
  - Smaller magnitude = stronger distinguisher
  - Example: -5.0 means correlation ≈ 2⁻⁵ = 0.03125
- **Sign:** Direction of correlation bias (+1 or -1)

## Output Files

The Python tools generate several output files:

### Generated Files

| File Type | Description | Format |
|-----------|-------------|--------|
| `output*.tex` | Distinguisher visualization | LaTeX source |
| `output*.pdf` | Compiled visualization | PDF (after running latexmk) |
| `output*.csv` | Correlation data and statistics | CSV |
| `minizinc-python.log` | MiniZinc solver logs | Text |
| `debug_output.txt` | Detailed CP solver output | Text |

### Compiling LaTeX Visualizations

To generate PDF from LaTeX output:

```bash
latexmk -pdf output_distinguisher.tex
```

Or manually:

```bash
pdflatex output_distinguisher.tex
```

**Requirements:** The LaTeX files require the following packages:
- `tikz`
- `xcolor`
- Standard LaTeX distributions (TeX Live, MiKTeX, MacTeX)

## Troubleshooting

### Common Issues

**Problem:** `difflin: command not found`
- **Solution:** Run `make` to compile the C tool first

**Problem:** MiniZinc solver not found
- **Solution:** Ensure MiniZinc is installed and in your system PATH
- **Check:** Run `minizinc --version`

**Problem:** Gurobi license error
- **Solution:** Verify Gurobi license is activated, or use `-de 0` to disable differential effect computation

**Problem:** Python import errors
- **Solution:** Activate virtual environment and reinstall packages:
  ```bash
  source .venv/bin/activate
  pip install minizinc pyyaml gurobipy
  ```

**Problem:** Solver times out
- **Solution:** Increase time limit with `-tl <seconds>` or try a simpler configuration (fewer rounds)

## License

This software is released under the GNU General Public License v3.0. See LICENSE file for details.

## Contact

For questions, issues, or contributions:
- Open an issue on GitHub
- Contact: Hosein Hadipour <hsn.hadipour@gmail.com>