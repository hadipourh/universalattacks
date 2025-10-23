# Information Gain Estimation for Orthros Key Recovery

This tool estimates the information gain from multiple differential-linear distinguishers for key recovery attacks on Orthros. It implements Monte Carlo simulation to measure how much key entropy is reduced when combining distinguishers that share overlapping key bits.

## Overview

The tool performs two main functions:

1. **Verification experiments**: Reproduce the 6-round key recovery attack verification from the paper
2. **Histogram generation**: Estimate combined information gain from multiple 5-round or 7-round distinguishers over 10,000 randomly sampled keys

## Dependencies

### Required
- C++ compiler with C++17 support (GCC or Clang)
- OpenMP (for parallel execution)
- Boost C++ Libraries (math/distributions)

### macOS Setup
```bash
brew install boost libomp
```

### Linux Setup
```bash
sudo apt-get install libboost-all-dev libomp-dev
```

## Building

**macOS:**
```bash
make
```

The build uses `makefile.opt` for platform-specific settings. For macOS, OpenMP is configured via Homebrew's libomp.

**Linux:**
```bash
# Edit makefile.opt or use default settings
make
```

**Clean build:**
```bash
make clean
```

## Usage

The tool provides three commands:

### 1. Verify 6-round attack
```bash
./main verify-6r
```
Reproduces the 6-round key recovery verification experiment from the paper.

### 2. Generate 5-round histogram
```bash
./main histogram-5r [--histogram <file>]
```
Runs Monte Carlo simulation over 10,000 keys to estimate information gain from four 5-round distinguishers. By default, samples are written to `histogram-5.csv`.

Options:
- `--histogram <file>` or `--histogram=<file>`: Write samples to the specified file
- `--no-histogram`: Skip writing samples to disk (only print statistics)

### 3. Generate 7-round histogram
```bash
./main histogram-7r [--histogram <file>]
```
Same as histogram-5r but for four 7-round distinguishers. Default output: `histogram-7.csv`.

### Help
```bash
./main help
```
Display usage information.

## Output

For verification experiments (verify-6r), the tool prints:
- Key index mappings for each distinguisher
- Weak-key class distributions
- Verification results

For histogram experiments (histogram-5r, histogram-7r), the tool outputs:
- Overlapping key-bit positions across distinguishers
- Statistics for each sampled key (progress shown every 100 samples)
- Summary statistics: mean, median, min, max information gain
- Percentage of keys classified as strong/weak
- Optional CSV file with all 10,000 information gain values (one per line)

## Makefile Targets

```
make          - Build the main executable
make dep      - Generate dependency file
make clean    - Remove binaries and object files
make tar      - Create tarball of source files
```

## Platform Notes

### macOS
The `makefile.opt` contains macOS-specific settings. If you encounter libomp version issues, update the library path:

```makefile
LIB=-L/opt/homebrew/opt/libomp/lib -lomp
```

### Linux
Standard OpenMP flags should work:

```makefile
OPT=-m64 -O2 -funroll-loops -fopenmp -std=c++17
```

## How It Works

The information gain estimation follows this algorithm from the paper:

1. **Identify overlapping bits**: Find key-bit positions that appear in multiple distinguishers
2. **For each of 10,000 randomly sampled keys**:
   - Initialize probability sum to 0
   - For each of the 2^k possible values of the overlapping bits:
     - Initialize conditional probability to 1
     - For each distinguisher:
       - Identify which weak-key class the sampled key belongs to
       - Count how many 16-bit assignments in that class match the fixed overlapping bits
       - Divide by the total number of 16-bit assignments (across all classes) that match the fixed overlapping bits to get a conditional probability
       - Multiply into the product of conditional probabilities
     - Add the product to the probability sum
   - Average the sum over all 2^k overlapping-bit values to obtain probability p
   - Record -log₂(p) as the information gain for this key

Once the overlapping bits are fixed, the distinguishers become conditionally independent because each then reads a disjoint set of remaining key bits. This allows us to multiply the conditional probabilities.

## Files

- `main.cpp`: Main implementation
- `main.h`: Header file with DISTINGUISHER and ORTHROS classes
- `orthros.cpp`: Orthros cipher implementation
- `Makefile`: Build configuration
- `makefile.opt`: Platform-specific compiler settings
- `plot_hist.py`: Python script to visualize histogram data

## Files

- `main.cpp`: Main implementation
- `main.h`: Header file with DISTINGUISHER and ORTHROS classes
- `orthros.cpp`: Orthros cipher implementation
- `Makefile`: Build configuration
- `makefile.opt`: Platform-specific compiler settings
- `plot_hist.py`: Python script to visualize histogram data

## License

This project is licensed under the GNU General Public License v3.0 - see the source files for details.

## Contact

For questions, issues, or contributions, please contact:
- **Hosein Hadipour**: [hsn.hadipour@gmail.com](mailto:hsn.hadipour@gmail.com)
- **GitHub Issues**: [https://github.com/hadipourh/universalattacks](https://github.com/hadipourh/universalattacks)

## Authors

Copyright (C) 2025 Yosuke Todo and Hosein Hadipour
