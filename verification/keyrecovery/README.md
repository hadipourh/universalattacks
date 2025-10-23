# Orthros Key Recovery Attack Tool

Tools for key recovery attacks on Orthros PRF using differential-linear cryptanalysis.

## Dependencies

### Required
- GCC or Clang with C99 support
- OpenMP (for parallel execution)
- Python 3.8+
- CUDA Toolkit (optional, for GPU acceleration)

### macOS OpenMP Setup
```bash
brew install libomp
```

### Python Packages
```bash
pip install numpy matplotlib scipy
```

### Optional (for GPU acceleration)
- NVIDIA CUDA Toolkit
- CUDA-capable GPU

## Building

The Makefile contains multiple pre-configured attack scenarios for different Orthros variants. Each configuration block in the Makefile defines parameters for a specific distinguisher.

### Configuration Blocks

The Makefile includes several commented configuration blocks. These are (r+1)-round attacks based on r-round distinguishers:
- **5-round_version-0**: 5-round attack (based on 4-round distinguisher) with correlation 2^{-2.7}
- **5-round_version-1**: 5-round attack (based on 4-round distinguisher) with correlation 2^{-2.3}
- **5-round_version-2**: 5-round attack (based on 4-round distinguisher) with correlation 2^{-5.7}
- **6-round_version-1**: 6-round attack (based on 5-round distinguisher, currently active)

**To switch configurations:**
1. Comment out the currently active configuration block
2. Uncomment the desired configuration block
3. Rebuild

### Build Targets

```bash
# Build CPU version (default)
make

# Build GPU version (requires CUDA)
make keyrecovery_gpu
```

The binaries are named `keyrecovery` (CPU) and `keyrecovery_gpu` (GPU).

## Usage

### Quick Start

```bash
# Run with default active configuration
make run

# Run with custom degree (data complexity)
make run DEG=20

# Run on GPU
make run-gpu

# Quick test with fewer samples
make test
```

### Direct Execution

The tool can be run directly with all parameters specified:

```bash
./keyrecovery <num_of_rounds> <rndOffset> <outputmask> <inputDiffPattern> \
              <keyValType> <plainNibValType> <DEG> <N> <trailFolder> \
              <dy1_hi> <dy1_lo> <dy2_hi> <dy2_lo> [num_threads]
```

**Parameters:**
- `num_of_rounds`: Total number of Orthros rounds
- `rndOffset`: Starting round offset
- `outputmask`: Output mask as 32 hex nibbles
- `inputDiffPattern`: Input difference pattern as comma-separated 0/1 values (32 values)
- `keyValType`: Key type (0=Weak, 1=Random, 2=Strong)
- `plainNibValType`: Plaintext nibble value type
- `DEG`: Data complexity as power of 2 (samples = 2^DEG)
- `N`: Number of independent experiments
- `trailFolder`: Name for results folder (e.g., "R5_V1")
- `dy1_hi`, `dy1_lo`: Left branch differential nibbles (hex)
- `dy2_hi`, `dy2_lo`: Right branch differential nibbles (hex)
- `num_threads`: Optional, number of OpenMP threads (default: auto-detect)

**Example:**
```bash
./keyrecovery 6 2 "00000000000000000000000000008088" \
  "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0" \
  0 0 28 4 "R5_V1" 0x8 0x2 0x4 0x4
```

### Configuration Parameters in Makefile

Each configuration block defines these variables:

```makefile
dy1 = 0x8 0x2              # Left branch differential (hi lo nibbles)
dy2 = 0x4 0x4              # Right branch differential (hi lo nibbles)
NAME="R5_V1"               # Trail/result folder name
NUM_OF_ROUNDS = 6          # Total rounds
RND_OFFSET = 2             # Starting offset
OUTPUT_MASK = "..."        # 32 hex nibbles
INPUT_DIFF_PATTERN = "..." # 32 comma-separated 0/1 values
KEY_VAL_TYPE = 0           # 0=Weak, 1=Random, 2=Strong
PLAIN_NIB_VAL_TYPE = 0     # Plaintext type
DEG = 28                   # Data complexity (2^28 samples)
N = 4                      # Number of experiments
CORR = 11                  # -log2(correlation)
```

### GPU Execution

For GPU acceleration:

```bash
# Build GPU version
make keyrecovery_gpu

# Run with custom CUDA configuration
make run-gpu CUDA_BLOCKS=512 CUDA_THREADS=512 DEG=25

# Or run directly
CUDA_THREADS_PER_BLOCK=256 CUDA_BLOCKS=256 ./keyrecovery_gpu <parameters...>
```

## Visualization

### Visualize Results

```bash
# Visualize latest run for current configuration
make viz NAME=R5_V1

# Visualize all experiments
make viz-all

# Compare multiple experiments
make viz-compare
```

The visualization tools generate:
- Correlation distribution plots (PDF/PNG)
- Statistical analysis
- Key recovery success probability

### Python Analysis Tool

For detailed key space analysis:

```bash
# Build analysis tool
make keyanalysis

# Run analysis
./keyanalysis
```

Or use the Python wrapper:

```bash
python3 keyanalysis.py
python3 visualize_keyanalysis.py
```

## Output Files

Results are stored in `data/<NAME>/`:
- `out_*.txt`: Raw correlation data
- `*.pdf`, `*.png`: Visualization plots
- `keyanalysis_results.txt`: Key space analysis

## Performance Notes

### Data Complexity
The `DEG` parameter controls data complexity as 2^DEG samples:
- DEG=18: 262K samples (~quick test)
- DEG=20: 1M samples (~fast)
- DEG=25: 33M samples (~standard)
- DEG=28: 268M samples (~high precision)

### CPU Performance
- Uses OpenMP for parallelization
- Automatically detects available CPU cores
- Typical: 200K-500K queries/second per core

### GPU Performance
- Significantly faster for large DEG values (≥25)
- Requires CUDA-capable GPU
- Configure `CUDA_BLOCKS` and `CUDA_THREADS` based on GPU

## Make Targets

```
make                 - Build CPU binary
make keyrecovery     - Build CPU binary (explicit)
make keyrecovery_gpu - Build GPU binary (requires CUDA)
make keyanalysis     - Build C++ analysis tool
make run             - Run CPU attack with current config
make run-cpu         - Same as run
make run-gpu         - Run GPU attack with current config
make test            - Quick test with DEG=8
make info            - Show current configuration
make perf-calc       - Show performance estimates
make viz             - Visualize results (requires NAME=...)
make viz-all         - Visualize all experiments
make viz-compare     - Compare all experiments
make clean           - Remove binaries
make deepclean       - Remove binaries and all data
make help            - Show all available targets
```

## Custom Configurations

To create a custom attack configuration:

1. Copy an existing configuration block in the Makefile
2. Modify the parameters for your attack scenario
3. Update the NAME variable to identify your configuration
4. Comment out other configurations
5. Run `make run`

## Tips

- Start with `make test` to verify the setup works
- Use `make info` to check current configuration before running
- Monitor memory usage for large DEG values (≥30)
- GPU version is most beneficial for DEG ≥ 25
- Use `make perf-calc` to estimate runtime before starting

## License

Distributed under the GNU General Public License v3.0. Refer to the `LICENSE` file for the full text.

## Contact

For questions, issues, or contributions, please contact:
- **Hosein Hadipour**: [hsn.hadipour@gmail.com](mailto:hsn.hadipour@gmail.com)
- **GitHub Issues**: [https://github.com/hadipourh/universalattacks](https://github.com/hadipourh/universalattacks)

## Acknowledgments

- Orthros PRF specification and design
- OpenMP community and runtime maintainers
- NVIDIA CUDA platform for GPU acceleration
