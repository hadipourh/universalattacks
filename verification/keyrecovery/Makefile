# Define the C and C++ compilers
CC = gcc
CXX = g++
NVCC = nvcc

# Define C source and object files
C_SRC = keyrecovery.c
TARGET = keyrecovery
EXEC = keyrecovery
GPU_TARGET = keyrecovery_gpu
GPU_SUPPORT_SRC = gpu_support.c
GPU_SUPPORT_OBJ = gpu_support.o

# Additional source files
UTILS_SRC = utils.c attackutils.c

# Detect platform and set OpenMP flags accordingly
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS: Use Homebrew libomp if available, otherwise warn
    OPENMP_FLAGS = -Xpreprocessor -fopenmp -lomp
    # Try to find Homebrew libomp installation
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(LIBOMP_PREFIX),)
        OPENMP_FLAGS += -I$(LIBOMP_PREFIX)/include -L$(LIBOMP_PREFIX)/lib
    endif
else
    # Linux: Standard OpenMP flag
    OPENMP_FLAGS = -fopenmp
endif

# Base compiler flags
CFLAGS  = -lm -O3 -march=native -Wall $(OPENMP_FLAGS)
HOST_ONLY_CFLAGS = -O3 -march=native -Wall
NVCC_FLAGS = -O3 -std=c++17 -lcurand
CXXFLAGS = -O3 -march=native -std=c++17

# Try to detect Boost headers (Homebrew on macOS) and optional BOOST_ROOT
BOOST_PREFIX := $(shell brew --prefix boost 2>/dev/null)
ifneq ($(BOOST_PREFIX),)
    CXXFLAGS += -I$(BOOST_PREFIX)/include
endif
ifdef BOOST_ROOT
    CXXFLAGS += -I$(BOOST_ROOT)/include
endif
# Common default include prefixes that often contain Boost headers
CXXFLAGS += -I/opt/homebrew/include -I/usr/local/include

CUDA_THREADS ?= 256
CUDA_BLOCKS ?= 256

# Optional: Comment out OPENMP_FLAGS above to disable parallelization


################################################################################################################################################
#  _  _         ____                       _                      _                    ___  
# | || |       |  _ \ ___  _   _ _ __   __| | __   _____ _ __ ___(_) ___  _ __        / _ \ 
# | || |_ _____| |_) / _ \| | | | '_ \ / _` | \ \ / / _ \ '__/ __| |/ _ \| '_ \ _____| | | |
# |__   _|_____|  _ < (_) | |_| | | | | (_| |  \ V /  __/ |  \__ \ | (_) | | | |_____| |_| |
#    |_|       |_| \_\___/ \__,_|_| |_|\__,_|___\_/ \___|_|  |___/_|\___/|_| |_|      \___/ 
#                                          |_____|                                          
# 
# 4-round_version-0
#corr=2^{-2.7}
################################################################################################################################################
# # # Input parameters for python (dy1 for left branch:nibbles from left to right)
#dy2 for rigjht branch-- nibbles from left to right
# dy1 = 0x2 0x8 0x1
# dy2 = 0x4 0x2 0x8

# NAME="R4_V0"
# ##Input parameter for C
# NUM_OF_ROUNDS = 5
# RND_OFFSET = 1

# OUTPUT_MASK = "000000000000000000000eee00000000" #4R version 1
# INPUT_DIFF_PATTERN = "0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0"
# KEY_VAL_TYPE = 0 # 0 for Weak, 1 for random, 2 for Strong
# PLAIN_NIB_VAL_TYPE = 0
# DEG =11 #iteration for each nibble choice
# N = 3 #number of experiments
#Threshold
################################################################################################################################################


################################################################################################################################################
#  _  _         ____                       _                      _                   _ 
# | || |       |  _ \ ___  _   _ _ __   __| | __   _____ _ __ ___(_) ___  _ __       / |
# | || |_ _____| |_) / _ \| | | | '_ \ / _` | \ \ / / _ \ '__/ __| |/ _ \| '_ \ _____| |
# |__   _|_____|  _ < (_) | |_| | | | | (_| |  \ V /  __/ |  \__ \ | (_) | | | |_____| |
#    |_|       |_| \_\___/ \__,_|_| |_|\__,_|___\_/ \___|_|  |___/_|\___/|_| |_|     |_|
#                                          |_____|                                      
# 
#4-round_version-1
#Corr= 2^{-2.3}
################################################################################################################################################
# # Input parameters for python
# dy1 = 0x4 0x2
# dy2 = 0x4 0x4

# NAME="R4_V1_Strong"
# ##Input parameter for C
# NUM_OF_ROUNDS = 5
# RND_OFFSET = 2
# DEG =12 #iteration for each nibble choice

# OUTPUT_MASK = "00000000000088800000000000000000" #4R version 1
# INPUT_DIFF_PATTERN = "0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0"
# KEY_VAL_TYPE = 0 # 0 for Weak, 1 for random, 2 for Strong
# PLAIN_NIB_VAL_TYPE = 0
# N = 1 #number of experiments

# CORR=0.203  #2^{-2.3}
#Threshold
################################################################################################################################################


################################################################################################################################################
#  _  _         ____                       _                      _                  ____  
# | || |       |  _ \ ___  _   _ _ __   __| | __   _____ _ __ ___(_) ___  _ __      |___ \ 
# | || |_ _____| |_) / _ \| | | | '_ \ / _` | \ \ / / _ \ '__/ __| |/ _ \| '_ \ _____ __) |
# |__   _|_____|  _ < (_) | |_| | | | | (_| |  \ V /  __/ |  \__ \ | (_) | | | |_____/ __/ 
#    |_|       |_| \_\___/ \__,_|_| |_|\__,_|___\_/ \___|_|  |___/_|\___/|_| |_|    |_____|
#                                          |_____|                                         
# 
#4-round_version-2
#corr 2^{-5.7}
################################################################################################################################################
# # Input parameters for python
# dy1 = 0x8 0x8
# dy2 = 0x8 0x2
# DY1_HI := $(word 1,$(dy1))
# DY1_LO := $(word 2,$(dy1))
# DY2_HI := $(word 1,$(dy2))
# DY2_LO := $(word 2,$(dy2))

# NAME="R4_V2"
# ##Input parameter for C
# NUM_OF_ROUNDS = 5
# RND_OFFSET = 2

# OUTPUT_MASK = "00000000222000000000000000000000"
# INPUT_DIFF_PATTERN = "0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"
# KEY_VAL_TYPE = 2 # 0 for Weak, 1 for random, 2 for Strong
# PLAIN_NIB_VAL_TYPE = 0
# DEG = 19 #iteration for each nibble choice
# N = 4 #number of experiments
# CORR = 8  # -log_2(correlation); e.g., 8 means correlation = 2^{-8}
################################################################################################################################################



################################################################################################################################################
#  ____        ____                       _                      _                   _ 
# | ___|      |  _ \ ___  _   _ _ __   __| | __   _____ _ __ ___(_) ___  _ __       / |
# |___ \ _____| |_) / _ \| | | | '_ \ / _` | \ \ / / _ \ '__/ __| |/ _ \| '_ \ _____| |
#  ___) |_____|  _ < (_) | |_| | | | | (_| |  \ V /  __/ |  \__ \ | (_) | | | |_____| |
# |____/      |_| \_\___/ \__,_|_| |_|\__,_|___\_/ \___|_|  |___/_|\___/|_| |_|     |_|
#                                         |_____|                                      
# 
#5-round_version-1
################################################################################################################################################
# Input parameters for python-- 5R version 1
dy1 = 0x8 0x2
dy2 = 0x4 0x4
DY1_HI := $(word 1,$(dy1))
DY1_LO := $(word 2,$(dy1))
DY2_HI := $(word 1,$(dy2))
DY2_LO := $(word 2,$(dy2))

NAME="R5_V1"
##Input parameter for C
NUM_OF_ROUNDS = 6
RND_OFFSET = 2

OUTPUT_MASK = "00000000000000000000000000008088"
INPUT_DIFF_PATTERN = "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0"
KEY_VAL_TYPE = 0 # 0 for Weak, 1 for random, 2 for Strong
PLAIN_NIB_VAL_TYPE = 0
DEG = 28 #iteration for each nibble choice
N = 4 #number of experiments
CORR = 11  # -log_2(correlation); e.g., 11 means correlation = 2^{-11}
################################################################################################################################################


# Data directory and visualization script
VISUALIZE_SCRIPT = scripts/visualize_correlations.py
DATA_DIR = data/$(NAME)
OUTPUT_FILE = $(DATA_DIR)/output.txt

# Enhanced run targets with GPU acceleration and custom parameters

# Allow DEG override from command line: make run DEG=18
DEG ?= 28

# CPU-only run target
.PHONY: run
run: $(EXEC)
	@echo "Running on CPU..."
	@echo "Total samples: 2^$(DEG) = $$(python3 -c "print(2**$(DEG))") per plaintext pair"
	./$(EXEC) $(NUM_OF_ROUNDS) $(RND_OFFSET) $(OUTPUT_MASK) $(INPUT_DIFF_PATTERN) $(KEY_VAL_TYPE) $(PLAIN_NIB_VAL_TYPE) $(DEG) $(N) $(NAME) $(DY1_HI) $(DY1_LO) $(DY2_HI) $(DY2_LO)

# CPU-only run (same as run)
run-cpu: $(EXEC)
	@echo "Running on CPU..."
	@echo "Total samples: 2^$(DEG) = $$(python3 -c "print(2**$(DEG))") per plaintext pair"
	./$(EXEC) $(NUM_OF_ROUNDS) $(RND_OFFSET) $(OUTPUT_MASK) $(INPUT_DIFF_PATTERN) $(KEY_VAL_TYPE) $(PLAIN_NIB_VAL_TYPE) $(DEG) $(N) $(NAME) $(DY1_HI) $(DY1_LO) $(DY2_HI) $(DY2_LO)

.PHONY: run-gpu
run-gpu: $(GPU_TARGET)
	@echo "Running on GPU..."
	@echo "Threads per block: $(CUDA_THREADS)"
	@echo "Blocks: $(CUDA_BLOCKS)"
	@echo "Per-thread samples: 2^$(DEG) = $$(python3 -c 'print(2**$(DEG))')"
	@echo "Total GPU samples: $$(python3 -c 'print($(CUDA_BLOCKS)*$(CUDA_THREADS)*(2**$(DEG)))')"
	@CUDA_THREADS_PER_BLOCK=$(CUDA_THREADS) CUDA_BLOCKS=$(CUDA_BLOCKS) ./$(GPU_TARGET) $(NUM_OF_ROUNDS) $(RND_OFFSET) $(OUTPUT_MASK) $(INPUT_DIFF_PATTERN) $(KEY_VAL_TYPE) $(PLAIN_NIB_VAL_TYPE) $(DEG) $(N) $(NAME) $(DY1_HI) $(DY1_LO) $(DY2_HI) $(DY2_LO)

# Quick test with smaller DEG
test: $(EXEC)
	@echo "Quick test with DEG=8 (256 samples per plaintext)..."
	./$(EXEC) $(NUM_OF_ROUNDS) $(RND_OFFSET) $(OUTPUT_MASK) $(INPUT_DIFF_PATTERN) $(KEY_VAL_TYPE) $(PLAIN_NIB_VAL_TYPE) 8 1 $(NAME) $(DY1_HI) $(DY1_LO) $(DY2_HI) $(DY2_LO)

# Compile the C program
all: keyrecovery

# Information and calculation targets
info:
	@echo "=== ORTHROS Key Recovery Configuration ==="
	@echo "Current DEG: $(DEG)"
	@echo "Samples per plaintext pair: 2^$(DEG) = $$(python3 -c "print(2**$(DEG))")"
	@echo "Number of experiments: $(N)"
	@echo "Trail folder: $(NAME)"
	@echo ""
	@echo "=== Usage Examples ==="
	@echo "Standard run:     make run"
	@echo "Custom DEG:       make run DEG=18"
	@echo "Force CPU:        make run-cpu DEG=15"
	@echo "Quick test:       make test"

# Calculate performance estimates
perf-calc:
	@echo "=== Performance Estimates (DEG=$(DEG)) ==="
	@echo "Total samples per plaintext: 2^$(DEG) = $$(echo "2^$(DEG)" | bc)"
	@echo "CPU (200K queries/sec): ~$$(echo "scale=1; 2^$(DEG) / 200000" | bc)s per plaintext"

keyrecovery:
	$(CC) -o $(EXEC) $(C_SRC) $(CFLAGS)

# Build the keyanalysis tool (C++ with Boost headers)
.PHONY: keyanalysis
keyanalysis:
	$(CXX) $(CXXFLAGS) -o keyanalysis keyanalysis.cpp

$(GPU_TARGET): keyrecovery.cu
	$(NVCC) $(NVCC_FLAGS) keyrecovery.cu -o $(GPU_TARGET)

# Visualization targets
viz: $(EXEC)
	@if [ -z "$(NAME)" ]; then echo "Error: NAME variable required. Usage: make viz NAME=weak_key_deg8"; exit 1; fi
	@FILE=$$(find data/$(NAME) -name "out_*.txt" 2>/dev/null | head -1); \
	if [ -z "$$FILE" ]; then echo "Error: No output file found for $(NAME)"; exit 1; fi; \
	echo "Visualizing $$FILE with label $(NAME)..."; \
	python3 $(VISUALIZE_SCRIPT) -d . -f "$$FILE" -l "$(NAME)" \
		--deg "$(DEG)" --rounds "$(NUM_OF_ROUNDS)" --offset "$(RND_OFFSET)" \
		--output-mask "$(OUTPUT_MASK)" --input-pattern "$(INPUT_DIFF_PATTERN)" \
		--trail-name "$(NAME)"

viz-all:
	python3 $(VISUALIZE_SCRIPT) -d . --batch-all

viz-compare:
	@echo "Comparing all generated experiments..."
	@FILES=$$(find data -name "out_*.txt" 2>/dev/null | tr '\n' ' '); \
	if [ -z "$$FILES" ]; then echo "Error: No output files found"; exit 1; fi; \
	python3 $(VISUALIZE_SCRIPT) -d . --compare $$FILES

# Clean up the generated files
clean:
	rm -f $(C_OBJ) $(EXEC) $(GPU_TARGET) $(OUTPUT_FILE)
	rm -f $(GPU_SUPPORT_OBJ) *.o *.pdf *.png

deepclean: clean
	rm -rf data/

# Help target to show available commands
help:
	@printf "Available targets (make <target>):\n"
	@printf "  %-12s %s\n" "all" "Build the CPU binary (default)"
	@printf "  %-12s %s\n" "keyrecovery" "Explicit CPU build alias"
	@printf "  %-12s %s\n" "run" "Run experiment on CPU with current parameters"
	@printf "  %-12s %s\n" "run-cpu" "Same as run (forces CPU path)"
	@printf "  %-12s %s\n" "run-gpu" "Run CUDA experiment (honours CUDA_* overrides)"
	@printf "  %-12s %s\n" "test" "Quick CPU sanity run with DEG=8"
	@printf "  %-12s %s\n" "info" "Show current configuration variables"
	@printf "  %-12s %s\n" "perf-calc" "Estimate runtime / sample counts"
	@printf "  %-12s %s\n" "viz" "Plot correlations for the latest run"
	@printf "  %-12s %s\n" "viz-all" "Batch-visualise every saved experiment"
	@printf "  %-12s %s\n" "viz-compare" "Overlay correlation curves from all outputs"
	@printf "  %-12s %s\n" "clean" "Delete binaries and generated artefacts"
	@printf "  %-12s %s\n" "deepclean" "Remove binaries plus all data/"
	@printf "  %-12s %s\n" "help" "Show this help screen"
	@printf "\nUsage examples:\n"
	@printf "  make run DEG=18\n"
	@printf "  make run-gpu CUDA_BLOCKS=512 CUDA_THREADS=512\n"
	@printf "  make viz NAME=$(NAME)\n"
.PHONY: all keyrecovery run run-cpu run-gpu test info perf-calc viz viz-all viz-compare clean deepclean help
