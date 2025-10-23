/*
 * CUDA Implementation of Orthros Differential-Linear Cryptanalysis Tool (DLCT)
 * 
 * Copyright (C) 2025 Hosein Hadipour
 * Email: hsn.hadipour@gmail.com
 * Date: September 30, 2025
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 * 
 * Features:
 * - GPU-accelerated differential-linear correlation computation for Orthros
 * - Support for reduced-round Orthros (configurable rounds and offset)
 * - Support for left/right branch modes and PRF mode  
 * - Self-test verification against Orthros test vectors
 * - Hardware random number generation (RDRAND)
 * - DLCT table generation for single-bit patterns
 * - Professional statistical analysis and visualization
 * - Command-line interface with offset and mode support
 */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "common_types.h"
// Define constants before including orthros_core.h to avoid extern warnings
#define ORTHROS_DEFINE_CONSTANTS
#include "orthros_core.h"
#include <curand_kernel.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include <string.h>
#include <stdbool.h>
#include <inttypes.h>

// Linux/Unix includes for RDRAND and timing
#include <sys/random.h>
#include <unistd.h>

// Define GRND_NONBLOCK if not available
#ifndef GRND_NONBLOCK
#define GRND_NONBLOCK 0x0001
#endif

// Architecture-specific RDRAND support
#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)) && defined(__RDRND__)
    #include <x86intrin.h>
    #include <cpuid.h>
    #define HAVE_X86_RDRAND 1
#else
    #define HAVE_X86_RDRAND 0
#endif

// Global counters for random source usage tracking
static struct {
    uint64_t rdrand_calls;
    uint64_t getrandom_calls; 
    uint64_t fallback_calls;
    uint64_t total_calls;
} random_stats = {0, 0, 0, 0};

// Structure for experiment key logging
// Global arrays for tracking experiments
static ExperimentKey *experiment_keys = NULL;
static double *correlation_values = NULL;

// External function declarations
void display_correlation_histogram(const double *correlations, int num_experiments, 
                                  Config *config);
void display_experiment_keys(ExperimentKey *keys, int num_experiments, 
                            Config *config);

// ================================
// Configuration Constants
// ================================
#define DEFAULT_BLOCKS          256
#define DEFAULT_THREADS         512
#define DEFAULT_EXPERIMENTS     16
#define SAMPLE_POWER_DEFAULT    18
#define MAX_EXPERIMENTS         1000
#define MAX_SAMPLE_POWER        32
#define DLCT_PRECISION_SCALE    100

// RDRAND constants
#define RDRAND_MASK            0x40000000
#define RETRY_LIMIT            10
#define RDRAND_SUCCESS         1
#define RDRAND_NOT_READY       0

// ================================
// Data Structures
// ================================

typedef struct {
    unsigned char *batch_keys;
    uint64_t *d_total_equal;
    bool *d_test_result;
    unsigned char *d_input_diff_left;
    unsigned char *d_input_diff_right;
    unsigned char *d_output_mask;
} MemoryContext;

// ================================
// Device Constants
// ================================
__device__ __constant__ unsigned char d_batch_round_keys[ORTHROS_BATCH_SIZE * ORTHROS_STATE_SIZE];

// ================================
// Error Handling Macros
// ================================
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(error)); \
        cleanup_memory_context(&memory_ctx); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

#define CHECK_ALLOC(ptr, name) do { \
    if (!(ptr)) { \
        fprintf(stderr, "Memory allocation failed for %s at %s:%d\n", \
                (name), __FILE__, __LINE__); \
        cleanup_memory_context(&memory_ctx); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

// Global memory context for cleanup
static MemoryContext memory_ctx = {0};

// ================================
// Memory Management Functions
// ================================

/**
 * @brief Initialize memory context with safe defaults
 */
void init_memory_context(MemoryContext *ctx) {
    memset(ctx, 0, sizeof(MemoryContext));
}

/**
 * @brief Cleanup all allocated memory safely
 */
void cleanup_memory_context(MemoryContext *ctx) {
    if (ctx->batch_keys) {
        free(ctx->batch_keys);
        ctx->batch_keys = NULL;
    }
    if (ctx->d_total_equal) {
        cudaFree(ctx->d_total_equal);
        ctx->d_total_equal = NULL;
    }
    if (ctx->d_test_result) {
        cudaFree(ctx->d_test_result);
        ctx->d_test_result = NULL;
    }
    if (ctx->d_input_diff_left) {
        cudaFree(ctx->d_input_diff_left);
        ctx->d_input_diff_left = NULL;
    }
    if (ctx->d_input_diff_right) {
        cudaFree(ctx->d_input_diff_right);
        ctx->d_input_diff_right = NULL;
    }
    if (ctx->d_output_mask) {
        cudaFree(ctx->d_output_mask);
        ctx->d_output_mask = NULL;
    }
}

// ================================
// Memory Management
// ================================

/**
 * @brief Update device memory with new config data
 */
void update_device_config(MemoryContext *ctx, const Config *config) {
    CUDA_CHECK(cudaMemcpy(ctx->d_input_diff_left, config->input_diff_left, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_input_diff_right, config->input_diff_right, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_output_mask, config->output_mask, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
}

/**
 * @brief Allocate device memory for DLCT computation  
 */
void allocate_correlation_memory(MemoryContext *ctx, const Config *config) {
    const size_t batch_size = ORTHROS_BATCH_SIZE * ORTHROS_STATE_SIZE * sizeof(unsigned char);

    ctx->batch_keys = (unsigned char*)malloc(batch_size);
    CHECK_ALLOC(ctx->batch_keys, "batch keys");

    CUDA_CHECK(cudaMalloc(&ctx->d_total_equal, sizeof(uint64_t)));
    CUDA_CHECK(cudaMalloc(&ctx->d_input_diff_left, ORTHROS_STATE_SIZE * sizeof(unsigned char)));
    CUDA_CHECK(cudaMalloc(&ctx->d_input_diff_right, ORTHROS_STATE_SIZE * sizeof(unsigned char)));
    CUDA_CHECK(cudaMalloc(&ctx->d_output_mask, ORTHROS_STATE_SIZE * sizeof(unsigned char)));
    
    // Copy the input arrays to device memory  
    CUDA_CHECK(cudaMemcpy(ctx->d_input_diff_left, config->input_diff_left, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_input_diff_right, config->input_diff_right, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(ctx->d_output_mask, config->output_mask, 
                         ORTHROS_STATE_SIZE * sizeof(unsigned char), cudaMemcpyHostToDevice));
}

// ================================
// Device Functions
// ================================

/**
 * @brief Self-test using Orthros test vectors
 */
__device__ bool orthros_self_test_device() {
    const unsigned char test_pt[ORTHROS_STATE_SIZE] = {
        0xa,0x9,0x4,0x7,0x4,0x3,0x6,0x7,0x1,0x0,0x9,0x2,0x4,0xc,0xc,0xd,
        0x4,0x7,0xf,0x2,0xd,0x5,0x7,0x1,0xd,0xe,0xe,0xa,0x8,0xf,0x0,0x5
    };
    const unsigned char test_key[ORTHROS_STATE_SIZE] = {
        0x4,0xa,0x2,0xb,0xe,0x6,0x0,0xe,0x3,0xd,0xb,0x6,0xa,0xb,0xe,0x0,
        0xc,0x0,0x3,0xe,0xa,0xe,0xc,0x6,0x6,0xf,0xd,0x0,0x5,0xd,0x0,0xc
    };
    const unsigned char expected_ct[ORTHROS_STATE_SIZE] = {
        0xf,0xb,0xf,0xd,0x0,0xb,0x2,0xd,0xc,0x5,0x8,0xc,0x8,0x3,0x5,0x2,
        0x7,0x4,0x4,0xb,0x9,0x7,0x6,0x1,0x9,0x1,0x7,0x6,0xe,0x8,0xd,0x4
    };
    
    unsigned char result[ORTHROS_STATE_SIZE];
    orthros_encrypt(test_pt, test_pt, result, test_key, 0, ORTHROS_MAX_ROUNDS, ORTHROS_MODE_LEFT);
    
    for (int i = 0; i < ORTHROS_STATE_SIZE; i++) {
        if (result[i] != expected_ct[i]) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Main DLCT computation kernel
 */
__device__ __forceinline__ void next_plaintext_state(curandStatePhilox4_32_10_t &state,
                                                     unsigned char *pt_left,
                                                     unsigned char *pt_right) {
    // Generate random state for both branches (32 nibbles each)
    for (int i = 0; i < ORTHROS_STATE_SIZE; i += 4) {
        uint4 rnd = curand4(&state);
        pt_left[i] = rnd.x & 0xF;
        pt_left[i+1] = rnd.y & 0xF;
        pt_left[i+2] = rnd.z & 0xF;
        pt_left[i+3] = rnd.w & 0xF;
    }
    for (int i = 0; i < ORTHROS_STATE_SIZE; i += 4) {
        uint4 rnd = curand4(&state);
        pt_right[i] = rnd.x & 0xF;
        pt_right[i+1] = rnd.y & 0xF;
        pt_right[i+2] = rnd.z & 0xF;
        pt_right[i+3] = rnd.w & 0xF;
    }
}

__global__ void orthros_dlct_kernel(
    const unsigned char *input_diff_left,
    const unsigned char *input_diff_right, 
    const unsigned char *output_mask,
    int rounds,
    int offset,
    int mode,
    uint64_t samples_per_thread,
    uint64_t *total_equal,
    uint64_t global_seed,
    int experiment_id,
    int total_threads
) {
    extern __shared__ uint64_t sdata[];

    const int lane = threadIdx.x;
    const int tid = blockIdx.x * blockDim.x + lane;
    
    // Ensure experiment_id is within bounds for the batch keys
    const int safe_exp_id = experiment_id % ORTHROS_BATCH_SIZE;
    const int key_offset = safe_exp_id * ORTHROS_STATE_SIZE;

    curandStatePhilox4_32_10_t rng_state;
    const unsigned long long sequence = static_cast<unsigned long long>(experiment_id) *
                                        static_cast<unsigned long long>(total_threads) +
                                        static_cast<unsigned long long>(tid);
    curand_init(static_cast<unsigned long long>(global_seed), sequence, 0ULL, &rng_state);

    unsigned char pt1_left[ORTHROS_STATE_SIZE], pt1_right[ORTHROS_STATE_SIZE];
    unsigned char pt2_left[ORTHROS_STATE_SIZE], pt2_right[ORTHROS_STATE_SIZE];
    unsigned char ct1[ORTHROS_STATE_SIZE], ct2[ORTHROS_STATE_SIZE];
    
    next_plaintext_state(rng_state, pt1_left, pt1_right);
    uint64_t counter = 0;

    for (uint64_t i = 0; i < samples_per_thread; ++i) {
        // Apply differential
        for (int j = 0; j < ORTHROS_STATE_SIZE; j++) {
            pt2_left[j] = pt1_left[j] ^ input_diff_left[j];
            pt2_right[j] = pt1_right[j] ^ input_diff_right[j];
        }
        
        // Encrypt both plaintexts
        orthros_encrypt(pt1_left, pt1_right, ct1, &d_batch_round_keys[key_offset], offset, rounds, static_cast<OrthrosMode>(mode));
        orthros_encrypt(pt2_left, pt2_right, ct2, &d_batch_round_keys[key_offset], offset, rounds, static_cast<OrthrosMode>(mode));
        
        // Compute linear approximation: compare dot products like OpenMP version
        const int parity_base = dot_product(ct1, output_mask);
        const int parity_diff = dot_product(ct2, output_mask);
        if (parity_base == parity_diff) {
            counter++;
        }
        
        next_plaintext_state(rng_state, pt1_left, pt1_right);
    }

    sdata[lane] = counter;
    __syncthreads();

    for (unsigned int stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
        if (lane < stride) {
            sdata[lane] += sdata[lane + stride];
        }
        __syncthreads();
    }

    if (lane == 0) {
        atomicAdd(reinterpret_cast<unsigned long long*>(total_equal),
                  static_cast<unsigned long long>(sdata[0]));
    }
}

/**
 * @brief Self-test kernel wrapper
 */
__global__ void orthros_self_test_kernel(bool *result) {
    *result = orthros_self_test_device();
}

// ================================
// Random Number Generation
// ================================

#if HAVE_X86_RDRAND
/**
 * @brief Check RDRAND support via CPUID
 */
static int check_rdrand_support() {
    uint32_t eax, ebx, ecx, edx;
    __cpuid(0, eax, ebx, ecx, edx);
    if (ebx != 0x756e6547 || edx != 0x49656e69 || ecx != 0x6c65746e) {
        return 0; // Not genuine Intel
    }
    __cpuid(1, eax, ebx, ecx, edx);
    return (ecx & RDRAND_MASK) ? 1 : 0;
}

/**
 * @brief Get random 64-bit value using RDRAND with tracking
 */
static int rdrand_64_impl(uint64_t *x) {
    static int supported = -1;
    if (supported < 0) {
        supported = check_rdrand_support();
    }
    
    random_stats.total_calls++;
    
    if (supported) {
        for (int i = 0; i < RETRY_LIMIT; i++) {
            unsigned long long tmp;
            if (_rdrand64_step(&tmp)) {
                *x = (uint64_t)tmp;
                random_stats.rdrand_calls++;
                return RDRAND_SUCCESS;
            }
        }
    }
    
    // Fallback to system random
    if (getrandom(x, sizeof(*x), 0) == sizeof(*x)) {
        random_stats.getrandom_calls++;
        return RDRAND_SUCCESS;
    }
    
    // Final fallback
    const uint64_t time_mix = ((uint64_t)time(NULL) << 32) ^ (uint64_t)clock();
    *x = time_mix ^ 0xA5A5A5A5A5A5A5A5ULL;
    random_stats.fallback_calls++;
    return RDRAND_SUCCESS;
}
#else
/**
 * @brief Fallback random generation for non-x86 platforms with tracking
 */
static int rdrand_64_impl(uint64_t *x) {
    random_stats.total_calls++;
    
    if (getrandom(x, sizeof(*x), 0) == sizeof(*x)) {
        random_stats.getrandom_calls++;
        return RDRAND_SUCCESS;
    }
    
    // Time-based fallback
    const uint64_t time_mix = ((uint64_t)time(NULL) << 32) ^ (uint64_t)clock();
    *x = time_mix ^ 0xA5A5A5A5A5A5A5A5ULL;
    random_stats.fallback_calls++;
    return RDRAND_SUCCESS;
}
#endif

/**
 * @brief Safe random number generation with fallback
 */
static bool get_random_uint64(uint64_t *value) {
    return rdrand_64_impl(value) == RDRAND_SUCCESS;
}

/**
 * @brief Generate random state for Orthros
 */
static void get_random_state(unsigned char *state) {
    uint64_t rnd[4];
    for (int i = 0; i < 4; i++) {
        if (!get_random_uint64(&rnd[i])) {
            rnd[i] = ((uint64_t)time(NULL) << 32) ^ ((uint64_t)(i + 1) * 0x9E3779B97F4A7C15ULL);
        }
    }
    
    // Convert to nibbles
    for (int i = 0; i < ORTHROS_STATE_SIZE; i++) {
        state[i] = (rnd[i / 16] >> ((i % 16) * 4)) & 0xF;
    }
}

/**
 * @brief Detect and report the specific random source being used
 */
static void report_random_source() {
#if HAVE_X86_RDRAND
    static int rdrand_available = -1;
    if (rdrand_available < 0) {
        rdrand_available = check_rdrand_support();
    }
    
    if (rdrand_available) {
        printf("Random source: CPU RDRAND (hardware entropy)\n");
        printf("  - Primary: Intel RDRAND instruction\n");
        printf("  - Fallback 1: getrandom() system call\n");
        printf("  - Fallback 2: time-based mixing\n");
    } else {
        printf("Random source: System entropy (software)\n");
        printf("  - Primary: getrandom() system call\n");
        printf("  - Fallback: time-based mixing\n");
        printf("  - Note: RDRAND not available on this CPU\n");
    }
#else
    // Test if getrandom is available
    uint64_t test_value;
    if (getrandom(&test_value, sizeof(test_value), GRND_NONBLOCK) == sizeof(test_value)) {
        printf("Random source: System entropy (software)\n");
        printf("  - Primary: getrandom() system call (kernel entropy pool)\n");
        printf("  - Fallback: time-based mixing\n");
        printf("  - Note: Non-x86 platform, RDRAND not supported\n");
    } else {
        printf("Random source: Time-based mixing (software fallback)\n");
        printf("  - Primary: time()/clock() mixing with constants\n");
        printf("  - Note: getrandom() not available or blocked\n");
        printf("  - Warning: Reduced entropy quality\n");
    }
#endif
    printf("  - Status: Testing actual usage during execution...\n");
    printf("\n");
}

/**
 * @brief Report actual random source usage statistics
 */
static void report_random_usage_stats() {
    if (random_stats.total_calls == 0) {
        return; // No random calls made
    }
    
    printf("Random Source Usage Statistics:\n");
    printf("================================\n");
    printf("Total random calls: %lu\n", random_stats.total_calls);
    
    if (random_stats.rdrand_calls > 0) {
        printf("RDRAND (hardware): %lu calls (%.1f%%)\n", 
               random_stats.rdrand_calls,
               100.0 * random_stats.rdrand_calls / random_stats.total_calls);
    }
    
    if (random_stats.getrandom_calls > 0) {
        printf("getrandom() (system): %lu calls (%.1f%%)\n", 
               random_stats.getrandom_calls,
               100.0 * random_stats.getrandom_calls / random_stats.total_calls);
    }
    
    if (random_stats.fallback_calls > 0) {
        printf("Time-based fallback: %lu calls (%.1f%%)\n", 
               random_stats.fallback_calls,
               100.0 * random_stats.fallback_calls / random_stats.total_calls);
    }
    
    // Report actual primary source used
    if (random_stats.rdrand_calls > 0) {
        printf("✓ Verified: RDRAND hardware entropy is working\n");
    } else if (random_stats.getrandom_calls > 0) {
        printf("✓ Verified: getrandom() system entropy is working\n");
    } else if (random_stats.fallback_calls > 0) {
        printf("⚠ Warning: Only time-based fallback was used\n");
    }
    
    printf("\n");
}

// ================================
// Timing Functions
// ================================

static clock_t start_time;

static void start_timer() {
    start_time = clock();
}

static double get_time_ms() {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC * 1000.0;
}

// ================================
// CUDA Operations
// ================================

/**
 * @brief Initialize CUDA device and constants
 */
static void initialize_cuda_device(const Config *config) {
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        fprintf(stderr, "Error: No CUDA devices found\n");
        exit(EXIT_FAILURE);
    }
    
    cudaDeviceProp props;
    CUDA_CHECK(cudaGetDeviceProperties(&props, 0));
    
    if (config->verbose) {
        printf("Using GPU: %s\n", props.name);
        printf("Compute capability: %d.%d\n", props.major, props.minor);
        printf("Multiprocessors: %d\n", props.multiProcessorCount);
        printf("Max threads per block: %d\n\n", props.maxThreadsPerBlock);
    }
    
    // Validate thread count
    if (config->threads > props.maxThreadsPerBlock) {
        fprintf(stderr, "Error: %d threads exceeds GPU capability (%d max).\n", 
                config->threads, props.maxThreadsPerBlock);
        exit(EXIT_FAILURE);
    }
    
    // Set device and initialize constants
    CUDA_CHECK(cudaSetDevice(0));
}

/**
 * @brief Run self-test on device
 */
static bool run_device_self_test() {
    // Allocate device memory for result
    CUDA_CHECK(cudaMalloc(&memory_ctx.d_test_result, sizeof(bool)));
    
    // Launch self-test kernel
    orthros_self_test_kernel<<<1, 1>>>(memory_ctx.d_test_result);
    CUDA_CHECK(cudaDeviceSynchronize());
    
    // Get result
    bool result = false;
    CUDA_CHECK(cudaMemcpy(&result, memory_ctx.d_test_result, sizeof(bool), cudaMemcpyDeviceToHost));
    
    return result;
}

// ================================
// Key Generation
// ================================

/**
 * @brief Generate batch of random keys with logging
 */
static void generate_experiment_keys(int num_experiments, int batch_start_id) {
    for (int i = 0; i < num_experiments; i++) {
        unsigned char key[ORTHROS_STATE_SIZE];
        get_random_state(key);
        
        // Log the key for this experiment
        int global_exp_id = batch_start_id + i;
        if (experiment_keys && global_exp_id < MAX_EXPERIMENTS) {
            experiment_keys[global_exp_id].experiment_id = global_exp_id;
            memcpy(experiment_keys[global_exp_id].key, key, ORTHROS_STATE_SIZE);
            experiment_keys[global_exp_id].correlation = 0.0; // Will be filled later
        }
        
        // Store key in batch
        memcpy(memory_ctx.batch_keys + i * ORTHROS_STATE_SIZE, key, ORTHROS_STATE_SIZE);
    }
    
    // Transfer to device constant memory
    const size_t transfer_size = num_experiments * ORTHROS_STATE_SIZE * sizeof(unsigned char);
    CUDA_CHECK(cudaMemcpyToSymbol(d_batch_round_keys, memory_ctx.batch_keys, transfer_size));
}

// ================================
// Correlation Computation
// ================================

/**
 * @brief Process a single correlation experiment
 */
static double process_single_experiment(const Config *config, int exp_id, int batch_relative_id, 
                                        int total_threads, uint64_t samples_per_thread) {
    uint64_t global_seed;
    if (!get_random_uint64(&global_seed)) {
        global_seed = ((uint64_t)time(NULL) << 32) ^ ((uint64_t)exp_id * 0x9E3779B97F4A7C15ULL);
    }

    CUDA_CHECK(cudaMemset(memory_ctx.d_total_equal, 0, sizeof(uint64_t)));

    start_timer();
    const size_t shared_bytes = static_cast<size_t>(config->threads) * sizeof(uint64_t);
    
    orthros_dlct_kernel<<<config->blocks, config->threads, shared_bytes>>>(
        memory_ctx.d_input_diff_left,
        memory_ctx.d_input_diff_right,
        memory_ctx.d_output_mask,
        config->rounds,
        config->offset,
        config->mode,
        samples_per_thread,
        memory_ctx.d_total_equal,
        global_seed,
        batch_relative_id,
        total_threads
    );
    CUDA_CHECK(cudaDeviceSynchronize());
    const double kernel_time = get_time_ms();

    uint64_t total_equal = 0;
    CUDA_CHECK(cudaMemcpy(&total_equal, memory_ctx.d_total_equal, sizeof(uint64_t), cudaMemcpyDeviceToHost));

    const uint64_t total_samples = static_cast<uint64_t>(config->blocks) * config->threads * samples_per_thread;
    const double correlation = (2.0 * total_equal) / total_samples - 1.0;

    if (config->verbose) {
        printf("Experiment %d: %.8f (%.2f ms, %lu/%lu equal)\n", 
               exp_id, correlation, kernel_time, total_equal, total_samples);
    }

    return correlation;
}

// ================================
// Main Correlation Computation
// ================================

/**
 * @brief Compute averaged correlation across multiple experiments
 */
static double compute_correlation(Config config) {
    // Check if memory is already allocated
    static bool memory_allocated = false;
    if (!memory_allocated) {
        allocate_correlation_memory(&memory_ctx, &config);
        memory_allocated = true;
    } else {
        // Update device memory with new config
        update_device_config(&memory_ctx, &config);
    }
    
    // Allocate tracking arrays
    experiment_keys = (ExperimentKey*)malloc(config.experiments * sizeof(ExperimentKey));
    correlation_values = (double*)malloc(config.experiments * sizeof(double));
    CHECK_ALLOC(experiment_keys, "experiment keys");
    CHECK_ALLOC(correlation_values, "correlation values");

    const uint64_t samples_per_thread = 1ULL << config.sample_power;
    const uint64_t total_samples = static_cast<uint64_t>(config.blocks) * config.threads * samples_per_thread;
    const int total_threads = config.blocks * config.threads;

    if (config.verbose) {
        printf("Starting %d experiments with 2^%d samples per thread\n",
               config.experiments, config.sample_power);
        printf("Total samples per experiment: %lu (%.2e)\n", 
               total_samples, (double)total_samples);
        printf("GPU configuration: %d blocks × %d threads = %d total\n\n",
               config.blocks, config.threads, total_threads);
    }

    double sum_correlations = 0.0;
    const int batch_size = (config.experiments < ORTHROS_BATCH_SIZE) ? config.experiments : ORTHROS_BATCH_SIZE;
    
    for (int batch_start = 0; batch_start < config.experiments; batch_start += batch_size) {
        const int current_batch_size = (batch_start + batch_size > config.experiments) ? 
                                       (config.experiments - batch_start) : batch_size;
        
        generate_experiment_keys(current_batch_size, batch_start);
        
        for (int i = 0; i < current_batch_size; i++) {
            const int global_exp_id = batch_start + i;
            const double correlation = process_single_experiment(&config, global_exp_id, i, 
                                                                 total_threads, samples_per_thread);
            sum_correlations += correlation;
            correlation_values[global_exp_id] = correlation;
            
            if (experiment_keys) {
                experiment_keys[global_exp_id].correlation = correlation;
            }
        }
    }

    // Display results
    if (config.verbose) {
        printf("\nExperiment Summary:\n");
        printf("===================\n");
        display_correlation_histogram(correlation_values, config.experiments, &config);
    }
    
    display_experiment_keys(experiment_keys, config.experiments, &config);
    
    // Cleanup tracking arrays
    if (experiment_keys) {
        free(experiment_keys);
        experiment_keys = NULL;
    }
    if (correlation_values) {
        free(correlation_values);
        correlation_values = NULL;
    }
    
    return sum_correlations / config.experiments;
}

// ================================
// DLCT Table Generation
// ================================

/**
 * @brief Generate full DLCT table for single-bit patterns
 */
static void generate_dlct_table(Config config) {
    printf("Generating DLCT table for %d-round Orthros (offset %d, mode %d)\n", 
           config.rounds, config.offset, static_cast<int>(config.mode));
    printf("Configuration: %d blocks × %d threads, 2^%d samples, %d experiments\n\n",
           config.blocks, config.threads, config.sample_power, config.experiments);
    
    // Prepare output files
    char csv_filename[256], dzn_filename[256];
    snprintf(csv_filename, sizeof(csv_filename), "orthros_%dr_offset%d_mode%d_dlct.csv", 
             config.rounds, config.offset, static_cast<int>(config.mode));
    snprintf(dzn_filename, sizeof(dzn_filename), "orthros_%dr_offset%d_mode%d_dlct.dzn", 
             config.rounds, config.offset, static_cast<int>(config.mode));
    
    FILE *csv_file = fopen(csv_filename, "w");
    FILE *dzn_file = fopen(dzn_filename, "w");
    
    if (!csv_file || !dzn_file) {
        fprintf(stderr, "Error: Could not create output files\n");
        if (csv_file) fclose(csv_file);
        if (dzn_file) fclose(dzn_file);
        return;
    }
    
    printf("Generating: %s (CSV format)\n", csv_filename);
    printf("Generating: %s (MiniZinc format)\n\n", dzn_filename);
    
    // Write headers
    fprintf(csv_file, "Input_Diff_Left_Bit,Input_Diff_Right_Bit,Output_Mask_Bit,Correlation_Pow2,Log2_Magnitude\n");
    fprintf(dzn_file, "%% Orthros %d-round DLCT (offset %d, mode %d, single-bit patterns)\n", 
            config.rounds, config.offset, static_cast<int>(config.mode));
    fprintf(dzn_file, "%% rounds=%d, offset=%d, mode=%d, blocks=%d, threads=%d, samples_per_thread=2^%d, total_samples=%llu, experiments=%d\n", 
        config.rounds, config.offset, static_cast<int>(config.mode), config.blocks, config.threads, config.sample_power, 
        (unsigned long long)config.blocks * config.threads * (1ULL << config.sample_power), config.experiments);
    
    // Storage for DLCT weights
    int dlct_weights[ORTHROS_STATE_SIZE][ORTHROS_STATE_SIZE][ORTHROS_STATE_SIZE];
    int max_weight = 0;
    
    printf("Progress: ");
    fflush(stdout);
    
    int progress_counter = 0;
    const int total_entries = ORTHROS_STATE_SIZE * ORTHROS_STATE_SIZE * ORTHROS_STATE_SIZE;
    
    // Generate DLCT entries for single-bit patterns
    for (int i = 0; i < ORTHROS_STATE_SIZE; i++) {
        for (int j = 0; j < ORTHROS_STATE_SIZE; j++) {
            for (int k = 0; k < ORTHROS_STATE_SIZE; k++) {
                if (progress_counter % (total_entries / 64) == 0) {
                    printf(".");
                    fflush(stdout);
                }
                progress_counter++;
                
                // Single-bit patterns
                Config temp_config = config;
                memset(temp_config.input_diff_left, 0, ORTHROS_STATE_SIZE);
                memset(temp_config.input_diff_right, 0, ORTHROS_STATE_SIZE);
                memset(temp_config.output_mask, 0, ORTHROS_STATE_SIZE);
                
                temp_config.input_diff_left[i] = 1;
                temp_config.input_diff_right[j] = 1;
                temp_config.output_mask[k] = 1;
                temp_config.verbose = false;
                
                const double correlation = compute_correlation(temp_config);
                const double corr_mag = fabs(correlation);
                
                int weight;
                if (corr_mag == 0.0) {
                    fprintf(csv_file, "%d,%d,%d,+2^(-inf),inf\n", i, j, k);
                    weight = -1;
                } else {
                    const double exponent = log2(corr_mag);
                    const char sign_char = (correlation >= 0.0) ? '+' : '-';
                    fprintf(csv_file, "%d,%d,%d,%c2^(%.2f),%.6f\n",
                            i, j, k, sign_char, exponent, -exponent);
                    
                    const double scaled = DLCT_PRECISION_SCALE * (-exponent);
                    weight = (int)round(scaled);
                    if (weight > max_weight) {
                        max_weight = weight;
                    }
                }
                
                dlct_weights[i][j][k] = weight;
            }
        }
    }
    
    printf(" done!\n\n");
    
    // Write MiniZinc format
    fprintf(dzn_file, "\n%% DLCT weights matrix [left_bit][right_bit][mask_bit]\n");
    fprintf(dzn_file, "array[0..%d, 0..%d, 0..%d] of int: dlct_weights = array3d(0..%d, 0..%d, 0..%d, [\n",
            ORTHROS_STATE_SIZE-1, ORTHROS_STATE_SIZE-1, ORTHROS_STATE_SIZE-1,
            ORTHROS_STATE_SIZE-1, ORTHROS_STATE_SIZE-1, ORTHROS_STATE_SIZE-1);
    
    for (int i = 0; i < ORTHROS_STATE_SIZE; i++) {
        for (int j = 0; j < ORTHROS_STATE_SIZE; j++) {
            for (int k = 0; k < ORTHROS_STATE_SIZE; k++) {
                fprintf(dzn_file, "%d", dlct_weights[i][j][k]);
                if (i < ORTHROS_STATE_SIZE-1 || j < ORTHROS_STATE_SIZE-1 || k < ORTHROS_STATE_SIZE-1) {
                    fprintf(dzn_file, ",");
                }
                if ((i * ORTHROS_STATE_SIZE * ORTHROS_STATE_SIZE + j * ORTHROS_STATE_SIZE + k + 1) % 8 == 0) {
                    fprintf(dzn_file, "\n");
                }
            }
        }
    }
    fprintf(dzn_file, "]);\n");
    fprintf(dzn_file, "\n%% Maximum weight found: %d\n", max_weight);
    
    fclose(csv_file);
    fclose(dzn_file);
    
    printf("DLCT table generation completed!\n");
    printf("Maximum weight: %d\n", max_weight);
    printf("Files written: %s, %s\n", csv_filename, dzn_filename);
}

// ================================
// Utility Functions for Hex Parsing
// ================================

static int hex_char_to_nibble(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static bool parse_hex_nibbles(const char *hex, unsigned char *dst) {
    const size_t len = strlen(hex);
    if (len != ORTHROS_STATE_SIZE) {
        return false;
    }
    for (size_t i = 0; i < len; ++i) {
        const int nib = hex_char_to_nibble(hex[i]);
        if (nib < 0) {
            return false;
        }
        dst[i] = static_cast<unsigned char>(nib & 0xF);
    }
    return true;
}

static void format_nibbles(const unsigned char *src, char *dst) {
    static const char HEX_DIGITS[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
    for (int i = 0; i < ORTHROS_STATE_SIZE; ++i) {
        dst[i] = HEX_DIGITS[src[i] & 0xF];
    }
    dst[ORTHROS_STATE_SIZE] = '\0';
}

// ================================
// Wildcard Mask Expansion Functions
// ================================

static int count_wildcards(const char *pattern) {
    int count = 0;
    for (size_t i = 0; i < strlen(pattern); i++) {
        if (pattern[i] == 'x' || pattern[i] == 'X') {
            count++;
        }
    }
    return count;
}

static void expand_wildcard_mask(const char *pattern, int hex_value, unsigned char *output_mask) {
    for (size_t i = 0; i < strlen(pattern); i++) {
        if (pattern[i] == 'x' || pattern[i] == 'X') {
            // All wildcards get the same hex value
            output_mask[i] = static_cast<unsigned char>(hex_value & 0xF);
        } else {
            // Use the fixed character from pattern
            int nib = hex_char_to_nibble(pattern[i]);
            if (nib >= 0) {
                output_mask[i] = static_cast<unsigned char>(nib & 0xF);
            } else {
                output_mask[i] = 0;
            }
        }
    }
}

static void run_wildcard_experiments(Config base_config) {
    const int num_wildcards = count_wildcards(base_config.wildcard_mask_pattern);
    const int total_combinations = 15;  // 15 hex values: 1-f (skip 0)
    
    printf("Wildcard mask analysis: %s\n", base_config.wildcard_mask_pattern);
    printf("Number of wildcards: %d\n", num_wildcards);
    printf("Total combinations: %d (all wildcards take same value: 1-f)\n", total_combinations);
    printf("Configuration:\n");
    printf("  Blocks: %d\n", base_config.blocks);
    printf("  Threads per block: %d\n", base_config.threads);
    printf("  Total threads: %d × %d = %d\n", base_config.blocks, base_config.threads, base_config.blocks * base_config.threads);
    printf("  Samples per thread: 2^%d = %d\n", base_config.sample_power, 1 << base_config.sample_power);
    printf("  Total samples: %lld ≈ 2^%.1f\n", 
           (long long)base_config.blocks * base_config.threads * (1LL << base_config.sample_power),
           log2((double)base_config.blocks * base_config.threads * (1LL << base_config.sample_power)));
    printf("  Experiments per combination: %d\n\n", base_config.experiments);
    
    // Store results for final summary
    double correlations[15];
    char masks[15][ORTHROS_STATE_SIZE + 1];
    
    for (int hex_val = 1; hex_val <= 15; hex_val++) {
        Config config = base_config;
        
        // Generate the specific mask for this hex value
        expand_wildcard_mask(base_config.wildcard_mask_pattern, hex_val, config.output_mask);
        
        // Format the mask for display
        char mask_hex[ORTHROS_STATE_SIZE + 1];
        format_nibbles(config.output_mask, mask_hex);
        strcpy(masks[hex_val - 1], mask_hex);
        
        printf("=== Combination %d/15: Hex=%x, Mask %s ===\n", hex_val, hex_val, mask_hex);
        
        // Run the experiment with this specific mask
        const double correlation = compute_correlation(config);
        const double corr_mag = fabs(correlation);
        correlations[hex_val - 1] = correlation;
        
        printf("Correlation: %+.12e (2^%.4f)\n", correlation, log2(corr_mag));
        printf("Security assessment: ");
        if (corr_mag > 0.1) {
            printf("BROKEN - Strong correlation detected\n");
        } else if (corr_mag > 0.01) {
            printf("WEAK - Significant correlation detected\n");
        } else if (corr_mag > 0.001) {
            printf("SUSPICIOUS - Detectable correlation\n");
        } else {
            printf("SECURE - No significant correlation\n");
        }
        printf("\n");
    }
    
    // Final summary table
    printf("==========================================================================\n");
    printf("WILDCARD ANALYSIS SUMMARY\n");
    printf("==========================================================================\n");
    printf("Pattern: %s\n", base_config.wildcard_mask_pattern);
    printf("┌─────┬─────┬────────────────────────────────┬─────────────────────┬─────────────┐\n");
    printf("│ Hex │ Dec │              Mask              │    Correlation      │   Log₂ |C|   │\n");
    printf("├─────┼─────┼────────────────────────────────┼─────────────────────┼─────────────┤\n");
    for (int i = 0; i < 15; i++) {
        const double corr_mag = fabs(correlations[i]);
        const char sign_char = (correlations[i] >= 0.0) ? '+' : '-';
        printf("│  %x  │ %2d  │ %s │ %c%.6e │   2^%7.2f │\n", 
               i + 1, i + 1, masks[i], sign_char, corr_mag, log2(corr_mag));
    }
    printf("└─────┴─────┴────────────────────────────────┴─────────────────────┴─────────────┘\n");
    
    printf("Wildcard analysis completed for all 15 hex values (1-f).\n");
}

static bool is_zero_state(const unsigned char *state) {
    for (int i = 0; i < ORTHROS_STATE_SIZE; ++i) {
        if (state[i] != 0) {
            return false;
        }
    }
    return true;
}

// ================================
// Command Line Interface
// ================================

/**
 * @brief Print usage information
 */
static void print_usage(const char *program_name) {
    printf("Orthros Differential-Linear Cryptanalysis Tool (CUDA)\n");
    printf("Usage: %s [options]\n\n", program_name);
    printf("Options:\n");
    printf("  -r, --rounds NUM           Number of rounds (1-%d) [default: 4]\n", ORTHROS_MAX_ROUNDS);
    printf("  -o, --offset NUM           Round offset (0-%d) [default: 4]\n", ORTHROS_MAX_ROUNDS - 1);
    printf("  -m, --mode NUM             Mode (0=left, 1=right, 2=PRF) [default: 0]\n");
    printf("  -d, --diff-left HEX        Left-branch input difference (32 hex nibbles)\n");
    printf("      --diff-right HEX       Right-branch input difference (32 hex nibbles)\n");
    printf("      --mask HEX             Output mask (32 hex nibbles)\n");
    printf("      --wildcard-mask HEX    Wildcard mask with 'x' for variable positions\n");
    printf("  -b, --blocks NUM           Number of blocks [default: %d]\n", DEFAULT_BLOCKS);
    printf("  -t, --threads NUM          Threads per block [default: %d]\n", DEFAULT_THREADS);
    printf("  -e, --experiments NUM      Number of experiments [default: %d, max: %d]\n",
           DEFAULT_EXPERIMENTS, MAX_EXPERIMENTS);
    printf("  -s, --sample-power NUM     Samples per thread = 2^NUM [default: %d, max: %d]\n",
           SAMPLE_POWER_DEFAULT, MAX_SAMPLE_POWER);
    printf("  -g, --gendlct              Generate full DLCT table (single-bit patterns)\n");
    printf("      --self-test            Run Orthros self-test and exit\n");
    printf("  -v, --verbose              Enable verbose output\n");
    printf("  -h, --help                 Show this help message\n\n");
    printf("Examples:\n");
    printf("  # 5-round PRF Distinguisher 0 (from paper)\n");
    printf("  %s -r 5 -o 3 -m 2 -d 00002020000000000000000000000000\\\n", program_name);
    printf("         --diff-right 00000000000044000000000000000000\\\n");
    printf("         --mask 00000000000000000000000000002022 -e 2 -s 10 -b 256 -t 512 -v\n\n");
    printf("  # Wildcard mask example (variable positions marked with 'x')\n");
    printf("  %s -r 5 -o 3 -m 0 -d 00002020000000000000000000000000\\\n", program_name);
    printf("         --diff-right 00000000000000000000000000000000\\\n");
    printf("         --wildcard-mask 0000000000000000000000000000x0xx -e 10 -s 10 -v\n\n");
    printf("  # Run self-test\n");
    printf("  %s --self-test\n\n", program_name);
    printf("Notes:\n");
    printf("  - Hex values must be exactly 32 nibbles (128 bits)\n");
    printf("  - offset + rounds must be ≤ %d\n", ORTHROS_MAX_ROUNDS);
    printf("  - Mode 0: Left branch, Mode 1: Right branch, Mode 2: PRF\n");
}

/**
 * @brief Parse command line arguments
 */
static Config parse_arguments(int argc, char **argv) {
    Config config = {0};
    
    // Set defaults
    config.blocks = DEFAULT_BLOCKS;
    config.threads = DEFAULT_THREADS;
    config.rounds = 4;
    config.offset = 4;
    config.experiments = DEFAULT_EXPERIMENTS;
    config.sample_power = SAMPLE_POWER_DEFAULT;
    config.mode = ORTHROS_MODE_LEFT;
    config.verbose = false;
    config.self_test_only = false;
    config.generate_dlct = false;
    config.use_wildcard_mask = false;
    
    // Initialize arrays to zero
    memset(config.input_diff_left, 0, ORTHROS_STATE_SIZE);
    memset(config.input_diff_right, 0, ORTHROS_STATE_SIZE);
    memset(config.output_mask, 0, ORTHROS_STATE_SIZE);
    memset(config.wildcard_mask_pattern, 0, sizeof(config.wildcard_mask_pattern));

    static struct option long_options[] = {
        {"rounds",       required_argument, 0, 'r'},
        {"offset",       required_argument, 0, 'o'},
        {"mode",         required_argument, 0, 'm'},
        {"diff-left",    required_argument, 0, 'd'},
        {"diff-right",   required_argument, 0, 2},
        {"mask",         required_argument, 0, 3},
        {"wildcard-mask", required_argument, 0, 4},
        {"blocks",       required_argument, 0, 'b'},
        {"threads",      required_argument, 0, 't'},
        {"experiments",  required_argument, 0, 'e'},
        {"sample-power", required_argument, 0, 's'},
        {"gendlct",      no_argument,       0, 'g'},
        {"self-test",    no_argument,       0, 1},
        {"verbose",      no_argument,       0, 'v'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    int c;
    while ((c = getopt_long(argc, argv, "r:o:m:d:b:t:e:s:gvh", long_options, NULL)) != -1) {
        switch (c) {
            case 'r':
                config.rounds = atoi(optarg);
                if (config.rounds < 1 || config.rounds > ORTHROS_MAX_ROUNDS) {
                    fprintf(stderr, "Error: rounds must be between 1 and %d\n", ORTHROS_MAX_ROUNDS);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                config.offset = atoi(optarg);
                if (config.offset < 0 || config.offset >= ORTHROS_MAX_ROUNDS) {
                    fprintf(stderr, "Error: offset must be between 0 and %d\n", ORTHROS_MAX_ROUNDS - 1);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                config.mode = static_cast<OrthrosMode>(atoi(optarg));
                if (config.mode < ORTHROS_MODE_LEFT || config.mode > ORTHROS_MODE_PRF) {
                    fprintf(stderr, "Error: mode must be 0, 1, or 2\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'd':
                if (!parse_hex_nibbles(optarg, config.input_diff_left)) {
                    fprintf(stderr, "Error: diff-left must be 32 hex nibbles\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'b':
                config.blocks = atoi(optarg);
                if (config.blocks < 1 || config.blocks > 65535) {
                    fprintf(stderr, "Error: blocks must be between 1 and 65535\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 't':
                config.threads = atoi(optarg);
                if (config.threads < 1 || config.threads > 1024) {
                    fprintf(stderr, "Error: threads must be between 1 and 1024\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'e':
                config.experiments = atoi(optarg);
                if (config.experiments < 1 || config.experiments > MAX_EXPERIMENTS) {
                    fprintf(stderr, "Error: experiments must be between 1 and %d\n", MAX_EXPERIMENTS);
                    exit(EXIT_FAILURE);
                }
                break;
            case 's':
                config.sample_power = atoi(optarg);
                if (config.sample_power < 2 || config.sample_power > MAX_SAMPLE_POWER) {
                    fprintf(stderr, "Error: sample power must be between 2 and %d\n", MAX_SAMPLE_POWER);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'g':
                config.generate_dlct = true;
                break;
            case 1: // --self-test
                config.self_test_only = true;
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            case 2: // --diff-right
                if (!parse_hex_nibbles(optarg, config.input_diff_right)) {
                    fprintf(stderr, "Error: diff-right must be 32 hex nibbles\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 3: // --mask
                if (!parse_hex_nibbles(optarg, config.output_mask)) {
                    fprintf(stderr, "Error: mask must be 32 hex nibbles\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 4: // --wildcard-mask
                if (strlen(optarg) != ORTHROS_STATE_SIZE) {
                    fprintf(stderr, "Error: wildcard-mask must be exactly 32 hex characters (with 'x' for wildcards)\n");
                    exit(EXIT_FAILURE);
                }
                strncpy(config.wildcard_mask_pattern, optarg, sizeof(config.wildcard_mask_pattern) - 1);
                config.use_wildcard_mask = true;
                break;
            case '?':
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    // Validation
    if (config.offset + config.rounds > ORTHROS_MAX_ROUNDS) {
        fprintf(stderr, "Error: offset + rounds must be ≤ %d\n", ORTHROS_MAX_ROUNDS);
        exit(EXIT_FAILURE);
    }
    // Skip difference validation for --gendlct and --self-test modes
    if (!config.generate_dlct && !config.self_test_only) {
        if (config.mode == ORTHROS_MODE_LEFT && is_zero_state(config.input_diff_left)) {
            fprintf(stderr, "Error: Provide a non-zero left difference for mode 0\n");
            exit(EXIT_FAILURE);
        }
        if (config.mode == ORTHROS_MODE_RIGHT && is_zero_state(config.input_diff_right)) {
            fprintf(stderr, "Error: Provide a non-zero right difference for mode 1\n");
            exit(EXIT_FAILURE);
        }
    }
    
    return config;
}

// ================================
// Display Functions (Stubs)
// ================================
// Main Function
// ================================

int main(int argc, char **argv) {
    printf("Orthros Differential-Linear Cryptanalysis Tool (CUDA)\n");
    printf("=====================================================================\n\n");
    
    // Initialize memory context
    init_memory_context(&memory_ctx);
    
    // Parse arguments
    const Config config = parse_arguments(argc, argv);
    
    // Report detailed randomness source
    report_random_source();
    
    // Initialize CUDA
    initialize_cuda_device(&config);
    
    // Handle self-test
    if (config.self_test_only) {
        printf("Running Orthros self-test...\n");
        const bool test_result = run_device_self_test();
        printf("Self-test: %s\n", test_result ? "PASSED" : "FAILED");
        cleanup_memory_context(&memory_ctx);
        return test_result ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    
    // Main computation
    if (config.generate_dlct) {
        generate_dlct_table(config);
    } else if (config.use_wildcard_mask) {
        // Run wildcard mask expansion experiments
        run_wildcard_experiments(config);
    } else {
        const double correlation = compute_correlation(config);
        const double corr_mag = fabs(correlation);
        
        char left_hex[ORTHROS_STATE_SIZE + 1];
        char right_hex[ORTHROS_STATE_SIZE + 1];
        char mask_hex[ORTHROS_STATE_SIZE + 1];
        format_nibbles(config.input_diff_left, left_hex);
        format_nibbles(config.input_diff_right, right_hex);
        format_nibbles(config.output_mask, mask_hex);
        
        printf("Results:\n");
        printf("========\n");
        printf("Configuration:\n");
        printf("  Blocks: %d\n", config.blocks);
        printf("  Threads per block: %d\n", config.threads);
        printf("  Total threads: %d × %d = %d\n", config.blocks, config.threads, config.blocks * config.threads);
        printf("  Samples per thread: 2^%d = %d\n", config.sample_power, 1 << config.sample_power);
        printf("  Total samples: %lld ≈ 2^%.1f\n", 
               (long long)config.blocks * config.threads * (1LL << config.sample_power),
               log2((double)config.blocks * config.threads * (1LL << config.sample_power)));
        printf("  Experiments: %d\n\n", config.experiments);
        printf("Analysis:\n");
        printf("Rounds: %d (offset %d)\n", config.rounds, config.offset);
        printf("Mode: %d\n", static_cast<int>(config.mode));
        printf("ΔL: %s\n", left_hex);
        printf("ΔR: %s\n", right_hex);
        printf("Output mask: %s\n", mask_hex);
        
        if (corr_mag > 0.0) {
            const double log_corr = log2(corr_mag);
            const char sign_char = (correlation >= 0.0) ? '+' : '-';
            printf("Average correlation: %c2^%.4f (%+.8e)\n", sign_char, log_corr, correlation);
        } else {
            printf("Average correlation: +2^(-inf) (0.0)\n");
        }
    }
    
    // Report actual random source usage
    report_random_usage_stats();
    
    // Cleanup and exit
    cleanup_memory_context(&memory_ctx);
    return EXIT_SUCCESS;
}