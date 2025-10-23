/*
 * Copyright (C) 2025 Hosein Hadipour
 * Email: hsn.hadipour@gmail.com
 * Date: September 24, 2025
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
 * OpenMP Implementation of Orthros Differential-Linear Cryptanalysis Tool (DLCT)
 *
 * This CPU/OpenMP driver mirrors the CUDA workflow to enable experiments
 * without GPU access. It reuses the shared Orthros core implementation and
 * provides the same command-line interface as the CUDA tool wherever feasible.
 */

#include "common_types.h"
#include "openmp/random_utils.h"
#include "orthros_core.h"

// Forward declarations for DLCT functions implemented in dlct_runner.cpp
double run_correlation_suite(Config config);
void generate_dlct_table_cpu(Config config);

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <getopt.h>
#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif


// ================================
// Configuration constants (shared with CUDA version)
// ================================
static constexpr int DEFAULT_BLOCKS        = 64;   // Reduced for OpenMP CPU performance
static constexpr int DEFAULT_THREADS       = 128;  // Reduced for OpenMP CPU performance
static constexpr int DEFAULT_EXPERIMENTS   = 16;
static constexpr int SAMPLE_POWER_DEFAULT  = 18;
static constexpr int MAX_EXPERIMENTS       = 1000;
static constexpr int MAX_SAMPLE_POWER      = 32;

// ================================
// Hex utilities (duplicated from CUDA driver for CLI parity)
// ================================

static int hex_char_to_nibble(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static bool parse_hex_nibbles(const char *hex, unsigned char *dst) {
    const size_t len = std::strlen(hex);
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
           std::log2((double)base_config.blocks * base_config.threads * (1LL << base_config.sample_power)));
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
        const double correlation = run_correlation_suite(config);
        const double corr_mag = std::fabs(correlation);
        correlations[hex_val - 1] = correlation;
        
        printf("Correlation: %+.12e (2^%.4f)\n", correlation, std::log2(corr_mag));
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
        const double corr_mag = std::fabs(correlations[i]);
        const char sign_char = (correlations[i] >= 0.0) ? '+' : '-';
        printf("│  %x  │ %2d  │ %s │ %c%.6e │   2^%7.2f │\n", 
               i + 1, i + 1, masks[i], sign_char, corr_mag, std::log2(corr_mag));
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

static bool orthros_run_self_test() {
    const unsigned char pt[ORTHROS_STATE_SIZE] = {
        0xa,0x9,0x4,0x7,0x4,0x3,0x6,0x7,0x1,0x0,0x9,0x2,0x4,0xc,0xc,0xd,
        0x4,0x7,0xf,0x2,0xd,0x5,0x7,0x1,0xd,0xe,0xe,0xa,0x8,0xf,0x0,0x5
    };
    const unsigned char key[ORTHROS_STATE_SIZE] = {
        0x4,0xa,0x2,0xb,0xe,0x6,0x0,0xe,0x3,0xd,0xb,0x6,0xa,0xb,0xe,0x0,
        0xc,0x0,0x3,0xe,0xa,0xe,0xc,0x6,0x6,0xf,0xd,0x0,0x5,0xd,0x0,0xc
    };
    const unsigned char expected[ORTHROS_STATE_SIZE] = {
        0xf,0xb,0xf,0xd,0x0,0xb,0x2,0xd,0xc,0x5,0x8,0xc,0x8,0x3,0x5,0x2,
        0x7,0x4,0x4,0xb,0x9,0x7,0x6,0x1,0x9,0x1,0x7,0x6,0xe,0x8,0xd,0x4
    };

    unsigned char ct[ORTHROS_STATE_SIZE];
    orthros_encrypt(pt, pt, ct, key, 0, ORTHROS_MAX_ROUNDS, ORTHROS_MODE_LEFT);
    for (int i = 0; i < ORTHROS_STATE_SIZE; ++i) {
        if (ct[i] != expected[i]) {
            return false;
        }
    }
    return true;
}

// ================================
// CLI helpers
// ================================

static void print_usage(const char *program_name) {
    printf("Orthros Differential-Linear Cryptanalysis Tool (OpenMP)\n");
    printf("Usage: %s [options]\n\n", program_name);
    printf("Options:\n");
    printf("  -r, --rounds NUM           Number of rounds (1-%d) [default: 4]\n", ORTHROS_MAX_ROUNDS);
    printf("  -o, --offset NUM           Round offset (0-%d) [default: 4]\n", ORTHROS_MAX_ROUNDS - 1);
    printf("  -m, --mode NUM             Mode (0=left, 1=right, 2=PRF) [default: 0]\n");
    printf("  -d, --diff-left HEX        Left-branch input difference (32 hex nibbles)\n");
    printf("      --diff-right HEX       Right-branch input difference (32 hex nibbles)\n");
    printf("      --mask HEX             Output mask (32 hex nibbles)\n");
    printf("      --wildcard-mask HEX    Expand wildcards ('x') to all hex combinations\n");
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
    printf("         --mask 00000000000000000000000000002022 -e 2 -s 24 -b 4 -t 4 -v\n\n");
    printf("  # Wildcard mask example (variable positions marked with 'x')\n");
    printf("  %s -r 5 -o 3 -m 0 -d 00002020000000000000000000000000\\\n", program_name);
    printf("         --diff-right 00000000000000000000000000000000\\\n");
    printf("         --wildcard-mask 0000000000000000000000000000x0xx -b 4 -t 8 -e 2 -s 15 -v\n\n");
    printf("  # Run self-test\n");
    printf("  %s --self-test\n\n", program_name);
    printf("Notes:\n");
    printf("  - Hex values must be exactly 32 nibbles (128 bits)\n");
    printf("  - offset + rounds must be ≤ %d\n", ORTHROS_MAX_ROUNDS);
    printf("  - Mode 0: Left branch, Mode 1: Right branch, Mode 2: PRF\n");
}

static void set_default_config(Config *config) {
    config->blocks = DEFAULT_BLOCKS;
    config->threads = DEFAULT_THREADS;
    config->rounds = 4;
    config->offset = 4;
    config->experiments = DEFAULT_EXPERIMENTS;
    config->sample_power = SAMPLE_POWER_DEFAULT;
    config->mode = ORTHROS_MODE_LEFT;
    config->verbose = false;
    config->self_test_only = false;
    config->generate_dlct = false;
    config->use_wildcard_mask = false;
    memset(config->wildcard_mask_pattern, 0, sizeof(config->wildcard_mask_pattern));

    parse_hex_nibbles("00000000000000000000110111100000", config->input_diff_left);
    parse_hex_nibbles("00000000000000008008000000000000", config->input_diff_right);
    parse_hex_nibbles("00000000ff0f00000000000000000000", config->output_mask);
}

static Config parse_arguments(int argc, char **argv) {
    Config config;
    set_default_config(&config);

    static struct option long_options[] = {
        {"rounds",       required_argument, 0, 'r'},
        {"offset",       required_argument, 0, 'o'},
        {"mode",         required_argument, 0, 'm'},
        {"diff-left",    required_argument, 0, 'd'},
        {"diff-right",   required_argument, 0,  2 },
        {"mask",         required_argument, 0,  3 },
        {"wildcard-mask", required_argument, 0, 4 },
        {"blocks",       required_argument, 0, 'b'},
        {"threads",      required_argument, 0, 't'},
        {"experiments",  required_argument, 0, 'e'},
        {"sample-power", required_argument, 0, 's'},
        {"gendlct",      no_argument,       0, 'g'},
        {"self-test",    no_argument,       0,  1 },
        {"verbose",      no_argument,       0, 'v'},
        {"help",         no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "r:o:m:d:b:t:e:s:g vh", long_options, &option_index)) != -1) {
        switch (c) {
            case 'r':
                config.rounds = std::atoi(optarg);
                if (config.rounds < 1 || config.rounds > ORTHROS_MAX_ROUNDS) {
                    fprintf(stderr, "Error: rounds must be between 1 and %d\n", ORTHROS_MAX_ROUNDS);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'o':
                config.offset = std::atoi(optarg);
                if (config.offset < 0 || config.offset >= ORTHROS_MAX_ROUNDS) {
                    fprintf(stderr, "Error: offset must be between 0 and %d\n", ORTHROS_MAX_ROUNDS - 1);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'm':
                config.mode = static_cast<OrthrosMode>(std::atoi(optarg));
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
                config.blocks = std::atoi(optarg);
                if (config.blocks < 1 || config.blocks > 65535) {
                    fprintf(stderr, "Error: blocks must be between 1 and 65535\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 't':
                config.threads = std::atoi(optarg);
                if (config.threads < 1 || config.threads > 1024) {
                    fprintf(stderr, "Error: threads must be between 1 and 1024\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 'e':
                config.experiments = std::atoi(optarg);
                if (config.experiments < 1 || config.experiments > MAX_EXPERIMENTS) {
                    fprintf(stderr, "Error: experiments must be between 1 and %d\n", MAX_EXPERIMENTS);
                    exit(EXIT_FAILURE);
                }
                break;
            case 's':
                config.sample_power = std::atoi(optarg);
                if (config.sample_power < 2 || config.sample_power > MAX_SAMPLE_POWER) {
                    fprintf(stderr, "Error: sample power must be between 2 and %d\n", MAX_SAMPLE_POWER);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'g':
                config.generate_dlct = true;
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_SUCCESS);
            case 1:
                config.self_test_only = true;
                break;
            case 2:
                if (!parse_hex_nibbles(optarg, config.input_diff_right)) {
                    fprintf(stderr, "Error: diff-right must be 32 hex nibbles\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 3:
                if (!parse_hex_nibbles(optarg, config.output_mask)) {
                    fprintf(stderr, "Error: mask must be 32 hex nibbles\n");
                    exit(EXIT_FAILURE);
                }
                break;
            case 4:
                // Wildcard mask pattern
                if (strlen(optarg) != ORTHROS_STATE_SIZE) {
                    fprintf(stderr, "Error: wildcard-mask must be exactly 32 characters\n");
                    exit(EXIT_FAILURE);
                }
                for (size_t i = 0; i < strlen(optarg); i++) {
                    char c = optarg[i];
                    if (!(c == 'x' || c == 'X' || (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F'))) {
                        fprintf(stderr, "Error: wildcard-mask must contain only hex characters (0-9, a-f, A-F) and wildcards (x, X)\n");
                        exit(EXIT_FAILURE);
                    }
                }
                strcpy(config.wildcard_mask_pattern, optarg);
                config.use_wildcard_mask = true;
                break;
            case '?':
            default:
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
        }
    }

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
// Entry point
// ================================

int main(int argc, char **argv) {
    printf("Orthros Differential-Linear Cryptanalysis Tool (OpenMP)\n");
    printf("=====================================================================\n\n");

#ifdef _OPENMP
    printf("OpenMP detected: using up to %d host threads\n\n", omp_get_max_threads());
#else
    printf("Warning: OpenMP not enabled; falling back to serial execution\n\n");
#endif

    Config config = parse_arguments(argc, argv);
    reset_random_usage();
    report_random_source();

    if (config.self_test_only) {
        printf("Running Orthros self-test...\n");
        const bool test_result = orthros_run_self_test();
        printf("Self-test: %s\n", test_result ? "PASSED" : "FAILED");
        report_random_usage_stats();
        return test_result ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    if (config.generate_dlct) {
        generate_dlct_table_cpu(config);
        report_random_usage_stats();
        return EXIT_SUCCESS;
    }

    if (config.use_wildcard_mask) {
        // Run wildcard mask expansion experiments
        run_wildcard_experiments(config);
        report_random_usage_stats();
        return EXIT_SUCCESS;
    }

    const double correlation = run_correlation_suite(config);
    const double corr_mag = std::fabs(correlation);

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
           std::log2((double)config.blocks * config.threads * (1LL << config.sample_power)));
    printf("  Experiments: %d\n\n", config.experiments);
    printf("Analysis:\n");
    printf("Rounds: %d (offset %d)\n", config.rounds, config.offset);
    printf("Mode: %d\n", static_cast<int>(config.mode));
    printf("ΔL: %s\n", left_hex);
    printf("ΔR: %s\n", right_hex);
    printf("Output mask: %s\n", mask_hex);

    if (corr_mag > 0.0) {
        const double log_corr = std::log2(corr_mag);
        const char sign_char = (correlation >= 0.0) ? '+' : '-';
        printf("Average correlation: %c2^%.4f (%+.8e)\n", sign_char, log_corr, correlation);
    } else {
        printf("Average correlation: +2^(-inf) (0.0)\n");
    }

    report_random_usage_stats();
    return EXIT_SUCCESS;
}
