/*
 * Key Analysis Tool for Orthros-PRF Differential-Linear Attacks
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
 * In case you use this tool please include the above copyright
 * information (name, contact, license).
 *
 * Description:
 * This program analyzes the valid key space for differential-linear attacks on Orthros-PRF.
 * It identifies weak-key classes and computes statistical parameters for key recovery attacks.
 * 
 * Note: This program is applicable to key recovery attacks with only 2 active S-boxes 
 * in the key recovery part.
 *
 * Key Concepts:
 * - Good Pair: A pair (x, y) that satisfies the differential property for given key nibbles
 * - Weak-Key Class: A set of key candidates sharing the same sorted list of good pairs
 * - Ordered vs Unordered Pairs:
 *   * Ordered: Both (x,y) and (y,x) are counted as separate pairs
 *   * Unordered: Only pairs where x < y are counted (matches Yosuke's convention)
 *
 * Compile (C++ with Boost.Math):
 *     g++ -O3 -march=native -std=c++17 -o keyanalysis keyanalysis.cpp
 *     make keyanalysis
 *     
 * Usage:
 *     ./keyanalysis --version <version_name> [--output <filename>]
 *     ./keyanalysis --custom --offset <n> --active <i1,i2> --dy1 <hex> --dy2 <hex> [--output <filename>]
 *     ./keyanalysis --version <version_name> --log2N <n> --Psucc <p> --capexp <e> [--unordered]
 * 
 * Examples:
 *     ./keyanalysis --version 6r-v2 --log2N 20 --Psucc 0.7 --capexp 18.09
 *     ./keyanalysis --version 6r-v2 --unordered --log2N 20 --Psucc 0.7 --capexp 18.09
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <boost/math/distributions/normal.hpp>

#define N 8
#define MAX_GOOD_PAIRS 1024

typedef struct {
    const char *name;
    int offset;
    int active_indices[2];
    uint8_t dy1;
    uint8_t dy2;
    const char *description;
} AttackVersion;

static const AttackVersion predefined_versions[] = {
    {"5r-v1", 2, {13, 28}, 0x42, 0x44, "5-round attack version 1"},
    {"5r-v2", 2, {6, 13}, 0x88, 0x82, "5-round attack version 2"},
    {"6r-v2", 2, {15, 18}, 0x82, 0x44, "6-round attack version 2"},
    {"8r-v1", 0, {27, 29}, 0x44, 0x81, "8-round attack version 1"},
    {"8r-v5", 0, {20, 21}, 0x12, 0x41, "8-round attack version 5"},
    {"8r-v6", 0, {6, 24}, 0x42, 0x41, "8-round attack version 6"},
    {"8r-v8", 0, {15, 18}, 0x82, 0x44, "8-round attack version 8"},
};

#define NUM_PREDEFINED_VERSIONS (sizeof(predefined_versions) / sizeof(AttackVersion))

static const uint8_t inv_kp_left[128] = {
    0, 11, 58, 29, 47, 36, 100, 49, 111, 33, 83, 73, 76, 68, 118, 65,
    20, 57, 63, 35, 80, 89, 4, 106, 12, 116, 115, 97, 42, 70, 75, 24,
    119, 26, 95, 81, 9, 39, 40, 113, 102, 105, 50, 67, 125, 54, 60, 28,
    7, 107, 74, 41, 92, 1, 16, 109, 85, 82, 37, 69, 127, 8, 31, 62,
    79, 30, 25, 13, 86, 18, 48, 123, 103, 3, 99, 55, 45, 72, 124, 108,
    53, 93, 44, 27, 90, 101, 122, 2, 114, 38, 77, 87, 51, 15, 22, 5,
    126, 59, 91, 6, 121, 34, 110, 17, 52, 96, 43, 61, 10, 64, 117, 21,
    19, 104, 56, 120, 88, 112, 78, 14, 98, 46, 23, 32, 84, 71, 94, 66
};

static const uint8_t inv_kp_right[128] = {
    57, 39, 6, 90, 51, 95, 101, 86, 63, 97, 75, 8, 33, 127, 14, 47,
    76, 37, 25, 124, 83, 71, 103, 56, 94, 123, 61, 98, 89, 20, 1, 4,
    60, 121, 105, 3, 18, 68, 66, 12, 116, 74, 48, 108, 53, 80, 5, 34,
    19, 24, 81, 36, 118, 2, 70, 91, 107, 22, 85, 113, 40, 28, 72, 122,
    111, 79, 38, 87, 119, 64, 32, 43, 82, 125, 26, 54, 0, 96, 112, 7,
    21, 44, 29, 109, 58, 27, 69, 11, 62, 55, 106, 13, 92, 115, 67, 77,
    41, 16, 117, 104, 50, 15, 99, 93, 46, 30, 102, 59, 88, 84, 10, 35,
    120, 52, 45, 23, 42, 78, 17, 114, 110, 65, 100, 73, 49, 9, 31, 126
};

static const uint8_t sb4[16] = {
    0x1, 0x0, 0x2, 0x4, 0x3, 0x8, 0x6, 0xd,
    0x9, 0xa, 0xb, 0xe, 0xf, 0xc, 0x7, 0x5
};

static uint8_t sb8[256];

/*
 * Data Structures:
 * 
 * GoodPair: Represents a pair (x, y) that satisfies the differential property
 *   - x, y: Byte values where S(x⊕k) ⊕ S(y⊕k) = dy for given key nibbles k
 * 
 * GoodPairList: Collection of good pairs for a specific key candidate
 * 
 * KeyPair: A candidate key (k1, k2) for the two active S-boxes
 * 
 * WeakKeyClass: Groups keys that share the same sorted list of good pairs
 *   - pairs[]: Sorted list of good pairs (class identifier)
 *   - keys[]: All key candidates in this class
 *   - Used for efficient lookup during key recovery attack
 */

// Store good pairs as (x, y) byte values for clarity and efficiency
typedef struct {
    uint8_t x;      // First value in the pair
    uint8_t y;      // Second value in the pair (y = x ^ dx for some dx)
} GoodPair;

typedef struct {
    GoodPair pairs[MAX_GOOD_PAIRS];
    int count;
} GoodPairList;

typedef struct {
    uint8_t k1;
    uint8_t k2;
} KeyPair;

typedef struct {
    GoodPair pairs[MAX_GOOD_PAIRS];
    int num_pairs;
    KeyPair *keys;
    int num_keys;
    int max_keys;
} WeakKeyClass;

typedef struct {
    int left;
    int right;
} CommonIndex;

void init_sb8(void) {
    for (int a = 0; a < 16; a++) {
        for (int b = 0; b < 16; b++) {
            sb8[a * 16 + b] = (sb4[a] << 4) | sb4[b];
        }
    }
}

void derive_involved_key_positions(int offset, const int *active_indices, int num_active,
                                   int *left_ik, int *right_ik, int *num_left, int *num_right,
                                   int *common_ik, int *num_common, CommonIndex *common_indices) {
    int active_bits[256];
    int num_active_bits = 0;
    
    for (int i = 0; i < num_active; i++) {
        for (int j = 0; j < 4; j++) {
            active_bits[num_active_bits++] = 4 * active_indices[i] + j;
        }
    }
    
    int inv_kp_left_at_roundr[128];
    int inv_kp_right_at_roundr[128];
    
    for (int i = 0; i < 128; i++) {
        inv_kp_left_at_roundr[i] = inv_kp_left[i];
        inv_kp_right_at_roundr[i] = inv_kp_right[i];
    }
    
    for (int r = 0; r < offset; r++) {
        int temp_left[128], temp_right[128];
        for (int i = 0; i < 128; i++) {
            temp_left[i] = inv_kp_left[inv_kp_left_at_roundr[i]];
            temp_right[i] = inv_kp_right[inv_kp_right_at_roundr[i]];
        }
        memcpy(inv_kp_left_at_roundr, temp_left, sizeof(temp_left));
        memcpy(inv_kp_right_at_roundr, temp_right, sizeof(temp_right));
    }
    
    *num_left = 0;
    *num_right = 0;
    for (int i = 0; i < num_active_bits; i++) {
        left_ik[*num_left] = inv_kp_left_at_roundr[active_bits[i]];
        right_ik[*num_right] = inv_kp_right_at_roundr[active_bits[i]];
        (*num_left)++;
        (*num_right)++;
    }
    
    *num_common = 0;
    int temp_common[256];
    for (int i = 0; i < *num_left; i++) {
        for (int j = 0; j < *num_right; j++) {
            if (left_ik[i] == right_ik[j]) {
                int already_added = 0;
                for (int k = 0; k < *num_common; k++) {
                    if (temp_common[k] == left_ik[i]) {
                        already_added = 1;
                        break;
                    }
                }
                if (!already_added) {
                    temp_common[*num_common] = left_ik[i];
                    common_indices[*num_common].left = i;
                    common_indices[*num_common].right = j;
                    (*num_common)++;
                }
                break;
            }
        }
    }
    
    for (int i = 0; i < *num_common; i++) {
        common_ik[i] = temp_common[i];
    }
}

void check_candidate(uint8_t k1, uint8_t k2, uint8_t dx, uint8_t dy1, uint8_t dy2,
                    GoodPairList *good_pairs, int use_unordered) {
    for (int x = 0; x < 256; x++) {
        uint8_t y = x ^ dx;
        
        // For unordered pairs, only check x < y to avoid counting both (x, y) and (y, x)
        if (use_unordered && x >= y) continue;
        
        // Check if this (x, y) pair satisfies the differential property
        if ((sb8[x ^ k1] ^ sb8[y ^ k1]) == dy1 &&
            (sb8[x ^ k2] ^ sb8[y ^ k2]) == dy2) {
            
            // Check if pair already exists (shouldn't happen, but be safe)
            int already_exists = 0;
            for (int i = 0; i < good_pairs->count; i++) {
                if (good_pairs->pairs[i].x == x && good_pairs->pairs[i].y == y) {
                    already_exists = 1;
                    break;
                }
            }
            
            if (!already_exists && good_pairs->count < MAX_GOOD_PAIRS) {
                good_pairs->pairs[good_pairs->count].x = x;
                good_pairs->pairs[good_pairs->count].y = y;
                good_pairs->count++;
            }
        }
    }
}

// Compare two good pairs for sorting and equality
int compare_good_pairs(const void *a, const void *b) {
    const GoodPair *pa = (const GoodPair *)a;
    const GoodPair *pb = (const GoodPair *)b;
    
    // First compare by x value
    if (pa->x != pb->x) {
        return pa->x - pb->x;
    }
    // If x values are equal, compare by y value
    return pa->y - pb->y;
}

// Compare two good pair lists for equality
int compare_good_pair_lists(const GoodPair *list1, int count1, const GoodPair *list2, int count2) {
    if (count1 != count2) return 0;
    
    for (int i = 0; i < count1; i++) {
        if (list1[i].x != list2[i].x || list1[i].y != list2[i].y) {
            return 0;
        }
    }
    return 1;
}

void print_help(const char *program_name) {
    printf("Usage: %s [OPTIONS]\n\n", program_name);
    printf("Analyze the valid key space for differential-linear attacks on Orthros-PRF.\n\n");
    printf("Options:\n");
    printf("  --version <name>       Use a predefined attack version\n");
    printf("  --list                 List all available predefined versions\n");
    printf("  --custom               Use custom attack parameters (requires --offset, --active, --dy1, --dy2)\n");
    printf("  --offset <n>           Offset value (for custom attacks)\n");
    printf("  --active <i1,i2>       Active S-box indices (comma-separated, for custom attacks)\n");
    printf("  --dy1 <hex>            Output difference in left branch (hex, for custom attacks)\n");
    printf("  --dy2 <hex>            Output difference in right branch (hex, for custom attacks)\n");
    printf("  --output <file>        Save detailed results to file (default: keyanalysis_results.txt)\n");
    printf("  --unordered            Use unordered pairs (default: ordered pairs)\n");
    printf("  --log2N <n>            Log2 of number of plaintext-ciphertext pairs (e.g., 20 for 2^20)\n");
    printf("  --Psucc <p>            Target detection probability (e.g., 0.7)\n");
    printf("  --capexp <e>           Capacity exponent (e.g., 18.09 for capacity = 2^(-18.09))\n");
    printf("  --help                 Display this help message\n\n");
    printf("Pair Counting Mode:\n");
    printf("  By default, the tool counts ordered pairs (x, y) where both (x,y) and (y,x) are counted.\n");
    printf("  Use --unordered to count only unordered pairs where x < y (matches Yosuke's convention).\n\n");
    printf("Statistical Analysis:\n");
    printf("  If --log2N, --Psucc, and --capexp are all provided, the program will:\n");
    printf("  - Compute capacity = 2^(-capexp)\n");
    printf("  - Compute sigma0 = sigma1 = sqrt(capacity) (for close distributions)\n");
    printf("  - Compute the expected number of false positives E\n\n");
    printf("Examples:\n");
    printf("  %s --version 5r-v1\n", program_name);
    printf("  %s --version 8r-v5 --output results.txt\n", program_name);
    printf("  %s --custom --offset 2 --active 13,28 --dy1 0x42 --dy2 0x44\n", program_name);
    printf("  %s --version 6r-v2 --log2N 20 --Psucc 0.7 --capexp 18.09\n", program_name);
    printf("  %s --version 6r-v2 --unordered --log2N 20 --Psucc 0.7 --capexp 18.09\n", program_name);
}

void list_versions(void) {
    printf("Available predefined attack versions:\n\n");
    for (size_t i = 0; i < NUM_PREDEFINED_VERSIONS; i++) {
        printf("  %-8s - %s\n", predefined_versions[i].name, predefined_versions[i].description);
        printf("           offset=%d, active=[%d, %d], dy1=0x%02x, dy2=0x%02x\n\n",
               predefined_versions[i].offset,
               predefined_versions[i].active_indices[0],
               predefined_versions[i].active_indices[1],
               predefined_versions[i].dy1,
               predefined_versions[i].dy2);
    }
}

const AttackVersion* find_version(const char *name) {
    for (size_t i = 0; i < NUM_PREDEFINED_VERSIONS; i++) {
        if (strcmp(predefined_versions[i].name, name) == 0) {
            return &predefined_versions[i];
        }
    }
    return NULL;
}

// Normal CDF and inverse CDF via Boost.Math (high-accuracy)
double norm_cdf(double x) {
    static const boost::math::normal dist(0.0, 1.0);
    return boost::math::cdf(dist, x);
}

double norm_ppf(double p) {
    static const boost::math::normal dist(0.0, 1.0);
    return boost::math::quantile(dist, p);
}

// Compute expected number of false positives E
double compute_E(int num_classes, int *class_sizes, double num_pairs, double P_succ,
                 double capacity, double sigma0, double sigma1) {
    // Symmetric approximation for close distributions (see paper Remark after Theorem 3.4)
    double mu0 = 0.5 * capacity;   // Under H1 (weak key)
    double mu1 = -0.5 * capacity;  // Under H0 (strong key)
    double E = 0.0;
    double Phi_inv_P_succ = norm_ppf(P_succ);
    
    for (int i = 0; i < num_classes; i++) {
        if (class_sizes[i] <= 0) continue;
        
        double M_i = class_sizes[i] * num_pairs;
        double sqrt_M_i = sqrt(M_i);
        
        // Threshold: tau_i = M_i * mu0 - sqrt(M_i) * sigma0 * Phi^{-1}(P_succ)
        double tau_i = M_i * mu0 - sqrt_M_i * sigma0 * Phi_inv_P_succ;
        
        // Type I error: P_FP,i = 1 - Phi((tau_i - M_i * mu1) / (sqrt(M_i) * sigma1))
        double z = (tau_i - M_i * mu1) / (sqrt_M_i * sigma1);
        double P_FP_i = 1.0 - norm_cdf(z);
        
        E += P_FP_i;
    }
    
    return E;
}

// Debug helper: print E breakdown aggregated by class size |X_i|
static void print_E_breakdown(FILE *out,
                              int num_classes, const int *class_sizes,
                              double N_val, double P_succ, double capacity) {
    if (!out) return;
    // Build histogram of class sizes
    int uniq_sizes[1024];
    int uniq_counts[1024];
    int uniq = 0;
    for (int i = 0; i < num_classes; i++) {
        int s = class_sizes[i];
        if (s <= 0) continue;
        int found = -1;
        for (int j = 0; j < uniq; j++) if (uniq_sizes[j] == s) { found = j; break; }
        if (found == -1) {
            if (uniq < (int)(sizeof(uniq_sizes)/sizeof(uniq_sizes[0]))) {
                uniq_sizes[uniq] = s;
                uniq_counts[uniq] = 1;
                uniq++;
            }
        } else {
            uniq_counts[found]++;
        }
    }

    // Compute contributions per size
    double gamma = norm_ppf(P_succ);
    double sumE = 0.0;
    fprintf(out, "E breakdown by class size |X_i| (debug):\n");
    fprintf(out, "  Using z = sqrt(N*|X_i|*capacity) - Phi^{-1}(P_succ)\n");
    fprintf(out, "  N = %.0f, capacity = %.10f, Phi^{-1}(P_succ) = %.6f\n",
            N_val, capacity, gamma);
    for (int j = 0; j < uniq; j++) {
        int s = uniq_sizes[j];
        int c = uniq_counts[j];
        double M = N_val * (double)s;
        double z = sqrt(M * capacity) - gamma;
        double pfp = 1.0 - norm_cdf(z);
        double contrib = c * pfp;
        sumE += contrib;
        fprintf(out, "  |X_i| = %d, count = %d, z = %.6f, PFP = %.6e, contrib = %.6f\n",
                s, c, z, pfp, contrib);
    }
    fprintf(out, "  Sum contributions = %.6f\n", sumE);
}

int main(int argc, char *argv[]) {
    init_sb8();
    
    // Default values
    int offset = 2;
    int active_indices[] = {13, 28};
    int num_active = 2;
    uint8_t dy1 = 0x42;
    uint8_t dy2 = 0x44;
    const char *output_filename = "keyanalysis_results.txt";
    const char *version_name = NULL;
    int use_custom = 0;
    int custom_offset = -1;
    int custom_active[2] = {-1, -1};
    int custom_dy1 = -1;
    int custom_dy2 = -1;
    int use_unordered = 0;     // Flag for unordered pairs (default: ordered)
    
    // Statistical parameters for E computation (optional)
    int log2_N = -1;           // log2 of number of plaintext-ciphertext pairs
    double P_succ = -1.0;      // Target detection probability
    double capacity_exp = -1.0; // Exponent for capacity (e.g., 18.09 for 2^(-18.09))
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--help") == 0) {
            print_help(argv[0]);
            return 0;
        } else if (strcmp(argv[i], "--list") == 0) {
            list_versions();
            return 0;
        } else if (strcmp(argv[i], "--version") == 0) {
            if (i + 1 < argc) {
                version_name = argv[++i];
            } else {
                fprintf(stderr, "Error: --version requires a version name\n");
                print_help(argv[0]);
                return 1;
            }
        } else if (strcmp(argv[i], "--output") == 0) {
            if (i + 1 < argc) {
                output_filename = argv[++i];
            } else {
                fprintf(stderr, "Error: --output requires a filename\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--custom") == 0) {
            use_custom = 1;
        } else if (strcmp(argv[i], "--offset") == 0) {
            if (i + 1 < argc) {
                custom_offset = atoi(argv[++i]);
            } else {
                fprintf(stderr, "Error: --offset requires a value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--active") == 0) {
            if (i + 1 < argc) {
                char *token = strtok(argv[++i], ",");
                if (token) custom_active[0] = atoi(token);
                token = strtok(NULL, ",");
                if (token) custom_active[1] = atoi(token);
            } else {
                fprintf(stderr, "Error: --active requires comma-separated indices\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--dy1") == 0) {
            if (i + 1 < argc) {
                custom_dy1 = (int)strtol(argv[++i], NULL, 16);
            } else {
                fprintf(stderr, "Error: --dy1 requires a hex value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--dy2") == 0) {
            if (i + 1 < argc) {
                custom_dy2 = (int)strtol(argv[++i], NULL, 16);
            } else {
                fprintf(stderr, "Error: --dy2 requires a hex value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--log2N") == 0) {
            if (i + 1 < argc) {
                log2_N = atoi(argv[++i]);
            } else {
                fprintf(stderr, "Error: --log2N requires an integer value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--Psucc") == 0) {
            if (i + 1 < argc) {
                P_succ = atof(argv[++i]);
            } else {
                fprintf(stderr, "Error: --Psucc requires a value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--capexp") == 0) {
            if (i + 1 < argc) {
                capacity_exp = atof(argv[++i]);
            } else {
                fprintf(stderr, "Error: --capexp requires a value\n");
                return 1;
            }
        } else if (strcmp(argv[i], "--unordered") == 0) {
            use_unordered = 1;
        }
    }
    
    // Apply version or custom parameters
    if (version_name) {
        const AttackVersion *version = find_version(version_name);
        if (version == NULL) {
            fprintf(stderr, "Error: Unknown version '%s'\n", version_name);
            fprintf(stderr, "Use --list to see available versions\n");
            return 1;
        }
        offset = version->offset;
        active_indices[0] = version->active_indices[0];
        active_indices[1] = version->active_indices[1];
        dy1 = version->dy1;
        dy2 = version->dy2;
        printf("Using attack version: %s - %s\n", version->name, version->description);
    } else if (use_custom) {
        if (custom_offset < 0 || custom_active[0] < 0 || custom_active[1] < 0 ||
            custom_dy1 < 0 || custom_dy2 < 0) {
            fprintf(stderr, "Error: Custom attack requires --offset, --active, --dy1, and --dy2\n");
            print_help(argv[0]);
            return 1;
        }
        offset = custom_offset;
        active_indices[0] = custom_active[0];
        active_indices[1] = custom_active[1];
        dy1 = (uint8_t)custom_dy1;
        dy2 = (uint8_t)custom_dy2;
        printf("Using custom attack parameters\n");
    } else {
        // Default to 5r-v1
        printf("Using default attack version: 5r-v1 (use --help for options)\n");
    }
    
    printf("Parameters: offset=%d, active=[%d, %d], dy1=0x%02x, dy2=0x%02x\n",
           offset, active_indices[0], active_indices[1], dy1, dy2);
    printf("Pair counting mode: %s\n", use_unordered ? "unordered (x < y only)" : "ordered (all pairs)");
    printf("Writing detailed results to: %s\n\n", output_filename);
    
    // Open output file
    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
        fprintf(stderr, "Error: Cannot open output file '%s'\n", output_filename);
        return 1;
    }
    
    fprintf(output_file, "Orthros-PRF Key Analysis Results\n");
    fprintf(output_file, "=================================\n");
    if (version_name) {
        fprintf(output_file, "Attack version: %s\n", version_name);
    } else {
        fprintf(output_file, "Attack version: custom\n");
    }
    fprintf(output_file, "Parameters: offset=%d, active=[%d, %d], dy1=0x%02x, dy2=0x%02x\n",
            offset, active_indices[0], active_indices[1], dy1, dy2);
    fprintf(output_file, "Pair counting mode: %s\n\n", use_unordered ? "unordered (x < y only)" : "ordered (all pairs)");
    
    int left_ik[256], right_ik[256], common_ik[256];
    int num_left, num_right, num_common;
    CommonIndex common_indices[256];
    
    derive_involved_key_positions(offset, active_indices, num_active,
                                 left_ik, right_ik, &num_left, &num_right,
                                 common_ik, &num_common, common_indices);
    
    uint64_t total_number_of_keys = 1ULL << (num_left + num_right - num_common);
    
    clock_t start_time = clock();
    
    int cnt = 0;
    WeakKeyClass *key_classes = (WeakKeyClass *)malloc(65536 * sizeof(WeakKeyClass));
    if (key_classes == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }
    int num_classes = 0;
    
    for (int k1 = 0; k1 < 256; k1++) {
        for (int k2 = 0; k2 < 256; k2++) {
            int flag = 1;
            
            for (int i = 0; i < num_common; i++) {
                int ls = 8 - common_indices[i].left - 1;
                int rs = 8 - common_indices[i].right - 1;
                if (((k1 >> ls) & 0x1) != ((k2 >> rs) & 0x1)) {
                    flag = 0;
                    break;
                }
            }
            
            if (flag) {
                GoodPairList good_pairs; good_pairs.count = 0;
                
                for (int dx = 1; dx < 256; dx++) {
                    check_candidate(k1, k2, dx, dy1, dy2, &good_pairs, use_unordered);
                }
                
                if (good_pairs.count > 0) {
                    // Write detailed output to file using the raw branch subkeys
                    fprintf(output_file,
                            "cnt = %5d, k_left = 0x%02x, k_right = 0x%02x, good pairs: [",
                            cnt, k1, k2);
                    for (int i = 0; i < good_pairs.count; i++) {
                        fprintf(output_file, "(0x%02x, 0x%02x)", 
                               good_pairs.pairs[i].x,
                               good_pairs.pairs[i].y);
                        if (i < good_pairs.count - 1) fprintf(output_file, ", ");
                    }
                    fprintf(output_file, "]\n");
                    
                    qsort(good_pairs.pairs, good_pairs.count, sizeof(GoodPair), compare_good_pairs);
                    
                    int class_idx = -1;
                    for (int i = 0; i < num_classes; i++) {
                        if (compare_good_pair_lists(key_classes[i].pairs, key_classes[i].num_pairs,
                                                    good_pairs.pairs, good_pairs.count)) {
                            class_idx = i;
                            break;
                        }
                    }
                    
                    if (class_idx == -1) {
                        class_idx = num_classes;
                        memcpy(key_classes[class_idx].pairs, good_pairs.pairs,
                               good_pairs.count * sizeof(GoodPair));
                        key_classes[class_idx].num_pairs = good_pairs.count;
                        key_classes[class_idx].num_keys = 0;
                        key_classes[class_idx].max_keys = 16;
                        key_classes[class_idx].keys = (KeyPair *)malloc(16 * sizeof(KeyPair));
                        num_classes++;
                    }
                    
                    if (key_classes[class_idx].num_keys >= key_classes[class_idx].max_keys) {
                        key_classes[class_idx].max_keys *= 2;
                        key_classes[class_idx].keys = (KeyPair *)realloc(key_classes[class_idx].keys,
                                                              key_classes[class_idx].max_keys * sizeof(KeyPair));
                    }
                    
                    key_classes[class_idx].keys[key_classes[class_idx].num_keys].k1 = k1;
                    key_classes[class_idx].keys[key_classes[class_idx].num_keys].k2 = k2;
                    key_classes[class_idx].num_keys++;
                    
                    cnt++;
                }
            }
        }
    }
    
    fprintf(output_file, "\n");
    fprintf(output_file, "Weak Key Classes:\n");
    fprintf(output_file, "=================\n");
    for (int i = 0; i < num_classes; i++) {
        fprintf(output_file, "Good pairs: (");
        for (int j = 0; j < key_classes[i].num_pairs; j++) {
            fprintf(output_file, "(0x%02x, 0x%02x)", 
                   key_classes[i].pairs[j].x,
                   key_classes[i].pairs[j].y);
            if (j < key_classes[i].num_pairs - 1) fprintf(output_file, ", ");
        }
        fprintf(output_file, "), Key candidates: {");
        for (int k = 0; k < key_classes[i].num_keys; k++) {
            fprintf(output_file, "(left=0x%02x, right=0x%02x)",
                   key_classes[i].keys[k].k1,
                   key_classes[i].keys[k].k2);
            if (k < key_classes[i].num_keys - 1) fprintf(output_file, ", ");
        }
        fprintf(output_file, "}\n");
    }
    
    double elapsed_time = (double)(clock() - start_time) / CLOCKS_PER_SEC;
    
    // Calculate union of all good pairs
    GoodPair all_good_pairs[MAX_GOOD_PAIRS * 100];  // Large enough buffer
    int total_unique_pairs = 0;
    
    for (int i = 0; i < num_classes; i++) {
        for (int j = 0; j < key_classes[i].num_pairs; j++) {
            // Check if this pair already exists in our union
            int already_exists = 0;
            for (int k = 0; k < total_unique_pairs; k++) {
                if (compare_good_pairs(&key_classes[i].pairs[j], &all_good_pairs[k]) == 0) {
                    already_exists = 1;
                    break;
                }
            }
            if (!already_exists && total_unique_pairs < MAX_GOOD_PAIRS * 100) {
                all_good_pairs[total_unique_pairs] = key_classes[i].pairs[j];
                total_unique_pairs++;
            }
        }
    }
    
    // Write summary to both file and terminal
    const char *separator = "#######################################################################################\n";
    fprintf(output_file, "\n%s", separator);
    fprintf(output_file, "Attack Summary:\n");
    fprintf(output_file, "===============\n");
    fprintf(output_file, "Elapsed time: %.2f seconds\n", elapsed_time);
    fprintf(output_file, "left involved key bits: [");
    for (int i = 0; i < num_left; i++) {
        fprintf(output_file, "%d", left_ik[i]);
        if (i < num_left - 1) fprintf(output_file, ", ");
    }
    fprintf(output_file, "]\n");
    fprintf(output_file, "right involved key bits: [");
    for (int i = 0; i < num_right; i++) {
        fprintf(output_file, "%d", right_ik[i]);
        if (i < num_right - 1) fprintf(output_file, ", ");
    }
    fprintf(output_file, "]\n");
    fprintf(output_file, "Common key bits: {");
    for (int i = 0; i < num_common; i++) {
        fprintf(output_file, "%d", common_ik[i]);
        if (i < num_common - 1) fprintf(output_file, ", ");
    }
    fprintf(output_file, "}\n");
    fprintf(output_file, "Total number of keys: %llu\n", (unsigned long long)total_number_of_keys);
    fprintf(output_file, "Number of weak keys: %d\n", cnt);
    fprintf(output_file, "Number of unique indices for the hash table (indices are sorted good pairs): %d\n", num_classes);
    fprintf(output_file, "Size of union of all good pairs: %d\n", total_unique_pairs);
    fprintf(output_file, "The conditions holds for %d keys out of %llu keys.\n", cnt, (unsigned long long)total_number_of_keys);
    
    if (cnt > 0) {
        double log2_weak_key_space = log2((double)cnt / total_number_of_keys);
        fprintf(output_file, "The attack works for 2^(%.2f) portion of the key space.\n", log2_weak_key_space);
    }
    
    if ((uint64_t)cnt < total_number_of_keys) {
        double learned_bits = log2((double)(total_number_of_keys - cnt) / total_number_of_keys);
        fprintf(output_file, "Number of learned bits in the worst case (when the key is not weak): %.2f\n", fabs(learned_bits));
    }
    
    // Compute expected number of false positives E if statistical parameters are provided
    if (log2_N >= 0 && P_succ > 0 && capacity_exp > 0) {
        // Calculate capacity and standard deviations from exponent
        double capacity = pow(2.0, -capacity_exp);
        // For close distributions: sigma_0^2 ≈ sigma_1^2 ≈ capacity
        double sigma0 = sqrt(capacity);
        double sigma1 = sqrt(capacity);
        
        // Collect class sizes (|X_i| = number of good pairs for each class)
        int *class_sizes = (int *)malloc(num_classes * sizeof(int));
        int sum_Xi = 0;
        for (int i = 0; i < num_classes; i++) {
            class_sizes[i] = key_classes[i].num_pairs;
            sum_Xi += class_sizes[i];
        }
        
        // Compute class size statistics
        int min_size = class_sizes[0];
        int max_size = class_sizes[0];
        double avg_size = (double)sum_Xi / num_classes;
        for (int i = 0; i < num_classes; i++) {
            if (class_sizes[i] < min_size) min_size = class_sizes[i];
            if (class_sizes[i] > max_size) max_size = class_sizes[i];
        }
        
        double N_val = pow(2.0, log2_N);
        double E = compute_E(num_classes, class_sizes, N_val, P_succ, capacity, sigma0, sigma1);
        
        fprintf(output_file, "\nStatistical Analysis:\n");
        fprintf(output_file, "=====================\n");
        fprintf(output_file, "Number of weak-key classes: %d\n", num_classes);
        fprintf(output_file, "Sum of |X_i| across all classes: %d (log2 ≈ %.2f)\n", sum_Xi, log2((double)sum_Xi));
        fprintf(output_file, "Class size distribution: min=%d, max=%d, avg=%.2f\n", min_size, max_size, avg_size);
        fprintf(output_file, "N = 2^%d = %.0f (number of plaintext-ciphertext pairs)\n", log2_N, N_val);
        fprintf(output_file, "P_succ = %.2f (target detection probability)\n", P_succ);
        fprintf(output_file, "capacity = 2^(-%.2f) = %.10f (distinguisher capacity)\n", capacity_exp, capacity);
        fprintf(output_file, "sigma0 = sigma1 = sqrt(capacity) = %.6f (standard deviations)\n", sigma0);
        fprintf(output_file, "Expected number of false positives: E = %.4f\n", E);
        // Debug breakdown to verify parameters and contributions
        print_E_breakdown(output_file, num_classes, class_sizes, N_val, P_succ, capacity);
        fprintf(output_file, "Type I error with threshold θ=1: P(D≥1) ≈ 1-exp(-E) = %.6f\n", 1.0 - exp(-E));
        fprintf(output_file, "Recommended balanced threshold: θ = ⌈E⌉ = %d\n", (int)ceil(E));
        fprintf(output_file, "Recommended conservative threshold: θ = ⌈E+√E⌉ = %d\n", (int)ceil(E + sqrt(E)));
        
        free(class_sizes);
    }
    
    fprintf(output_file, "%s", separator);
    
    // Print summary to terminal
    printf("\n%s", separator);
    printf("Attack Summary:\n");
    printf("===============\n");
    printf("Elapsed time: %.2f seconds\n", elapsed_time);
    printf("left involved key bits: [");
    for (int i = 0; i < num_left; i++) {
        printf("%d", left_ik[i]);
        if (i < num_left - 1) printf(", ");
    }
    printf("]\n");
    printf("right involved key bits: [");
    for (int i = 0; i < num_right; i++) {
        printf("%d", right_ik[i]);
        if (i < num_right - 1) printf(", ");
    }
    printf("]\n");
    printf("Common key bits: {");
    for (int i = 0; i < num_common; i++) {
        printf("%d", common_ik[i]);
        if (i < num_common - 1) printf(", ");
    }
    printf("}\n");
    printf("Total number of keys: %llu\n", (unsigned long long)total_number_of_keys);
    printf("Number of weak keys: %d\n", cnt);
    printf("Number of unique indices for the hash table (indices are sorted good pairs): %d\n", num_classes);
    printf("Size of union of all good pairs: %d\n", total_unique_pairs);
    printf("The conditions holds for %d keys out of %llu keys.\n", cnt, (unsigned long long)total_number_of_keys);
    
    if (cnt > 0) {
        double log2_weak_key_space = log2((double)cnt / total_number_of_keys);
        printf("The attack works for 2^(%.2f) portion of the key space.\n", log2_weak_key_space);
    }
    
    if ((uint64_t)cnt < total_number_of_keys) {
        double learned_bits = log2((double)(total_number_of_keys - cnt) / total_number_of_keys);
        printf("Number of learned bits in the worst case (when the key is not weak): %.2f\n", fabs(learned_bits));
    }
    
    // Print statistical analysis if parameters were provided
    if (log2_N >= 0 && P_succ > 0 && capacity_exp > 0) {
        // Calculate capacity and standard deviations from exponent
        double capacity = pow(2.0, -capacity_exp);
        double sigma0 = sqrt(capacity);
        double sigma1 = sqrt(capacity);
        
        // Collect class sizes (|X_i| = number of good pairs for each class)
        int *class_sizes = (int *)malloc(num_classes * sizeof(int));
        for (int i = 0; i < num_classes; i++) {
            class_sizes[i] = key_classes[i].num_pairs;
        }
        
        double N_val = pow(2.0, log2_N);
        double E = compute_E(num_classes, class_sizes, N_val, P_succ, capacity, sigma0, sigma1);
        
        printf("\nStatistical Analysis:\n");
        printf("=====================\n");
        printf("N = 2^%d = %.0f (number of plaintext-ciphertext pairs)\n", log2_N, N_val);
        printf("P_succ = %.2f (target detection probability)\n", P_succ);
        printf("capacity = 2^(-%.2f) = %.10f (distinguisher capacity)\n", capacity_exp, capacity);
        printf("sigma0 = sigma1 = sqrt(capacity) = %.6f (standard deviations)\n", sigma0);
        printf("Expected number of false positives: E = %.4f\n", E);
        // Mirror debug breakdown to terminal
        print_E_breakdown(stdout, num_classes, class_sizes, N_val, P_succ, capacity);
        printf("Type I error with threshold θ=1: P(D≥1) ≈ 1-exp(-E) = %.6f\n", 1.0 - exp(-E));
        printf("Recommended balanced threshold: θ = ⌈E⌉ = %d\n", (int)ceil(E));
        printf("Recommended conservative threshold: θ = ⌈E+√E⌉ = %d\n", (int)ceil(E + sqrt(E)));
        
        free(class_sizes);
    }
    
    printf("%s", separator);
    printf("\nDetailed results saved to: %s\n", output_filename);
    
    for (int i = 0; i < num_classes; i++) {
        free(key_classes[i].keys);
    }
    free(key_classes);
    fclose(output_file);
    
    return 0;
}
