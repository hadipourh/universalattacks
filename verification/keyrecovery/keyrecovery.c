/*
* Copyright (C) 2025 Hosein Hadipour and Mostafizar Rahman
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
*/

#include<time.h>
#include<math.h>
#include<stdint.h>
#include<stdio.h>     // For snprintf, printf, etc.
#include<stdlib.h>    // For malloc, atoi
#include <sys/stat.h>  // For mkdir
#include <errno.h> 
#include <string.h>    // For strrchr, strcpy

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    uint8_t left;
    uint8_t right;
} SubKeyPair;

typedef struct {
    uint8_t x;
    uint8_t y;
} BytePair;

typedef struct {
    SubKeyPair *keys;
    size_t key_count;
    BytePair *pairs;
    size_t pair_count;
    uint64_t total_key_space;
} KeyAnalysisData;

#define KEYRECOVERY_SHARED_TYPES

int initialize_key_analysis_data(int offset,
                                 const int *active_indices,
                                 int num_active,
                                 uint8_t dy1,
                                 uint8_t dy2,
                                 int use_unordered,
                                 KeyAnalysisData *out_data);
void free_key_analysis_data(KeyAnalysisData *data);

void kr_set_enumeration_data(const SubKeyPair *keys,
                             size_t key_count,
                             const BytePair *pairs,
                             size_t pair_count,
                             const int *active_indices,
                             int num_active);
void kr_release_enumeration_data(void);
void kr_set_current_mode(int mode);
void kr_set_dy_bytes(uint8_t dy1, uint8_t dy2);

// Function to create directories recursively
int mkdir_p(const char *path) {
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
        if (*p == '/') {
            *p = 0;
            if (mkdir(tmp, 0755) == -1 && errno != EEXIST)
                return -1;
            *p = '/';
        }
    if (mkdir(tmp, 0755) == -1 && errno != EEXIST)
        return -1;
    return 0;
}

#include "orthrosmodkey.c"
#include "attackutils.c"

const int STATE_SIZE=32;


int main(int argc, char *argv[]) {
    if (argc < 14 || argc > 15) {
        fprintf(stderr,
                "Usage: %s num_of_rounds rndOffset outputmask inputDiffPattern "
                "keyValType plainNibValType DEG N trailFolder dy1_hi dy1_lo dy2_hi dy2_lo [num_threads]\n",
                argv[0]);
        fprintf(stderr, "  num_threads: optional, defaults to auto-detection (use all available cores)\n");
        return 1;
    }
    
    clock_t start_time = clock();
    int num_of_rounds = atoi(argv[1]);
    int rndOffset = atoi(argv[2]);
    char outputmask[33];
    strncpy(outputmask, argv[3], 32);
    outputmask[32] = '\0';

    int inputDiffPattern[32] = {0};
    char *token = strtok(argv[4], ",");
    for (int i = 0; i < 32 && token != NULL; ++i) {
        inputDiffPattern[i] = atoi(token);
        token = strtok(NULL, ",");
    }

    int keyValType = atoi(argv[5]);
    int plainNibValType = atoi(argv[6]);  // Currently unused but parsed for compatibility
    (void)plainNibValType;
    int DEG = atoi(argv[7]);
    int N = atoi(argv[8]);

    char trailFolder[400];  // Large buffer for trail folder name
    strncpy(trailFolder, argv[9], sizeof(trailFolder) - 1);
    trailFolder[sizeof(trailFolder) - 1] = '\0';  // Ensure null termination

    int dy1_hi = (int)strtol(argv[10], NULL, 16);
    int dy1_lo = (int)strtol(argv[11], NULL, 16);
    int dy2_hi = (int)strtol(argv[12], NULL, 16);
    int dy2_lo = (int)strtol(argv[13], NULL, 16);

    if ((dy1_hi < 0 || dy1_hi > 0xF) || (dy1_lo < 0 || dy1_lo > 0xF) ||
        (dy2_hi < 0 || dy2_hi > 0xF) || (dy2_lo < 0 || dy2_lo > 0xF)) {
        fprintf(stderr, "Error: dy values must be 4-bit hex numbers (0x0 .. 0xF)\n");
        return 1;
    }

    uint8_t dy1_byte = (uint8_t)(((dy1_hi & 0xF) << 4) | (dy1_lo & 0xF));
    uint8_t dy2_byte = (uint8_t)(((dy2_hi & 0xF) << 4) | (dy2_lo & 0xF));

    // Handle optional num_threads parameter
    int requested_threads = 0;  // 0 means auto-detect
    if (argc == 15) {
        requested_threads = atoi(argv[14]);
        if (requested_threads < 0) {
            fprintf(stderr, "Error: num_threads must be positive (or 0 for auto)\n");
            return 1;
        }
    }
    
    // Configure OpenMP threads
    #ifdef _OPENMP
    if (requested_threads > 0) {
        omp_set_num_threads(requested_threads);
        printf("OpenMP: Using %d threads (user-specified)\n", requested_threads);
    } else {
        // Auto-detect: use all available processors
        int num_procs = omp_get_num_procs();
        omp_set_num_threads(num_procs);
        printf("OpenMP: Using %d threads (auto-detected)\n", num_procs);
    }
    #else
    if (requested_threads > 0) {
        printf("Warning: num_threads specified but OpenMP not enabled. Running serially.\n");
    }
    printf("Running in serial mode (OpenMP not enabled)\n");
    #endif

    int active_indices[32];
    int num_active = 0;
    for (int i = 0; i < 32; ++i) {
        if (inputDiffPattern[i] == 1) {
            if (num_active < 32) {
                active_indices[num_active++] = i;
            }
        }
    }

    if (num_active < 1) {
        fprintf(stderr, "Error: input difference pattern must activate at least one S-box\n");
        return 1;
    }

    KeyAnalysisData analysis_data = {0};
    if (initialize_key_analysis_data(rndOffset,
                                     active_indices,
                                     num_active,
                                     dy1_byte,
                                     dy2_byte,
                                     /*use_unordered=*/1,
                                     &analysis_data) != 0) {
        fprintf(stderr, "Error: failed to generate weak-key information\n");
        return 1;
    }

    printf("Key-analysis summary: weak keys = %zu, union of good pairs = %zu, total key space = %llu\n",
           analysis_data.key_count,
           analysis_data.pair_count,
           (unsigned long long)analysis_data.total_key_space);

    kr_set_enumeration_data(analysis_data.keys,
                            analysis_data.key_count,
                            analysis_data.pairs,
                            analysis_data.pair_count,
                            active_indices,
                            num_active);
    kr_set_dy_bytes(dy1_byte, dy2_byte);

    // Create directory name
    char dir_name[512];  // Large buffer to prevent any truncation warnings
    snprintf(dir_name, sizeof(dir_name), "data/%s/output%d", trailFolder, DEG);

    // Create the directory recursively
    if (mkdir_p(dir_name) == -1) {
        perror("Error creating directory");
        kr_release_enumeration_data();
        free_key_analysis_data(&analysis_data);
        return 1;
    }

    int activeNib = num_active;
    int (*keyBitPos)[8] = malloc(activeNib * sizeof(*keyBitPos));
    if (!keyBitPos) {
        fprintf(stderr, "Error: memory allocation failed for keyBitPos\n");
        kr_release_enumeration_data();
        free_key_analysis_data(&analysis_data);
        return 1;
    }

    keyNibs = malloc(activeNib * sizeof(*keyNibs));
    if (!keyNibs) {
        fprintf(stderr, "Error: memory allocation failed for keyNibs\n");
        free(keyBitPos);
        kr_release_enumeration_data();
        free_key_analysis_data(&analysis_data);
        return 1;
    }

    const int mode = 2; //0 for left, 1 for right, 2 for prf

    uint64_t N1 = 1ULL << DEG; // Number of queries:  N1 = 2^(DEG)
    srand(time(NULL));
    (void)init_prng(((unsigned int)rand() << 16) | rand());  // Initialize PRNG, return value unused
    getKeyBitPosNew(inputDiffPattern, keyBitPos, activeNib, rndOffset);

    unsigned char outMask[32];
    hex_string_to_bytes(outputmask,outMask,sizeof(outMask));
    if((mode==modeLeftBr) && (mode==modeRightBr)){
        printf("Error: Set PRF");
        return 1;
    } 

    // Init the file name for storing the data based on the key choice
    //used only for good key nibbles
    char file_name_corr[100];
    

    for(int num_of_experiments=0;num_of_experiments<N;num_of_experiments++){
        printf("\n----------------------------------------------------------------------------------------------------");
        printf("\n Exp No: %d\n", num_of_experiments);
        strcpy(file_name_corr, dir_name);
        strcat(file_name_corr, "/out");
        (void)init_prng(233);  // Initialize PRNG, return value unused
        unsigned char key[32];
        kr_set_current_mode(keyValType);
        fixKey(key, keyBitPos, keyValType, activeNib, file_name_corr);
        runExpGoodPairs(key,outMask,activeNib,N1,inputDiffPattern,rndOffset,num_of_rounds, file_name_corr);
    }

    //free dynamic memories
    kr_release_enumeration_data();
    free_key_analysis_data(&analysis_data);
    free(keyNibs);
    free(keyBitPos);

    clock_t end_time = clock();
    // Calculate the runtime in seconds
double runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\nRuntime: %f seconds\n", runtime);
    // Normal successful termination
    return 0;
}

// -----------------------------------------------------------------------------
// Key analysis helpers (adapted from keyanalysis.cpp)
// -----------------------------------------------------------------------------

#define MAX_GOOD_PAIRS 1024

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

static const uint8_t sb4_table[16] = {
    0x1, 0x0, 0x2, 0x4, 0x3, 0x8, 0x6, 0xd,
    0x9, 0xa, 0xb, 0xe, 0xf, 0xc, 0x7, 0x5
};

static uint8_t sb8_table[256];
static int sb8_initialized = 0;

typedef struct {
    int left;
    int right;
} CommonIndex;

typedef struct {
    BytePair pairs[MAX_GOOD_PAIRS];
    int count;
} GoodPairList;

static void init_sb8(void) {
    if (sb8_initialized) return;
    for (int a = 0; a < 16; a++) {
        for (int b = 0; b < 16; b++) {
            sb8_table[a * 16 + b] = (uint8_t)((sb4_table[a] << 4) | sb4_table[b]);
        }
    }
    sb8_initialized = 1;
}

static void derive_involved_key_positions(int offset,
                                          const int *active_indices,
                                          int num_active,
                                          int *left_ik,
                                          int *right_ik,
                                          int *num_left,
                                          int *num_right,
                                          int *common_ik,
                                          int *num_common,
                                          CommonIndex *common_indices) {
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
        int temp_left[128];
        int temp_right[128];
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

static void check_candidate(uint8_t k1,
                            uint8_t k2,
                            uint8_t dx,
                            uint8_t dy1,
                            uint8_t dy2,
                            GoodPairList *good_pairs,
                            int use_unordered) {
    for (int x = 0; x < 256; x++) {
        uint8_t y = (uint8_t)(x ^ dx);

        if (use_unordered && x >= y) {
            continue;
        }

        if ((sb8_table[x ^ k1] ^ sb8_table[y ^ k1]) == dy1 &&
            (sb8_table[x ^ k2] ^ sb8_table[y ^ k2]) == dy2) {
            int duplicate = 0;
            for (int i = 0; i < good_pairs->count; i++) {
                if (good_pairs->pairs[i].x == (uint8_t)x &&
                    good_pairs->pairs[i].y == (uint8_t)y) {
                    duplicate = 1;
                    break;
                }
            }
            if (!duplicate && good_pairs->count < MAX_GOOD_PAIRS) {
                good_pairs->pairs[good_pairs->count].x = (uint8_t)x;
                good_pairs->pairs[good_pairs->count].y = (uint8_t)y;
                good_pairs->count++;
            }
        }
    }
}

static int compare_byte_pairs(const void *a, const void *b) {
    const BytePair *pa = (const BytePair *)a;
    const BytePair *pb = (const BytePair *)b;
    if (pa->x != pb->x) return (int)pa->x - (int)pb->x;
    return (int)pa->y - (int)pb->y;
}

static int compare_subkey_pairs(const void *a, const void *b) {
    const SubKeyPair *ka = (const SubKeyPair *)a;
    const SubKeyPair *kb = (const SubKeyPair *)b;
    if (ka->left != kb->left) return (int)ka->left - (int)kb->left;
    return (int)ka->right - (int)kb->right;
}

int initialize_key_analysis_data(int offset,
                                 const int *active_indices,
                                 int num_active,
                                 uint8_t dy1,
                                 uint8_t dy2,
                                 int use_unordered,
                                 KeyAnalysisData *out_data) {
    if (!out_data) return -1;
    memset(out_data, 0, sizeof(*out_data));

    init_sb8();

    if (num_active != 2) {
        fprintf(stderr, "Error: initialize_key_analysis_data expects exactly 2 active S-boxes (got %d)\n", num_active);
        return -1;
    }

    int left_ik[256];
    int right_ik[256];
    int common_ik[256];
    CommonIndex common_indices[256];
    int num_left = 0;
    int num_right = 0;
    int num_common = 0;

    derive_involved_key_positions(offset,
                                  active_indices,
                                  num_active,
                                  left_ik,
                                  right_ik,
                                  &num_left,
                                  &num_right,
                                  common_ik,
                                  &num_common,
                                  common_indices);

    uint64_t total_keys = 1ULL << (num_left + num_right - num_common);
    out_data->total_key_space = total_keys;

    size_t key_capacity = 0;
    size_t pair_capacity = 0;
    uint8_t pair_seen[256][256] = {{0}};

    for (int k1 = 0; k1 < 256; k1++) {
        for (int k2 = 0; k2 < 256; k2++) {
            int valid = 1;
            for (int i = 0; i < num_common; i++) {
                int ls = 8 - common_indices[i].left - 1;
                int rs = 8 - common_indices[i].right - 1;
                if (((k1 >> ls) & 0x1) != ((k2 >> rs) & 0x1)) {
                    valid = 0;
                    break;
                }
            }
            if (!valid) continue;

            GoodPairList good_pairs;
            good_pairs.count = 0;
            for (int dx = 1; dx < 256; dx++) {
                check_candidate((uint8_t)k1,
                                (uint8_t)k2,
                                (uint8_t)dx,
                                dy1,
                                dy2,
                                &good_pairs,
                                use_unordered);
            }

            if (good_pairs.count > 0) {
                if (out_data->key_count == key_capacity) {
                    size_t new_cap = key_capacity ? key_capacity * 2 : 128;
                    SubKeyPair *tmp = (SubKeyPair *)realloc(out_data->keys, new_cap * sizeof(SubKeyPair));
                    if (!tmp) {
                        free(out_data->keys);
                        free(out_data->pairs);
                        memset(out_data, 0, sizeof(*out_data));
                        return -1;
                    }
                    out_data->keys = tmp;
                    key_capacity = new_cap;
                }
                out_data->keys[out_data->key_count].left = (uint8_t)k1;
                out_data->keys[out_data->key_count].right = (uint8_t)k2;
                out_data->key_count++;

                for (int i = 0; i < good_pairs.count; i++) {
                    uint8_t x = good_pairs.pairs[i].x;
                    uint8_t y = good_pairs.pairs[i].y;
                    if (!pair_seen[x][y]) {
                        if (out_data->pair_count == pair_capacity) {
                            size_t new_cap = pair_capacity ? pair_capacity * 2 : 256;
                            BytePair *tmp = (BytePair *)realloc(out_data->pairs, new_cap * sizeof(BytePair));
                            if (!tmp) {
                                free(out_data->keys);
                                free(out_data->pairs);
                                memset(out_data, 0, sizeof(*out_data));
                                return -1;
                            }
                            out_data->pairs = tmp;
                            pair_capacity = new_cap;
                        }
                        out_data->pairs[out_data->pair_count] = good_pairs.pairs[i];
                        out_data->pair_count++;
                        pair_seen[x][y] = 1;
                    }
                }
            }
        }
    }

    if (out_data->key_count > 1) {
        qsort(out_data->keys, out_data->key_count, sizeof(SubKeyPair), compare_subkey_pairs);
    }
    if (out_data->pair_count > 1) {
        qsort(out_data->pairs, out_data->pair_count, sizeof(BytePair), compare_byte_pairs);
    }

    return 0;
}

void free_key_analysis_data(KeyAnalysisData *data) {
    if (!data) return;
    free(data->keys);
    free(data->pairs);
    data->keys = NULL;
    data->pairs = NULL;
    data->key_count = 0;
    data->pair_count = 0;
    data->total_key_space = 0;
}
