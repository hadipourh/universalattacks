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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __APPLE__
    #include <sys/random.h>
#elif defined(__linux__)
    #include <sys/random.h>
#else
    #include <fcntl.h>
    #include <unistd.h>
#endif

#ifndef KEYRECOVERY_SHARED_TYPES
#define KEYRECOVERY_SHARED_TYPES
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
#endif

typedef struct {
    uint64_t state[4];  // 256-bit state
} thread_rng_state;

#ifndef HAVE_SYSTEM_ENTROPY_HELPERS
#define HAVE_SYSTEM_ENTROPY_HELPERS
static void get_system_entropy(void *buffer, size_t length) {
#ifdef __APPLE__
    if (getentropy(buffer, length) != 0) {
        FILE *urandom = fopen("/dev/urandom", "rb");
        if (urandom) {
            size_t read_bytes = fread(buffer, 1, length, urandom);
            (void)read_bytes;
            fclose(urandom);
        }
    }
#elif defined(__linux__)
    if (getrandom(buffer, length, 0) < 0) {
        FILE *urandom = fopen("/dev/urandom", "rb");
        if (urandom) {
            size_t read_bytes = fread(buffer, 1, length, urandom);
            (void)read_bytes;
            fclose(urandom);
        }
    }
#else
    FILE *urandom = fopen("/dev/urandom", "rb");
    if (urandom) {
        size_t read_bytes = fread(buffer, 1, length, urandom);
        (void)read_bytes;
        fclose(urandom);
    }
#endif
}
#endif

static void init_thread_rng(thread_rng_state *rng, unsigned int seed) {
    get_system_entropy(rng->state, sizeof(rng->state));
    rng->state[0] ^= (uint64_t)seed;
#ifdef _OPENMP
    rng->state[1] ^= (uint64_t)omp_get_thread_num();
#endif
    rng->state[2] ^= (uint64_t)time(NULL);
    rng->state[3] ^= (uint64_t)clock();
}

static inline uint32_t thread_rand(thread_rng_state *rng) {
    const uint64_t result = rng->state[0] + rng->state[3];
    const uint64_t t = rng->state[1] << 17;

    rng->state[2] ^= rng->state[0];
    rng->state[3] ^= rng->state[1];
    rng->state[1] ^= rng->state[2];
    rng->state[0] ^= rng->state[3];
    rng->state[2] ^= t;

    rng->state[3] = (rng->state[3] << 45) | (rng->state[3] >> 19);

    return (uint32_t)(result >> 32);
}

const int GOOD = 0;
const int RANDOM = 1;
const int BAD = 2;

unsigned char (*keyNibs)[2] = NULL;

static const SubKeyPair *g_weak_keys = NULL;
static size_t g_weak_key_count = 0;
static uint8_t *g_weak_key_bitmap = NULL;  // 256 * 256 bitmap

static const BytePair *g_good_pairs = NULL;
static size_t g_good_pair_count = 0;

static int g_active_nibbles[32];
static int g_num_active_nibbles = 0;

static SubKeyPair g_current_key_pair = {0, 0};
static int g_current_key_mode = RANDOM;
static int g_should_emit_good_pairs = 0;
static BytePair *g_selected_good_pairs = NULL;
static size_t g_selected_good_pair_count = 0;
static uint8_t g_dy1_byte = 0;
static uint8_t g_dy2_byte = 0;

static const uint8_t sb4_table_local[16] = {
    0x1, 0x0, 0x2, 0x4, 0x3, 0x8, 0x6, 0xd,
    0x9, 0xa, 0xb, 0xe, 0xf, 0xc, 0x7, 0x5
};
static uint8_t sb8_table_local[256];
static int sb8_local_initialized = 0;

static void ensure_sb8_local(void) {
    if (sb8_local_initialized) return;
    for (int a = 0; a < 16; a++) {
        for (int b = 0; b < 16; b++) {
            sb8_table_local[a * 16 + b] = (uint8_t)((sb4_table_local[a] << 4) | sb4_table_local[b]);
        }
    }
    sb8_local_initialized = 1;
}

static void clear_selected_pairs(void) {
    free(g_selected_good_pairs);
    g_selected_good_pairs = NULL;
    g_selected_good_pair_count = 0;
    g_should_emit_good_pairs = 0;
}

static void populate_selected_pairs(uint8_t k_left, uint8_t k_right) {
    clear_selected_pairs();
    ensure_sb8_local();

    size_t capacity = 64;
    BytePair *pairs = (BytePair *)malloc(capacity * sizeof(BytePair));
    if (!pairs) {
        return;
    }

    size_t count = 0;
    for (int dx = 1; dx < 256; dx++) {
        for (int x = 0; x < 256; x++) {
            uint8_t y = (uint8_t)(x ^ dx);
            if (x >= y) continue;

            uint8_t left_diff = (uint8_t)(sb8_table_local[x ^ k_left] ^ sb8_table_local[y ^ k_left]);
            uint8_t right_diff = (uint8_t)(sb8_table_local[x ^ k_right] ^ sb8_table_local[y ^ k_right]);

            if (left_diff == g_dy1_byte && right_diff == g_dy2_byte) {
                if (count == capacity) {
                    capacity *= 2;
                    BytePair *tmp = (BytePair *)realloc(pairs, capacity * sizeof(BytePair));
                    if (!tmp) {
                        free(pairs);
                        return;
                    }
                    pairs = tmp;
                }
                pairs[count].x = (uint8_t)x;
                pairs[count].y = y;
                count++;
            }
        }
    }

    if (count == 0) {
        free(pairs);
        return;
    }

    g_selected_good_pairs = pairs;
    g_selected_good_pair_count = count;
    g_should_emit_good_pairs = 1;
}

static inline uint8_t nib_high(uint8_t value) {
    return (uint8_t)((value >> 4) & 0xF);
}

static inline uint8_t nib_low(uint8_t value) {
    return (uint8_t)(value & 0xF);
}

void kr_set_enumeration_data(const SubKeyPair *keys,
                             size_t key_count,
                             const BytePair *pairs,
                             size_t pair_count,
                             const int *active_indices,
                             int num_active) {
    g_weak_keys = keys;
    g_weak_key_count = key_count;
    g_good_pairs = pairs;
    g_good_pair_count = pair_count;
    g_num_active_nibbles = num_active;

    if (num_active > 0 && active_indices) {
        memcpy(g_active_nibbles, active_indices, (size_t)num_active * sizeof(int));
    } else {
        g_num_active_nibbles = 0;
    }

    free(g_weak_key_bitmap);
    g_weak_key_bitmap = NULL;
    if (g_weak_key_count > 0) {
        g_weak_key_bitmap = (uint8_t *)calloc(256 * 256, sizeof(uint8_t));
        if (g_weak_key_bitmap) {
            for (size_t i = 0; i < g_weak_key_count; ++i) {
                const SubKeyPair *kp = &g_weak_keys[i];
                g_weak_key_bitmap[((size_t)kp->left << 8) | kp->right] = 1;
            }
        }
    }
}

void kr_set_current_mode(int mode) {
    g_current_key_mode = mode;
}

void kr_set_dy_bytes(uint8_t dy1, uint8_t dy2) {
    g_dy1_byte = dy1;
    g_dy2_byte = dy2;
    ensure_sb8_local();
}

void kr_release_enumeration_data(void) {
    g_weak_keys = NULL;
    g_weak_key_count = 0;
    g_good_pairs = NULL;
    g_good_pair_count = 0;
    g_num_active_nibbles = 0;
    free(g_weak_key_bitmap);
    g_weak_key_bitmap = NULL;
    g_current_key_mode = RANDOM;
    clear_selected_pairs();
}

void fixKey(unsigned char key[32],
            int (*keyBitPos)[8],
            int keyValType,
            int numActNib,
            char *fileName) {
    for (int i = 0; i < 32; i++) {
        key[i] = (unsigned char)(rand() & 0xF);
    }

    SubKeyPair pair = { (uint8_t)(rand() & 0xFF), (uint8_t)(rand() & 0xFF) };

    if (keyValType == GOOD && g_weak_key_count > 0) {
        size_t idx = (size_t)(rand() % (int)g_weak_key_count);
        pair = g_weak_keys[idx];
    } else if (keyValType == BAD && g_weak_key_count > 0 && g_weak_key_bitmap) {
        do {
            pair.left = (uint8_t)(rand() & 0xFF);
            pair.right = (uint8_t)(rand() & 0xFF);
        } while (g_weak_key_bitmap[((size_t)pair.left << 8) | pair.right]);
    }

    int pair_is_weak = (g_weak_key_bitmap && g_weak_key_bitmap[((size_t)pair.left << 8) | pair.right]) ? 1 : 0;
    if (keyValType == GOOD || (keyValType == RANDOM && pair_is_weak)) {
        populate_selected_pairs(pair.left, pair.right);
    } else {
        clear_selected_pairs();
    }

    if (keyNibs) {
        for (int idx = 0; idx < numActNib; ++idx) {
            uint8_t leftNib;
            uint8_t rightNib;
            if (idx == 0) {
                leftNib = nib_high(pair.left);
                rightNib = nib_high(pair.right);
            } else if (idx == 1) {
                leftNib = nib_low(pair.left);
                rightNib = nib_low(pair.right);
            } else {
                leftNib = (uint8_t)(rand() & 0xF);
                rightNib = (uint8_t)(rand() & 0xF);
            }
            keyNibs[idx][0] = leftNib;
            keyNibs[idx][1] = rightNib;
        }
    }

    if (fileName && keyNibs) {
        for (int index = 0; index < numActNib; index++) {
            char buffer[20];
            snprintf(buffer, sizeof(buffer), "_%d-%d", keyNibs[index][0], keyNibs[index][1]);
            strcat(fileName, buffer);
        }
        strcat(fileName, ".txt");
    }

    setNibbles(key, keyBitPos, keyNibs, numActNib);

    g_current_key_pair = pair;
    g_should_emit_good_pairs = (g_selected_good_pair_count > 0);
}

void keyExpansion(unsigned char key[32], unsigned char expdKey1[][32], unsigned char expdKey2[][32], int num){
    unsigned char key1[32], key2[32];
    for(int nibIndex=0;nibIndex<32;nibIndex++){
        key1[nibIndex]=key[nibIndex];
        key2[nibIndex]=key[nibIndex];
    }
    for(int rndNum=0;rndNum<num;rndNum++){
        for(int nibIndex=0;nibIndex<32;nibIndex++){
            expdKey1[rndNum][nibIndex]=key1[nibIndex];
            expdKey2[rndNum][nibIndex]=key2[nibIndex];
        }
        KEYSCHEDULING(expdKey1[rndNum],1);
        KEYSCHEDULING(expdKey2[rndNum],2);
        for(int nibIndex=0;nibIndex<32;nibIndex++){
            key1[nibIndex]=expdKey1[rndNum][nibIndex];
            key2[nibIndex]=expdKey2[rndNum][nibIndex];
        }
    }
}

void runExpGoodPairs(unsigned char key[32],
                     unsigned char outMask[32],
                     int numActNib,
                     uint64_t DEGREE,
                     int inpDifffPattern[32],
                     int roundOffset,
                     int num_of_rounds,
                     char *fileNameCorr){
    int debug=0;
    uint64_t CORR;
    FILE *fileCorrData = fopen(fileNameCorr, "w");
    if (fileCorrData == NULL) {
        printf("Failed to open the file.\n");
        return;
    }

    if (g_should_emit_good_pairs && g_selected_good_pairs && g_selected_good_pair_count > 0) {
        fprintf(fileCorrData, "# KEY left=0x%02x right=0x%02x\n",
                g_current_key_pair.left,
                g_current_key_pair.right);
        fprintf(fileCorrData, "# GOOD_PAIRS_BEGIN count=%zu\n", g_selected_good_pair_count);
        for (size_t i = 0; i < g_selected_good_pair_count; ++i) {
            fprintf(fileCorrData, "# GOOD_PAIR 0x%02x 0x%02x\n",
                    g_selected_good_pairs[i].x,
                    g_selected_good_pairs[i].y);
        }
        fprintf(fileCorrData, "# GOOD_PAIRS_END\n");
    }

    if (g_good_pair_count == 0) {
        fclose(fileCorrData);
        printf("Warning: no good pairs available. Skipping experiments.\n");
        return;
    }

    unsigned char expdKey1[roundOffset+num_of_rounds+1][32],expdKey2[roundOffset+num_of_rounds+1][32];
    keyExpansion(key,expdKey1,expdKey2,roundOffset+num_of_rounds+1);
    printf("\ntot=%zu", g_good_pair_count);
    printf("\nProcessing pairs: ");
    fflush(stdout);

    for(size_t index=0; index<g_good_pair_count; index++){
        // Progress reporting every 10% or every 100 pairs, whichever is more frequent
        if (index % ((g_good_pair_count / 10) + 1) == 0 || index % 100 == 0) {
            printf("Processing pairs: %zu/%zu (%.1f%%)\n", index, g_good_pair_count, 
                   (double)index * 100.0 / g_good_pair_count);
            fflush(stdout);
        }
        
        if(debug){printf("\nIteration %zu\n", index);}
        unsigned char plaintext1[32] = {0};
        unsigned char plaintext2[32] = {0};

        BytePair pair = g_good_pairs[index];

        for (int act = 0; act < g_num_active_nibbles && act < numActNib; ++act) {
            int pos = g_active_nibbles[act];
            uint8_t px = 0;
            uint8_t py = 0;
            if (act == 0) {
                px = nib_high(pair.x);
                py = nib_high(pair.y);
            } else if (act == 1) {
                px = nib_low(pair.x);
                py = nib_low(pair.y);
            }
            plaintext1[pos] = px;
            plaintext2[pos] = py;
        }


        uint64_t counter0 = 0ULL;
        uint64_t counter1 = 0ULL; 

        #pragma omp parallel
        {
            thread_rng_state local_rng;
            uint64_t local_counter0 = 0;
            uint64_t local_counter1 = 0;
            unsigned char local_plaintext1[32];
            unsigned char local_plaintext2[32];
            unsigned char local_ciphertext1[32];
            unsigned char local_ciphertext2[32];

            #ifdef _OPENMP
            init_thread_rng(&local_rng, (unsigned int)time(NULL) + omp_get_thread_num() * 12345);
            #else
            init_thread_rng(&local_rng, (unsigned int)time(NULL));
            #endif

            for(int i=0; i<32; i++){
                if(inpDifffPattern[i]==1){
                    local_plaintext1[i] = plaintext1[i];
                    local_plaintext2[i] = plaintext2[i];
                }
            }

            #pragma omp for schedule(dynamic, 1024)
            for (uint64_t loopcnt = 0; loopcnt < DEGREE; loopcnt++){
                for(int i=0; i<32; i++){
                    if(inpDifffPattern[i]==0){
                        unsigned char random_nib = (unsigned char)(thread_rand(&local_rng) & 0xF);
                        local_plaintext1[i] = random_nib;
                        local_plaintext2[i] = random_nib;
                    }
                }

                ORTHROS(local_plaintext1, local_ciphertext1, expdKey1, expdKey2, roundOffset, num_of_rounds);
                ORTHROS(local_plaintext2, local_ciphertext2, expdKey1, expdKey2, roundOffset, num_of_rounds);

                if(dot_prod(local_ciphertext1, outMask, 32) == dot_prod(local_ciphertext2, outMask, 32)){
                    local_counter0++;
                } else {
                    local_counter1++;
                }
            }

            #pragma omp atomic
            counter0 += local_counter0;

            #pragma omp atomic
            counter1 += local_counter1;
        }

        if (counter0 > counter1)
            CORR = counter0 - counter1;
        else
            CORR = counter1 - counter0;

        if(debug){printf("Correlation counters: %llu %llu\n", (unsigned long long)counter0, (unsigned long long)counter1);}

        char strData[100]="";
        for(int actNibIndex=0;actNibIndex<numActNib && actNibIndex < g_num_active_nibbles;actNibIndex++){
            char temp[10];
            uint8_t px = (actNibIndex == 0) ? nib_high(pair.x) : nib_low(pair.x);
            uint8_t py = (actNibIndex == 0) ? nib_high(pair.y) : nib_low(pair.y);
            snprintf(temp, sizeof(temp),"%d %d ", px, py);
            strncat(strData, temp, sizeof(strData) - strlen(strData) - 1);
        }

        char corr_str[21];
        snprintf(corr_str, sizeof(corr_str), "%llu", (unsigned long long)CORR);
        strncat(strData, corr_str, sizeof(strData) - strlen(strData) - 1);

        write_string_to_file(fileCorrData,strData);

        if(debug){
            printf("Pair %zu written\n", index);
        }
    }
    printf("%zu/%zu (100.0%%) - Complete!\n", g_good_pair_count, g_good_pair_count);
    fclose(fileCorrData);
}

void getKeyBitPos(int *inpDiffPat,  int (*keyBitPos)[8]){
    int iter=0;
    for(int i=0;i<32;i++){
        if(inpDiffPat[i]==1){
            for(int j=0;j<4;j++){
                keyBitPos[iter][j]=invPermK1[4*i+j];
                keyBitPos[iter][4+j]=invPermK2[4*i+j];
            }
            iter++;
        }
    }
}

void getKeyBitPosNew(int *inpDiffPat,  int (*keyBitPos)[8], int numActiveNib, int rndOffset){
    int iter=0;
    int branch=2;
    int stateBits[numActiveNib][branch][128];
    int newStateBits[numActiveNib][branch][128];
    
    for(int i=0;i<numActiveNib;i++){
        for(int j=0;j<branch;j++){
            for(int k=0;k<128;k++){
                stateBits[i][j][k]=0;
                newStateBits[i][j][k]=0;
            }
        }
    }
    for(int i=0;i<32;i++){
        if(inpDiffPat[i]==1){
            for(int j=0;j<branch;j++){
                for(int k=0;k<128;k++){
                    if(j==0){
                        stateBits[iter][j][k]=(invPermK1[4*i]==k) || (invPermK1[4*i+1]==k) || (invPermK1[4*i+2]==k) || (invPermK1[4*i+3]==k);
                        if(stateBits[iter][j][k]==1){
                            stateBits[iter][j][k]=1+(k==invPermK1[4*i+0]?0:(k==invPermK1[4*i+1]?1:(k==invPermK1[4*i+2]?2:3)));
                        }
                    }
                    else{
                        stateBits[iter][j][k]=(invPermK2[4*i]==k) || (invPermK2[4*i+1]==k) || (invPermK2[4*i+2]==k) || (invPermK2[4*i+3]==k);
                        if(stateBits[iter][j][k]==1){
                            stateBits[iter][j][k]=1+(k==invPermK2[4*i+0]?0:(k==invPermK2[4*i+1]?1:(k==invPermK2[4*i+2]?2:3)));
                        }
                    }
                }
            }
            iter++;
        }
    }

    for(int rnd=0;rnd<rndOffset;rnd++){
        for(int i=0;i<numActiveNib;i++){
            for(int br=0;br<branch;br++){
                for(int bit=0;bit<128;bit++){
                    if(stateBits[i][br][bit]!=0){
                        if(br==0){
                            newStateBits[i][br][invPermK1[bit]]=stateBits[i][br][bit];
                        }
                        else{
                            newStateBits[i][br][invPermK2[bit]]=stateBits[i][br][bit];
                        }
                    }
                }
            }
        }
        for(int i=0;i<numActiveNib;i++){
            for(int br=0;br<branch;br++){
                for(int bit=0;bit<128;bit++){
                    stateBits[i][br][bit]=newStateBits[i][br][bit];
                    newStateBits[i][br][bit]=0;
                }
            }
        }
    }

    for(int nibNum=0;nibNum<numActiveNib;nibNum++){
        for(int bitNum=0;bitNum<128;bitNum++){
            if(stateBits[nibNum][0][bitNum]!=0){
                int order=stateBits[nibNum][0][bitNum]-1;
                keyBitPos[nibNum][order]=bitNum;
            }
        }
        for(int bitNum=0;bitNum<128;bitNum++){
            if(stateBits[nibNum][1][bitNum]!=0){
                int order=4+stateBits[nibNum][1][bitNum]-1;
                keyBitPos[nibNum][order]=bitNum;
            }
        }
    }
}
