#include <cuda_runtime.h>
#include <curand_kernel.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <sys/stat.h>

#ifdef __APPLE__
#include <sys/random.h>
#elif defined(__linux__)
#include <sys/random.h>
#else
#include <fcntl.h>
#include <unistd.h>
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

int mkdir_p(const char *path);

#include "orthrosmodkey.c"
#include "attackutils.c"

inline void cuda_check(cudaError_t err, const char *file, int line) {
	if (err != cudaSuccess) {
		fprintf(stderr, "CUDA error %s:%d: %s\n", file, line, cudaGetErrorString(err));
		exit(EXIT_FAILURE);
	}
}

#define CUDA_CHECK(expr) cuda_check((expr), __FILE__, __LINE__)

#define MAX_TOTAL_ROUNDS 64

__constant__ unsigned char d_RC_1[12][32];
__constant__ unsigned char d_RC_2[12][32];
__constant__ unsigned char d_permN1[32];
__constant__ unsigned char d_permN2[32];
__constant__ unsigned int d_permB1[128];
__constant__ unsigned int d_permB2[128];
__constant__ unsigned char d_SBOX[16];
__constant__ unsigned char d_expdKey1[MAX_TOTAL_ROUNDS][32];
__constant__ unsigned char d_expdKey2[MAX_TOTAL_ROUNDS][32];
__constant__ int d_inpDiffPattern[32];
__constant__ unsigned char d_outMask[32];

__device__ void NIBBLEtoBIT_device(const unsigned char *nibble, unsigned char *bit) {
	for (int i = 0; i < 128; ++i) {
		bit[i] = 0;
	}
	for (int i = 0; i < 32; ++i) {
		unsigned char v = nibble[i] & 0xF;
		bit[(4 * i) + 0] = (v >> 3) & 0x1;
		bit[(4 * i) + 1] = (v >> 2) & 0x1;
		bit[(4 * i) + 2] = (v >> 1) & 0x1;
		bit[(4 * i) + 3] = v & 0x1;
	}
}

__device__ void BITtoNIBBLE_device(unsigned char *nibble, const unsigned char *bit) {
	for (int i = 0; i < 32; ++i) {
		nibble[i] = (bit[(4 * i) + 0] << 3) ^ (bit[(4 * i) + 1] << 2) ^ (bit[(4 * i) + 2] << 1) ^ (bit[(4 * i) + 3] << 0);
	}
}

__device__ void BITPERM_device(unsigned char *state, int branch) {
	unsigned char bit[128];
	unsigned char bit_buf[128];
	NIBBLEtoBIT_device(state, bit);
	if (branch == 1) {
		for (int i = 0; i < 128; ++i) {
			bit_buf[d_permB1[i]] = bit[i];
		}
	} else {
		for (int i = 0; i < 128; ++i) {
			bit_buf[d_permB2[i]] = bit[i];
		}
	}
	BITtoNIBBLE_device(state, bit_buf);
}

__device__ void NIBBLEPERM_device(unsigned char *state, int branch) {
	unsigned char buf[32];
	for (int i = 0; i < 32; ++i) {
		buf[i] = state[i];
	}
	if (branch == 1) {
		for (int i = 0; i < 32; ++i) {
			state[d_permN1[i]] = buf[i];
		}
	} else {
		for (int i = 0; i < 32; ++i) {
			state[d_permN2[i]] = buf[i];
		}
	}
}

__device__ void ProcessBranch_device(unsigned char *state, const unsigned char rk[MAX_TOTAL_ROUNDS][32], int branch, int offset, int roundNum) {
	for (int i = 0; i < 32; ++i) {
		state[i] ^= rk[offset][i];
	}
	for (int r = offset; r < offset + roundNum; ++r) {
		for (int i = 0; i < 32; ++i) {
			state[i] = d_SBOX[state[i]];
		}

		if (r < 4) {
			BITPERM_device(state, branch);
		} else {
			NIBBLEPERM_device(state, branch);
		}

		unsigned char buf[32];
		for (int i = 0; i < 32; ++i) {
			buf[i] = state[i];
		}
		for (int i = 0; i < 8; ++i) {
			state[(i * 4) + 0] = buf[(i * 4) + 1] ^ buf[(i * 4) + 2] ^ buf[(i * 4) + 3];
			state[(i * 4) + 1] = buf[(i * 4) + 0] ^ buf[(i * 4) + 2] ^ buf[(i * 4) + 3];
			state[(i * 4) + 2] = buf[(i * 4) + 0] ^ buf[(i * 4) + 1] ^ buf[(i * 4) + 3];
			state[(i * 4) + 3] = buf[(i * 4) + 0] ^ buf[(i * 4) + 1] ^ buf[(i * 4) + 2];
		}

		if (branch == 1) {
			for (int i = 0; i < 32; ++i) {
				state[i] ^= rk[r + 1][i] ^ d_RC_1[r][i];
			}
		} else {
			for (int i = 0; i < 32; ++i) {
				state[i] ^= rk[r + 1][i] ^ d_RC_2[r][i];
			}
		}
	}
}

__device__ void ORTHROS_device(const unsigned char *plaintext, unsigned char *ciphertext, int offset, int roundNum) {
	unsigned char X1[32];
	unsigned char X2[32];

	for (int i = 0; i < 32; ++i) {
		X1[i] = plaintext[i];
		X2[i] = plaintext[i];
	}

	ProcessBranch_device(X1, d_expdKey1, 1, offset, roundNum);
	ProcessBranch_device(X2, d_expdKey2, 2, offset, roundNum);

	for (int i = 0; i < 32; ++i) {
		ciphertext[i] = X1[i] ^ X2[i];
	}
}

__device__ int dot_prod_device(const unsigned char *A, const unsigned char *B) {
	int output = 0;
	for (int i = 0; i < 32; ++i) {
		unsigned char C = A[i] & B[i];
		while (C > 0) {
			output ^= (C & 0x1);
			C >>= 1;
		}
	}
	return output;
}

__global__ void init_rng_kernel(curandStatePhilox4_32_10_t *rng_states,
				uint64_t seed,
				uint64_t totalThreads) {
	uint64_t global_tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
	if (global_tid >= totalThreads) {
		return;
	}
	curand_init(static_cast<unsigned long long>(seed),
		    static_cast<unsigned long long>(global_tid),
		    0ULL,
		    &rng_states[global_tid]);
}

__global__ void sampling_kernel(curandStatePhilox4_32_10_t *rng_states,
				const unsigned char *base_plaintext1,
				const unsigned char *base_plaintext2,
				int roundOffset,
				int numRounds,
				uint64_t totalThreads,
				uint64_t totalSamples,
				unsigned long long *counter0,
				unsigned long long *counter1) {
	uint64_t global_tid = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
	if (global_tid >= totalThreads) {
		return;
	}

	curandStatePhilox4_32_10_t local_rng = rng_states[global_tid];

	unsigned long long local_counter0 = 0;
	unsigned long long local_counter1 = 0;

	unsigned char local_plaintext1[32];
	unsigned char local_plaintext2[32];
	unsigned char local_ciphertext1[32];
	unsigned char local_ciphertext2[32];

	for (uint64_t iter = global_tid; iter < totalSamples; iter += totalThreads) {
		for (int i = 0; i < 32; ++i) {
			local_plaintext1[i] = base_plaintext1[i];
			local_plaintext2[i] = base_plaintext2[i];
		}

		for (int i = 0; i < 32; ++i) {
			if (d_inpDiffPattern[i] == 0) {
				unsigned char random_nib = static_cast<unsigned char>(curand(&local_rng) & 0xF);
				local_plaintext1[i] = random_nib;
				local_plaintext2[i] = random_nib;
			}
		}

		ORTHROS_device(local_plaintext1, local_ciphertext1, roundOffset, numRounds);
		ORTHROS_device(local_plaintext2, local_ciphertext2, roundOffset, numRounds);

		if (dot_prod_device(local_ciphertext1, d_outMask) == dot_prod_device(local_ciphertext2, d_outMask)) {
			++local_counter0;
		} else {
			++local_counter1;
		}
	}

	rng_states[global_tid] = local_rng;

	atomicAdd(counter0, local_counter0);
	atomicAdd(counter1, local_counter1);
}

static void ensure_device_constants_loaded() {
	static bool loaded = false;
	if (loaded) {
		return;
	}
	CUDA_CHECK(cudaMemcpyToSymbol(d_RC_1, RC_1, sizeof(RC_1)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_RC_2, RC_2, sizeof(RC_2)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_permN1, permN1, sizeof(permN1)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_permN2, permN2, sizeof(permN2)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_permB1, permB1, sizeof(permB1)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_permB2, permB2, sizeof(permB2)));
	CUDA_CHECK(cudaMemcpyToSymbol(d_SBOX, SBOX, sizeof(SBOX)));
	loaded = true;
}

int mkdir_p(const char *path) {
	if (!path || *path == '\0') {
		return -1;
	}

	size_t len = strlen(path);
	if (len == 0) {
		return -1;
	}

	char *tmp = (char *)malloc(len + 1);
	if (!tmp) {
		errno = ENOMEM;
		return -1;
	}

	memcpy(tmp, path, len + 1);
	if (tmp[len - 1] == '/') {
		tmp[len - 1] = '\0';
	}

	for (char *p = tmp + 1; *p; ++p) {
		if (*p == '/') {
			*p = '\0';
			if (mkdir(tmp, 0755) == -1 && errno != EEXIST) {
				int saved_errno = errno;
				free(tmp);
				errno = saved_errno;
				return -1;
			}
			*p = '/';
		}
	}

	if (mkdir(tmp, 0755) == -1 && errno != EEXIST) {
		int saved_errno = errno;
		free(tmp);
		errno = saved_errno;
		return -1;
	}

	free(tmp);
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

	static const uint8_t sb4_table_local_keyrecovery[16] = {
		0x1, 0x0, 0x2, 0x4, 0x3, 0x8, 0x6, 0xd,
		0x9, 0xa, 0xb, 0xe, 0xf, 0xc, 0x7, 0x5
	};

	static uint8_t sb8_table_local_keyrecovery[256];
	static int sb8_initialized_local_keyrecovery = 0;

	typedef struct {
		int left;
		int right;
	} CommonIndex;

	typedef struct {
		BytePair pairs[MAX_GOOD_PAIRS];
		int count;
	} GoodPairList;

	static void init_sb8_local_keyrecovery(void) {
		if (sb8_initialized_local_keyrecovery) return;
		for (int a = 0; a < 16; a++) {
			for (int b = 0; b < 16; b++) {
				sb8_table_local_keyrecovery[a * 16 + b] = (uint8_t)((sb4_table_local_keyrecovery[a] << 4) | sb4_table_local_keyrecovery[b]);
			}
		}
		sb8_initialized_local_keyrecovery = 1;
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

			if ((sb8_table_local_keyrecovery[x ^ k1] ^ sb8_table_local_keyrecovery[y ^ k1]) == dy1 &&
			    (sb8_table_local_keyrecovery[x ^ k2] ^ sb8_table_local_keyrecovery[y ^ k2]) == dy2) {
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

		init_sb8_local_keyrecovery();

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

static uint32_t getenv_u32(const char *name, uint32_t default_val) {
	const char *env = getenv(name);
	if (!env || strlen(env) == 0) {
		return default_val;
	}
	char *endptr = NULL;
	long val = strtol(env, &endptr, 10);
	if (endptr == env || val <= 0) {
		return default_val;
	}
	return static_cast<uint32_t>(val);
}


static void runExpGoodPairsGPU(unsigned char key[32],
							   unsigned char outMask[32],
							   int numActNib,
							   uint64_t DEGREE,
							   int inpDiffPattern[32],
							   int roundOffset,
							   int num_of_rounds,
							   uint64_t requested_total_threads,
							   const char *fileNameCorr) {
	ensure_device_constants_loaded();

	FILE *fileCorrData = fopen(fileNameCorr, "w");
	if (!fileCorrData) {
		fprintf(stderr, "Failed to open %s for writing\n", fileNameCorr);
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

	const int totalRounds = roundOffset + num_of_rounds + 1;
	if (totalRounds > MAX_TOTAL_ROUNDS) {
		fprintf(stderr, "Error: total rounds %d exceeds MAX_TOTAL_ROUNDS %d\n", totalRounds, MAX_TOTAL_ROUNDS);
		fclose(fileCorrData);
		return;
	}

	unsigned char expdKey1_host[MAX_TOTAL_ROUNDS][32];
	unsigned char expdKey2_host[MAX_TOTAL_ROUNDS][32];
	keyExpansion(key, expdKey1_host, expdKey2_host, totalRounds);

	CUDA_CHECK(cudaMemcpyToSymbol(d_inpDiffPattern, inpDiffPattern, sizeof(int) * 32));
	CUDA_CHECK(cudaMemcpyToSymbol(d_outMask, outMask, sizeof(unsigned char) * 32));
	CUDA_CHECK(cudaMemcpyToSymbol(d_expdKey1, expdKey1_host, sizeof(unsigned char) * totalRounds * 32));
	CUDA_CHECK(cudaMemcpyToSymbol(d_expdKey2, expdKey2_host, sizeof(unsigned char) * totalRounds * 32));

	uint32_t threads_per_block = getenv_u32("CUDA_THREADS_PER_BLOCK", 256);
	uint32_t requested_threads_per_block = threads_per_block;
	uint32_t num_blocks = getenv_u32("CUDA_BLOCKS", 256);

	uint32_t requested_blocks = num_blocks;
	uint64_t desired_total_threads = (uint64_t)requested_threads_per_block * (uint64_t)requested_blocks;
	if (requested_total_threads > 0) {
		desired_total_threads = requested_total_threads;
	}

	int device = 0;
	CUDA_CHECK(cudaGetDevice(&device));
	cudaDeviceProp device_prop;
	CUDA_CHECK(cudaGetDeviceProperties(&device_prop, device));
	cudaFuncAttributes func_attr;
	CUDA_CHECK(cudaFuncGetAttributes(&func_attr, sampling_kernel));

	uint32_t kernel_max_threads = (uint32_t)func_attr.maxThreadsPerBlock;
	uint32_t hardware_max_threads = (uint32_t)device_prop.maxThreadsPerBlock;
	uint32_t absolute_max_threads = kernel_max_threads < hardware_max_threads ? kernel_max_threads : hardware_max_threads;
	if (absolute_max_threads == 0) {
		absolute_max_threads = hardware_max_threads;
	}

	if (threads_per_block > absolute_max_threads) {
		printf("CUDA: clamping threads per block from %u to %u (device limit)\n", threads_per_block, absolute_max_threads);
		threads_per_block = absolute_max_threads;
	}

	int min_grid_size = 0;
	int optimal_block_size = 0;
	CUDA_CHECK(cudaOccupancyMaxPotentialBlockSize(&min_grid_size, &optimal_block_size, sampling_kernel, 0, 0));
	if (optimal_block_size > 0 && threads_per_block > (uint32_t)optimal_block_size) {
		printf("CUDA: clamping threads per block to %d based on occupancy (requested %u)\n", optimal_block_size, threads_per_block);
		threads_per_block = (uint32_t)optimal_block_size;
	}

	if (threads_per_block == 0) {
		threads_per_block = 1;
	}

	if (threads_per_block != requested_threads_per_block || requested_total_threads > 0) {
		uint64_t adjusted_blocks = (desired_total_threads + threads_per_block - 1) / threads_per_block;
		if (adjusted_blocks == 0) {
			adjusted_blocks = 1;
		}
		if (adjusted_blocks > (uint64_t)device_prop.maxGridSize[0]) {
			printf("CUDA: clamping grid dimension from %llu to %d (device limit)\n", (unsigned long long)adjusted_blocks, device_prop.maxGridSize[0]);
			adjusted_blocks = device_prop.maxGridSize[0];
		}
		num_blocks = (uint32_t)adjusted_blocks;
		printf("CUDA: adjusted grid to %u blocks of %u threads (target total %llu)\n",
			num_blocks,
			threads_per_block,
			(unsigned long long)desired_total_threads);
	}

	uint64_t totalThreads = (uint64_t)threads_per_block * (uint64_t)num_blocks;
	if (totalThreads == 0) {
		fprintf(stderr, "Invalid CUDA launch configuration: zero threads\n");
		fclose(fileCorrData);
		return;
	}
	if (requested_total_threads > 0) {
		printf("CUDA: configured %u blocks x %u threads (requested total %llu, actual %llu)\n",
			   num_blocks,
			   threads_per_block,
			   (unsigned long long)requested_total_threads,
			   (unsigned long long)totalThreads);
	}

	curandStatePhilox4_32_10_t *d_rng_states = NULL;
	unsigned char *d_base_plaintext1 = NULL;
	unsigned char *d_base_plaintext2 = NULL;
	unsigned long long *d_counter0 = NULL;
	unsigned long long *d_counter1 = NULL;

	CUDA_CHECK(cudaMalloc(&d_rng_states, totalThreads * sizeof(curandStatePhilox4_32_10_t)));
	CUDA_CHECK(cudaMalloc(&d_base_plaintext1, 32 * sizeof(unsigned char)));
	CUDA_CHECK(cudaMalloc(&d_base_plaintext2, 32 * sizeof(unsigned char)));
	CUDA_CHECK(cudaMalloc(&d_counter0, sizeof(unsigned long long)));
	CUDA_CHECK(cudaMalloc(&d_counter1, sizeof(unsigned long long)));

	printf("\ntot=%zu", g_good_pair_count);
	printf("\nProcessing pairs: ");
	fflush(stdout);

	unsigned char plaintext1[32];
	unsigned char plaintext2[32];

	for (size_t index = 0; index < g_good_pair_count; ++index) {
		// Progress reporting every 10% or every 100 pairs, whichever is more frequent
		if (index % ((g_good_pair_count / 10) + 1) == 0 || index % 100 == 0) {
			printf("Processing pairs: %zu/%zu (%.1f%%)\n", index, g_good_pair_count, 
				   (double)index * 100.0 / g_good_pair_count);
			fflush(stdout);
		}
		for (int i = 0; i < 32; ++i) {
			plaintext1[i] = 0;
			plaintext2[i] = 0;
		}

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

		uint64_t base_seed = (uint64_t)time(NULL) ^ (uint64_t)clock() ^ index;

		init_rng_kernel<<<num_blocks, threads_per_block>>>(d_rng_states,
								  base_seed,
								  totalThreads);
		CUDA_CHECK(cudaGetLastError());
		CUDA_CHECK(cudaDeviceSynchronize());

		CUDA_CHECK(cudaMemcpy(d_base_plaintext1, plaintext1, 32 * sizeof(unsigned char), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemcpy(d_base_plaintext2, plaintext2, 32 * sizeof(unsigned char), cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaMemset(d_counter0, 0, sizeof(unsigned long long)));
		CUDA_CHECK(cudaMemset(d_counter1, 0, sizeof(unsigned long long)));

		sampling_kernel<<<num_blocks, threads_per_block>>>(d_rng_states,
						   d_base_plaintext1,
						   d_base_plaintext2,
						   roundOffset,
						   num_of_rounds,
						   totalThreads,
						   DEGREE,
						   d_counter0,
						   d_counter1);
		CUDA_CHECK(cudaGetLastError());
		CUDA_CHECK(cudaDeviceSynchronize());

		unsigned long long counter0 = 0ULL;
		unsigned long long counter1 = 0ULL;
		CUDA_CHECK(cudaMemcpy(&counter0, d_counter0, sizeof(unsigned long long), cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaMemcpy(&counter1, d_counter1, sizeof(unsigned long long), cudaMemcpyDeviceToHost));

		unsigned long long total_count = counter0 + counter1;
		if (total_count != DEGREE) {
			fprintf(stderr,
				"Warning: GPU sample tally mismatch (got %llu, expected %llu)\n",
				(unsigned long long)total_count,
				(unsigned long long)DEGREE);
		}

		uint64_t CORR = (counter0 > counter1) ? (counter0 - counter1) : (counter1 - counter0);

		char strData[128];
		strData[0] = '\0';
		for (int actNibIndex = 0; actNibIndex < numActNib && actNibIndex < g_num_active_nibbles; ++actNibIndex) {
			char temp[16];
			uint8_t px = (actNibIndex == 0) ? nib_high(pair.x) : nib_low(pair.x);
			uint8_t py = (actNibIndex == 0) ? nib_high(pair.y) : nib_low(pair.y);
			snprintf(temp, sizeof(temp), "%d %d ", px, py);
			strncat(strData, temp, sizeof(strData) - strlen(strData) - 1);
		}
		char corr_str[32];
		snprintf(corr_str, sizeof(corr_str), "%llu", (unsigned long long)CORR);
		strncat(strData, corr_str, sizeof(strData) - strlen(strData) - 1);

		write_string_to_file(fileCorrData, strData);
	}

	printf("\nComplete!\n");
	fflush(stdout);

	fclose(fileCorrData);

	cudaFree(d_rng_states);
	cudaFree(d_base_plaintext1);
	cudaFree(d_base_plaintext2);
	cudaFree(d_counter0);
	cudaFree(d_counter1);
}

int main(int argc, char *argv[]) {
	if (argc < 14 || argc > 15) {
		fprintf(stderr,
				"Usage: %s num_of_rounds rndOffset outputmask inputDiffPattern "
				"keyValType plainNibValType DEG N trailFolder dy1_hi dy1_lo dy2_hi dy2_lo [num_threads]\n",
				argv[0]);
		fprintf(stderr, "  num_threads: optional, defaults to auto-detection (use all available GPU threads).\n");
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
	int plainNibValType = atoi(argv[6]);
	(void)plainNibValType;
	int DEG = atoi(argv[7]);
	int N = atoi(argv[8]);

	char trailFolder[400];
	strncpy(trailFolder, argv[9], sizeof(trailFolder) - 1);
	trailFolder[sizeof(trailFolder) - 1] = '\0';

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

	uint64_t requested_total_threads = 0;
	if (argc == 15) {
		long tmp = atol(argv[14]);
		if (tmp < 0) {
			fprintf(stderr, "Error: num_threads must be positive (or 0 for auto)\n");
			return 1;
		}
		requested_total_threads = (uint64_t)tmp;
	}

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

	char dir_name[512];
	snprintf(dir_name, sizeof(dir_name), "data/%s/output%d", trailFolder, DEG);
	if (mkdir_p(dir_name) == -1) {
		perror("Error creating directory");
		kr_release_enumeration_data();
		free_key_analysis_data(&analysis_data);
		return 1;
	}

	int activeNib = num_active;
	int (*keyBitPos)[8] = (int (*)[8])malloc(activeNib * sizeof(*keyBitPos));
	if (!keyBitPos) {
		fprintf(stderr, "Error: memory allocation failed for keyBitPos\n");
		kr_release_enumeration_data();
		free_key_analysis_data(&analysis_data);
		return 1;
	}

	keyNibs = (unsigned char (*)[2])malloc(activeNib * sizeof(*keyNibs));
	if (!keyNibs) {
		fprintf(stderr, "Error: memory allocation failed for keyNibs\n");
		free(keyBitPos);
		kr_release_enumeration_data();
		free_key_analysis_data(&analysis_data);
		return 1;
	}

	uint64_t N1 = 1ULL << DEG;
	srand(time(NULL));
	(void)init_prng(((unsigned int)rand() << 16) | rand());
	getKeyBitPosNew(inputDiffPattern, keyBitPos, activeNib, rndOffset);

	unsigned char outMask[32];
	if (hex_string_to_bytes(outputmask, outMask, sizeof(outMask)) != 0) {
		fprintf(stderr, "Invalid output mask\n");
		free(keyNibs);
		free(keyBitPos);
		kr_release_enumeration_data();
		free_key_analysis_data(&analysis_data);
		return 1;
	}

	char file_name_corr[768];
	size_t dir_len = strlen(dir_name);
	size_t suffix_allowance = (size_t)activeNib * 6U + 4U;
	if (dir_len + strlen("/out") + suffix_allowance + 1U > sizeof(file_name_corr)) {
		fprintf(stderr, "Error: output path is too long\n");
		free(keyNibs);
		free(keyBitPos);
		kr_release_enumeration_data();
		free_key_analysis_data(&analysis_data);
		return 1;
	}

	for (int num_of_experiments = 0; num_of_experiments < N; ++num_of_experiments) {
		printf("\n----------------------------------------------------------------------------------------------------");
		printf("\n Exp No: %d\n", num_of_experiments);
		int written = snprintf(file_name_corr, sizeof(file_name_corr), "%s/out", dir_name);
		if (written < 0 || (size_t)written >= sizeof(file_name_corr)) {
			fprintf(stderr, "Error: failed to build output filename\n");
			free(keyNibs);
			free(keyBitPos);
			kr_release_enumeration_data();
			free_key_analysis_data(&analysis_data);
			return 1;
		}
		(void)init_prng(233);
		unsigned char key[32];
		kr_set_current_mode(keyValType);
		fixKey(key, keyBitPos, keyValType, activeNib, file_name_corr);
		runExpGoodPairsGPU(key,
					outMask,
					activeNib,
					N1,
					inputDiffPattern,
					rndOffset,
					num_of_rounds,
					requested_total_threads,
					file_name_corr);
	}

	kr_release_enumeration_data();
	free_key_analysis_data(&analysis_data);
	free(keyNibs);
	free(keyBitPos);

	clock_t end_time = clock();
	double runtime = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	printf("\nRuntime: %f seconds\n", runtime);
	return 0;
}
