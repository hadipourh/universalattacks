#include "common_types.h"
#include "openmp/random_utils.h"
#include "orthros_core.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <inttypes.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

extern void display_correlation_histogram(const double *correlations, int num_experiments, Config *config);
extern void display_experiment_keys(ExperimentKey *keys, int num_experiments, Config *config);

namespace {

constexpr int DEFAULT_BATCH_SIZE = 64;
constexpr int MAX_EXPERIMENTS = 1000;
constexpr int DLCT_PRECISION_SCALE = 100;

struct ExperimentWorkspace {
	std::vector<uint64_t> thread_states;
	std::vector<ExperimentKey> experiment_keys;
	std::vector<double> correlation_values;
	bool seeds_initialized = false;

	void ensure(int total_threads, int experiments) {
		if (static_cast<int>(thread_states.size()) != total_threads) {
			thread_states.assign(total_threads, 0ULL);
			seeds_initialized = false;
		}
		if (static_cast<int>(experiment_keys.size()) != experiments) {
			experiment_keys.resize(experiments);
		}
		if (static_cast<int>(correlation_values.size()) != experiments) {
			correlation_values.assign(experiments, 0.0);
		} else {
			std::fill(correlation_values.begin(), correlation_values.begin() + experiments, 0.0);
		}
	}
};

ExperimentWorkspace &workspace() {
	static ExperimentWorkspace instance;
	return instance;
}

uint64_t mix_time_seed(int experiment_id, int salt) {
	const uint64_t t = static_cast<uint64_t>(time(nullptr));
	const uint64_t c = static_cast<uint64_t>(clock());
	return (t << 32) ^ c ^ (static_cast<uint64_t>(experiment_id + 1) * 0x9E3779B97F4A7C15ULL) ^
		   (static_cast<uint64_t>(salt + 37) * 0xBF58476D1CE4E5B9ULL);
}

void initialize_thread_states(ExperimentWorkspace &ws, int total_threads) {
	uint64_t base_seed = 0;
	if (!get_random_uint64(base_seed)) {
		base_seed = mix_time_seed(0, 0);
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
	for (int tid = 0; tid < total_threads; ++tid) {
		int host_tid = 0;
#ifdef _OPENMP
		host_tid = omp_get_thread_num();
#endif
		ws.thread_states[tid] = derive_thread_seed(base_seed, static_cast<uint64_t>(tid), host_tid, 0ULL);
	}

	ws.seeds_initialized = true;
}

struct ExperimentTimer {
	void start() {
		start_point = std::chrono::steady_clock::now();
	}

	double elapsed_ms() const {
		const auto now = std::chrono::steady_clock::now();
		return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(now - start_point).count();
	}

	std::chrono::steady_clock::time_point start_point;
};

void generate_random_key(unsigned char *key, int experiment_id) {
	int nib = 0;
	while (nib < ORTHROS_STATE_SIZE) {
		uint64_t value;
		if (!get_random_uint64(value)) {
			value = mix_time_seed(experiment_id, nib);
		}
		for (int j = 0; j < 16 && nib < ORTHROS_STATE_SIZE; ++j) {
			key[nib++] = static_cast<unsigned char>((value >> (j * 4)) & 0xF);
		}
	}
}

void assign_experiment_key(ExperimentWorkspace &ws, int experiment_id) {
	ExperimentKey &entry = ws.experiment_keys[experiment_id];
	entry.experiment_id = experiment_id;
	generate_random_key(entry.key, experiment_id);
	entry.correlation = 0.0;
}

void fill_random_state(uint64_t &seed, unsigned char *state) {
	const uint64_t lo = splitmix64(seed);
	const uint64_t hi = splitmix64(seed);
	for (int i = 0; i < 16; ++i) {
		state[i] = static_cast<unsigned char>((lo >> (i * 4)) & 0xF);
		state[16 + i] = static_cast<unsigned char>((hi >> (i * 4)) & 0xF);
	}
}

void xor_with_diff(const unsigned char *base, const unsigned char *diff, unsigned char *out) {
	for (int i = 0; i < ORTHROS_STATE_SIZE; ++i) {
		out[i] = static_cast<unsigned char>((base[i] ^ diff[i]) & 0xF);
	}
}

double process_single_experiment(const Config &config,
								 ExperimentWorkspace &ws,
								 int experiment_id,
								 int batch_relative_id,
								 uint64_t samples_per_thread,
								 int total_threads) {
	unsigned char key_local[ORTHROS_STATE_SIZE];
	std::memcpy(key_local, ws.experiment_keys[experiment_id].key, ORTHROS_STATE_SIZE);

	ExperimentTimer timer;
	timer.start();

	uint64_t total_equal = 0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_equal) schedule(static)
#endif
	for (int tid = 0; tid < total_threads; ++tid) {
		uint64_t seed = ws.thread_states[tid];
#ifdef _OPENMP
		seed ^= 0xD1B54A32D192ED03ULL * static_cast<uint64_t>(omp_get_thread_num() + 1);
#endif

		unsigned char base_left[ORTHROS_STATE_SIZE];
		unsigned char base_right[ORTHROS_STATE_SIZE];
		unsigned char diff_left[ORTHROS_STATE_SIZE];
		unsigned char diff_right[ORTHROS_STATE_SIZE];
		unsigned char ct_base[ORTHROS_STATE_SIZE];
		unsigned char ct_diff[ORTHROS_STATE_SIZE];

		uint64_t local_equal = 0;

		for (uint64_t sample = 0; sample < samples_per_thread; ++sample) {
			fill_random_state(seed, base_left);
			std::memcpy(base_right, base_left, ORTHROS_STATE_SIZE);
			xor_with_diff(base_left, config.input_diff_left, diff_left);
			xor_with_diff(base_left, config.input_diff_right, diff_right);

			orthros_encrypt(base_left, base_right, ct_base, key_local,
							config.offset, config.rounds, config.mode);
			orthros_encrypt(diff_left, diff_right, ct_diff, key_local,
							config.offset, config.rounds, config.mode);

			const int parity_base = dot_product(ct_base, config.output_mask);
			const int parity_diff = dot_product(ct_diff, config.output_mask);
			if (parity_base == parity_diff) {
				local_equal++;
			}
		}

		ws.thread_states[tid] = seed;
		total_equal += local_equal;
	}

	const double elapsed_ms = timer.elapsed_ms();
	const uint64_t total_samples = static_cast<uint64_t>(total_threads) * samples_per_thread;
	const int64_t signed_diff = static_cast<int64_t>(total_equal) * 2 - static_cast<int64_t>(total_samples);
	const double correlation = static_cast<double>(signed_diff) / static_cast<double>(total_samples);

	if (config.verbose) {
		const double corr_mag = std::fabs(correlation);
		const char sign_char = (correlation >= 0.0) ? '+' : '-';
		if (corr_mag > 0.0) {
			const double log_corr = std::log2(corr_mag);
			std::printf("Experiment %d: Equal pairs: %" PRIu64 "/%" PRIu64 ", Correlation: %c2^%.2f, Time: %.2f ms\n",
						experiment_id + 1, total_equal, total_samples, sign_char, log_corr, elapsed_ms);
		} else {
			std::printf("Experiment %d: Equal pairs: %" PRIu64 "/%" PRIu64 ", Correlation: %c2^(-inf), Time: %.2f ms\n",
						experiment_id + 1, total_equal, total_samples, sign_char, elapsed_ms);
		}
	}

	ws.experiment_keys[experiment_id].correlation = correlation;
	return correlation;
}

void zero_state(unsigned char *state) {
	std::memset(state, 0, ORTHROS_STATE_SIZE);
}

void set_single_bit(unsigned char *state, int bit_index) {
	zero_state(state);
	const int nib_index = bit_index / 4;
	const int bit_pos = 3 - (bit_index % 4);
	state[nib_index] = static_cast<unsigned char>(1u << bit_pos);
}

void format_state(const unsigned char *state, char *out) {
	static const char HEX_DIGITS[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};
	for (int i = 0; i < ORTHROS_STATE_SIZE; ++i) {
		out[i] = HEX_DIGITS[state[i] & 0xF];
	}
	out[ORTHROS_STATE_SIZE] = '\0';
}

void validate_config(Config &config) {
	if (config.experiments > MAX_EXPERIMENTS) {
		std::printf("Warning: Limiting experiments from %d to %d\n", config.experiments, MAX_EXPERIMENTS);
		config.experiments = MAX_EXPERIMENTS;
	}
}

} // namespace

double run_correlation_suite(Config config) {
	validate_config(config);

	const int total_threads = config.blocks * config.threads;
	const uint64_t samples_per_thread = 1ULL << config.sample_power;
	const uint64_t total_samples = static_cast<uint64_t>(total_threads) * samples_per_thread;

	if (config.verbose) {
	char left_hex[ORTHROS_STATE_SIZE + 1];
	char right_hex[ORTHROS_STATE_SIZE + 1];
	char mask_hex[ORTHROS_STATE_SIZE + 1];
	format_state(config.input_diff_left, left_hex);
	format_state(config.input_diff_right, right_hex);
	format_state(config.output_mask, mask_hex);

		std::printf("Computing correlation with:\n");
		std::printf("  Rounds: %d (offset %d)\n", config.rounds, config.offset);
		std::printf("  Mode: %d\n", static_cast<int>(config.mode));
		std::printf("  ΔL: %s\n", left_hex);
		std::printf("  ΔR: %s\n", right_hex);
		std::printf("  Output mask: %s\n", mask_hex);
		std::printf("  Grid: %d blocks × %d threads = %d total threads\n",
					config.blocks, config.threads, total_threads);
		std::printf("  Samples per thread: 2^%d = %" PRIu64 "\n", config.sample_power, samples_per_thread);
		std::printf("  Total samples per experiment: %" PRIu64 "\n", total_samples);
		std::printf("  Experiments: %d\n\n", config.experiments);
#ifdef _OPENMP
		std::printf("  OpenMP threads available: %d\n\n", omp_get_max_threads());
#endif
	}

	ExperimentWorkspace &ws = workspace();
	ws.ensure(total_threads, config.experiments);
	if (!ws.seeds_initialized) {
		initialize_thread_states(ws, total_threads);
	}

	double sum_correlations = 0.0;
	const int batch_count = (config.experiments + DEFAULT_BATCH_SIZE - 1) / DEFAULT_BATCH_SIZE;

	for (int batch = 0; batch < batch_count; ++batch) {
		const int exp_start = batch * DEFAULT_BATCH_SIZE;
		const int exp_end = std::min(exp_start + DEFAULT_BATCH_SIZE, config.experiments);
		const int experiments_in_batch = exp_end - exp_start;

		for (int i = 0; i < experiments_in_batch; ++i) {
			assign_experiment_key(ws, exp_start + i);
		}

		if (config.verbose && batch == 0) {
			std::printf("Processing %d batches of up to %d experiments each\n\n", batch_count, DEFAULT_BATCH_SIZE);
		}

		for (int exp = exp_start; exp < exp_end; ++exp) {
			const double corr = process_single_experiment(config,
														  ws,
														  exp,
														  exp - exp_start,
														  samples_per_thread,
														  total_threads);
			ws.correlation_values[exp] = corr;
			sum_correlations += corr;
		}
	}

	if (config.experiments > 1) {
		display_correlation_histogram(ws.correlation_values.data(), config.experiments, &config);
	}
	display_experiment_keys(ws.experiment_keys.data(), config.experiments, &config);

	return sum_correlations / config.experiments;
}

void generate_dlct_table_cpu(Config config) {
	constexpr int BIT_COUNT = ORTHROS_STATE_SIZE * 4;

	std::printf("Generating DLCT table for %d-round Orthros (offset %d, mode %d)\n",
				config.rounds, config.offset, static_cast<int>(config.mode));
	std::printf("Configuration: %d blocks × %d threads, 2^%d samples, %d experiments\n\n",
				config.blocks, config.threads, config.sample_power, config.experiments);

	char csv_filename[256];
	char dzn_filename[256];
	std::snprintf(csv_filename, sizeof(csv_filename), "orthros_%dr_dlct_cpu.csv", config.rounds);
	std::snprintf(dzn_filename, sizeof(dzn_filename), "orthros_%dr_dlct_cpu.dzn", config.rounds);

	FILE *csv_file = std::fopen(csv_filename, "w");
	FILE *dzn_file = std::fopen(dzn_filename, "w");

	if (!csv_file || !dzn_file) {
		std::fprintf(stderr, "Error: Could not create output files\n");
		if (csv_file) std::fclose(csv_file);
		if (dzn_file) std::fclose(dzn_file);
		return;
	}

	std::fprintf(csv_file, "Input_Bit,Output_Bit,Diff_Left_Hex,Diff_Right_Hex,Output_Hex,Correlation_Pow2,Log2_Magnitude\n");
	std::fprintf(dzn_file, "%% Orthros %d-round DLCT (single-bit patterns, CPU OpenMP)\n", config.rounds);
	const unsigned long long total_samples = static_cast<unsigned long long>(config.blocks) * config.threads * (1ULL << config.sample_power);
	std::fprintf(dzn_file,
				 "%% rounds=%d, offset=%d, mode=%d, blocks=%d, threads=%d, samples_per_thread=2^%d, total_samples=%llu, experiments=%d\n",
				 config.rounds, config.offset, static_cast<int>(config.mode),
				 config.blocks, config.threads, config.sample_power, total_samples, config.experiments);

	ExperimentWorkspace &ws = workspace();
	ws.seeds_initialized = false;

	int dlct_weights[BIT_COUNT][BIT_COUNT];
	int max_weight = 0;

	std::printf("Progress: ");
	std::fflush(stdout);

	unsigned char single_left[ORTHROS_STATE_SIZE];
	unsigned char single_right[ORTHROS_STATE_SIZE];
	unsigned char single_mask[ORTHROS_STATE_SIZE];

	for (int input_bit = 0; input_bit < BIT_COUNT; ++input_bit) {
		if (input_bit % 16 == 0) {
			std::printf(".");
			std::fflush(stdout);
		}

		set_single_bit(single_left, input_bit);
		zero_state(single_right);

		for (int output_bit = 0; output_bit < BIT_COUNT; ++output_bit) {
			set_single_bit(single_mask, output_bit);

			Config temp = config;
			std::memcpy(temp.input_diff_left, single_left, ORTHROS_STATE_SIZE);
			std::memcpy(temp.input_diff_right, single_right, ORTHROS_STATE_SIZE);
			std::memcpy(temp.output_mask, single_mask, ORTHROS_STATE_SIZE);
			temp.verbose = false;

			const double correlation = run_correlation_suite(temp);
			const double corr_mag = std::fabs(correlation);

			char left_hex[ORTHROS_STATE_SIZE + 1];
			char right_hex[ORTHROS_STATE_SIZE + 1];
			char mask_hex[ORTHROS_STATE_SIZE + 1];
			format_state(single_left, left_hex);
			format_state(single_right, right_hex);
			format_state(single_mask, mask_hex);

			int weight;
			if (corr_mag == 0.0) {
				std::fprintf(csv_file, "%d,%d,%s,%s,%s,+2^(-inf),inf\n",
							 input_bit, output_bit, left_hex, right_hex, mask_hex);
				weight = -1;
			} else {
				const double exponent = std::log2(corr_mag);
				const char sign_char = (correlation >= 0.0) ? '+' : '-';
				std::fprintf(csv_file, "%d,%d,%s,%s,%s,%c2^(%.2f),%.6f\n",
							 input_bit, output_bit, left_hex, right_hex, mask_hex,
							 sign_char, exponent, -exponent);
				const double scaled = DLCT_PRECISION_SCALE * (-exponent);
				weight = static_cast<int>(std::llround(scaled));
				if (weight > max_weight) {
					max_weight = weight;
				}
			}
			dlct_weights[input_bit][output_bit] = weight;
		}
	}

	std::printf(" Done!\n\n");

	const int sentinel_weight = max_weight + 100;
	std::fprintf(dzn_file, "array[1..%d,1..%d] of int: DLCT =\n[|\n", BIT_COUNT, BIT_COUNT);
	for (int i = 0; i < BIT_COUNT; ++i) {
		if (i > 0) std::fprintf(dzn_file, " |\n");
		for (int j = 0; j < BIT_COUNT; ++j) {
			const int value = (dlct_weights[i][j] < 0) ? sentinel_weight : dlct_weights[i][j];
			std::fprintf(dzn_file, "%d", value);
			if (j < BIT_COUNT - 1) std::fprintf(dzn_file, ", ");
		}
	}
	std::fprintf(dzn_file, " |\n];\n\n");
	std::fprintf(dzn_file, "int: DLCT_PRECISION_SCALE = %d;\n", DLCT_PRECISION_SCALE);
	std::fprintf(dzn_file, "int: DLCT_SENTINEL_WEIGHT = %d;\n", sentinel_weight);

	std::fclose(csv_file);
	std::fclose(dzn_file);

	std::printf("DLCT table generation completed!\n");
	std::printf("Files generated:\n");
	std::printf("  - %s (CSV format for analysis)\n", csv_filename);
	std::printf("  - %s (MiniZinc format for optimization)\n", dzn_filename);
}

