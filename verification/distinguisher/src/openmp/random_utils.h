#pragma once

#include <cstdint>

struct RandomStats {
    uint64_t rdrand_calls = 0;
    uint64_t getrandom_calls = 0;
    uint64_t fallback_calls = 0;
    uint64_t total_calls = 0;
};

void reset_random_usage();
RandomStats current_random_stats();

bool get_random_uint64(uint64_t &value);
inline bool get_random_uint64(uint64_t *value) {
    if (!value) {
        return false;
    }
    return get_random_uint64(*value);
}

uint64_t splitmix64(uint64_t &state);
uint64_t derive_thread_seed(uint64_t global_seed,
                            uint64_t logical_thread,
                            int omp_thread_id,
                            uint64_t experiment_id);

void report_random_source();
void report_random_usage_stats();
