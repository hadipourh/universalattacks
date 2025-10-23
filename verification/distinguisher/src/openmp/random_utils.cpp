#include "openmp/random_utils.h"

#include <cstdio>
#include <cstdint>
#include <ctime>
#include <inttypes.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __linux__
#include <sys/random.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <Security/SecRandom.h>
#endif

#if (defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)) && defined(__RDRND__)
#include <x86intrin.h>
#include <cpuid.h>
#define HAVE_X86_RDRAND 1
#else
#define HAVE_X86_RDRAND 0
#endif

namespace {
RandomStats g_random_stats;

uint64_t mix_entropy(uint64_t base, uint64_t salt) {
    base ^= 0x9E3779B97F4A7C15ULL * (salt + 0x9E3779B97F4A7C15ULL);
    base = (base ^ (base >> 30)) * 0xBF58476D1CE4E5B9ULL;
    base = (base ^ (base >> 27)) * 0x94D049BB133111EBULL;
    return base ^ (base >> 31);
}

#if HAVE_X86_RDRAND
int check_rdrand_support() {
    uint32_t eax, ebx, ecx, edx;
    __cpuid(0, eax, ebx, ecx, edx);
    if (ebx != 0x756e6547 || edx != 0x49656e69 || ecx != 0x6c65746e) {
        return 0;
    }
    __cpuid(1, eax, ebx, ecx, edx);
    return (ecx & 0x40000000) ? 1 : 0;
}
#endif

int rdrand_64_impl(uint64_t &x) {
#if HAVE_X86_RDRAND
    static int supported = -1;
    if (supported < 0) {
        supported = check_rdrand_support();
    }
    if (supported) {
        for (int i = 0; i < 10; ++i) {
            unsigned long long tmp;
            if (_rdrand64_step(&tmp)) {
                x = static_cast<uint64_t>(tmp);
                g_random_stats.total_calls++;
                g_random_stats.rdrand_calls++;
                return 1;
            }
        }
    }
#endif
#ifdef __linux__
    if (getrandom(&x, sizeof(x), 0) == sizeof(x)) {
        g_random_stats.total_calls++;
        g_random_stats.getrandom_calls++;
        return 1;
    }
#endif
#ifdef __APPLE__
    if (SecRandomCopyBytes(kSecRandomDefault, sizeof(x), &x) == errSecSuccess) {
        g_random_stats.total_calls++;
        g_random_stats.getrandom_calls++;
        return 1;
    }
#endif
    const uint64_t time_mix = (static_cast<uint64_t>(time(nullptr)) << 32) ^ static_cast<uint64_t>(clock());
    x = time_mix ^ 0xA5A5A5A5A5A5A5A5ULL;
    g_random_stats.total_calls++;
    g_random_stats.fallback_calls++;
    return 1;
}

} // namespace

void reset_random_usage() {
    g_random_stats = {};
}

RandomStats current_random_stats() {
    return g_random_stats;
}

bool get_random_uint64(uint64_t &value) {
    return rdrand_64_impl(value) == 1;
}

uint64_t splitmix64(uint64_t &state) {
    state += 0x9E3779B97F4A7C15ULL;
    uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

uint64_t derive_thread_seed(uint64_t global_seed,
                            uint64_t logical_thread,
                            int omp_thread_id,
                            uint64_t experiment_id) {
    uint64_t seed = global_seed ^ mix_entropy(global_seed, logical_thread + 0x1000U);
    seed ^= mix_entropy(seed, experiment_id + 0xB5AD4ECEDA1CE2A9ULL);
    seed ^= mix_entropy(seed, static_cast<uint64_t>(omp_thread_id) + 0xD1B54A32D192ED03ULL);
    seed ^= (logical_thread + 1) * 0x94D049BB133111EBULL;
    return seed;
}

void report_random_source() {
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
#elif defined(__linux__)
    printf("Random source: System entropy (software)\n");
    printf("  - Primary: getrandom() system call\n");
    printf("  - Fallback: time-based mixing\n");
#elif defined(__APPLE__)
    printf("Random source: macOS Security framework\n");
    printf("  - Primary: SecRandomCopyBytes()\n");
    printf("  - Fallback: time-based mixing\n");
#else
    printf("Random source: Time-based mixing (software fallback)\n");
#endif
    printf("  - Status: Testing actual usage during execution...\n\n");
}

void report_random_usage_stats() {
    if (g_random_stats.total_calls == 0) {
        return;
    }
    printf("Random Source Usage Statistics:\n");
    printf("================================\n");
    printf("Total random calls: %" PRIu64 "\n", g_random_stats.total_calls);
    if (g_random_stats.rdrand_calls > 0) {
        printf("RDRAND (hardware): %" PRIu64 " calls (%.1f%%)\n",
               g_random_stats.rdrand_calls,
               100.0 * g_random_stats.rdrand_calls / g_random_stats.total_calls);
    }
    if (g_random_stats.getrandom_calls > 0) {
        printf("getrandom() (system): %" PRIu64 " calls (%.1f%%)\n",
               g_random_stats.getrandom_calls,
               100.0 * g_random_stats.getrandom_calls / g_random_stats.total_calls);
    }
    if (g_random_stats.fallback_calls > 0) {
        printf("Time-based fallback: %" PRIu64 " calls (%.1f%%)\n",
               g_random_stats.fallback_calls,
               100.0 * g_random_stats.fallback_calls / g_random_stats.total_calls);
    }
    if (g_random_stats.rdrand_calls > 0) {
        printf("\u2713 Verified: RDRAND hardware entropy is working\n");
    } else if (g_random_stats.getrandom_calls > 0) {
        printf("\u2713 Verified: getrandom() system entropy is working\n");
    } else if (g_random_stats.fallback_calls > 0) {
        printf("\u26A0 Warning: Only time-based fallback was used\n");
    }
    printf("\n");
}
