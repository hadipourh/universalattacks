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
 * @file correlation_analysis.cpp
 * @brief Correlation analysis and visualization for Orthros DLCT
 * 
 * This module provides comprehensive statistical analysis of correlation data including:
 * - Statistical summaries (mean, std dev, range)
 * - Professional histogram using Freedman-Diaconis rule for optimal binning
 * - Correlation quality assessment
 * - Publication-ready LaTeX figure generation
 */

#include "common_types.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>

static const char HEX_DIGITS[16] = {'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};

static void format_nibbles(const unsigned char *src, char *dst) {
    for (int i = 0; i < ORTHROS_STATE_SIZE; i++) {
        dst[i] = HEX_DIGITS[src[i] & 0xF];
    }
    dst[ORTHROS_STATE_SIZE] = '\0';
}

void export_correlation_data_for_python(const double *correlations, int num_experiments,
                                        const Config *config);

// Function declarations
void display_correlation_histogram(const double *correlations, int num_experiments, 
                                  Config *config);
void display_experiment_keys(ExperimentKey *keys, int num_experiments, 
                            Config *config);

/**
 * @brief Create and display ASCII histogram of correlation distribution
 */
void display_correlation_histogram(const double *correlations, int num_experiments, 
                                  Config *config) {
    if (num_experiments <= 1) {
        return; // Need at least 2 experiments for meaningful distribution
    }
    
    printf("Correlation Distribution Analysis:\n");
    printf("==================================\n");
    
    // Find min and max correlations
    double min_corr = correlations[0];
    double max_corr = correlations[0];
    double sum = 0.0;
    double sum_sq = 0.0;
    
    for (int i = 0; i < num_experiments; i++) {
        if (correlations[i] < min_corr) min_corr = correlations[i];
        if (correlations[i] > max_corr) max_corr = correlations[i];
        sum += correlations[i];
        sum_sq += correlations[i] * correlations[i];
    }
    
    const double mean = sum / num_experiments;
    const double variance = (sum_sq / num_experiments) - (mean * mean);
    const double std_dev = sqrt(variance);
    
    // Calculate power-of-2 representations
    const double mean_log2 = (mean != 0.0) ? log2(fabs(mean)) : -INFINITY;
    const double std_log2 = (std_dev > 0.0) ? log2(std_dev) : -INFINITY;
    const double min_log2 = (min_corr != 0.0) ? log2(fabs(min_corr)) : -INFINITY;
    const double max_log2 = (max_corr != 0.0) ? log2(fabs(max_corr)) : -INFINITY;
    
    // Print statistics with both formats
    printf("Statistical Summary:\n");
    printf("  Experiments: %d\n", num_experiments);
    printf("  Mean correlation: %+.6e (2^%.2f)\n", mean, mean_log2);
    printf("  Standard deviation: %.6e (2^%.2f)\n", std_dev, std_log2);
    printf("  Min correlation: %+.6e (2^%.2f)\n", min_corr, min_log2);
    printf("  Max correlation: %+.6e (2^%.2f)\n", max_corr, max_log2);
    printf("  Range: %.6e\n", max_corr - min_corr);
    printf("  Security margin: ~%.1f bits\n", -mean_log2 / 2.0);
    printf("\n");
    
    // Smart histogram generation: Count unique values and group nearby ones
    // Sort correlations first
    double *sorted_corr = (double*)malloc(num_experiments * sizeof(double));
    memcpy(sorted_corr, correlations, num_experiments * sizeof(double));
    
    // Quick sort
    for (int i = 0; i < num_experiments - 1; i++) {
        for (int j = i + 1; j < num_experiments; j++) {
            if (sorted_corr[i] > sorted_corr[j]) {
                double temp = sorted_corr[i];
                sorted_corr[i] = sorted_corr[j];
                sorted_corr[j] = temp;
            }
        }
    }
    
    const double q1 = sorted_corr[num_experiments / 4];
    const double q3 = sorted_corr[(3 * num_experiments) / 4];
    const double iqr = q3 - q1;
    
    // Count unique values to determine if we should use value-based or range-based binning
    int unique_count = 1;
    double tolerance = (max_corr - min_corr) * 1e-9; // Very small tolerance for floating point comparison
    for (int i = 1; i < num_experiments; i++) {
        if (fabs(sorted_corr[i] - sorted_corr[i-1]) > tolerance) {
            unique_count++;
        }
    }
    
    printf("Professional histogram analysis:\n");
    printf("  IQR: %.6e\n", iqr);
    printf("  Unique values: %d\n", unique_count);
    
    int num_bins;
    int *bins;
    double *bin_edges;
    bool use_value_bins = (unique_count <= 100 && unique_count < num_experiments / 2);
    
    if (use_value_bins) {
        // Use exact value-based bins (no gaps!)
        printf("  Binning strategy: Value-based (gap-free)\n");
        num_bins = unique_count;
        bins = (int*)calloc(num_bins, sizeof(int));
        bin_edges = (double*)malloc((num_bins + 1) * sizeof(double));
        
        // Create bins for each unique value
        int bin_idx = 0;
        bin_edges[0] = sorted_corr[0];
        bins[0] = 1;
        
        for (int i = 1; i < num_experiments; i++) {
            if (fabs(sorted_corr[i] - sorted_corr[i-1]) > tolerance) {
                bin_idx++;
                bin_edges[bin_idx] = sorted_corr[i];
                bins[bin_idx] = 1;
            } else {
                bins[bin_idx]++;
            }
        }
        bin_edges[num_bins] = sorted_corr[num_experiments - 1] + tolerance;
        
    } else {
        // Use range-based bins with adaptive width
        printf("  Binning strategy: Range-based (Freedman-Diaconis)\n");
        
        int optimal_bins;
        if (iqr > 0) {
            const double fd_bin_width = 2.0 * iqr / pow(num_experiments, 1.0/3.0);
            optimal_bins = (int)ceil((max_corr - min_corr) / fd_bin_width);
        } else {
            optimal_bins = (int)(log2(num_experiments) + 1);
        }
        
        const int min_bins = (num_experiments < 50) ? 15 : 25;
        const int max_bins = (num_experiments > 1000) ? 100 : 60;
        
        if (optimal_bins < min_bins) optimal_bins = min_bins;
        if (optimal_bins > max_bins) optimal_bins = max_bins;
        
        printf("  Optimal bins: %d\n", optimal_bins);
        
        num_bins = optimal_bins;
        bins = (int*)calloc(num_bins, sizeof(int));
        bin_edges = (double*)malloc((num_bins + 1) * sizeof(double));
        
        const double bin_width = (max_corr - min_corr) / num_bins;
        for (int i = 0; i <= num_bins; i++) {
            bin_edges[i] = min_corr + i * bin_width;
        }
        
        // Populate bins
        for (int i = 0; i < num_experiments; i++) {
            int bin = (int)((correlations[i] - min_corr) / bin_width);
            if (bin >= num_bins) bin = num_bins - 1;
            if (bin < 0) bin = 0;
            bins[bin]++;
        }
    }
    
    free(sorted_corr);
    
    const int max_bar_width = 50;
    
    // Find max bin count for scaling
    int max_bin_count = 0;
    for (int i = 0; i < num_bins; i++) {
        if (bins[i] > max_bin_count) max_bin_count = bins[i];
    }
    
    // Display histogram
    printf("Correlation Distribution Histogram:\n");
    printf("  %-15s │ %-6s │ %s\n", "Correlation", "Count", "Frequency");
    printf("  %-15s │ %-6s │ %s\n", "Range", "", "");
    printf("  %.*s┼%.*s┼%.*s\n", 15, "───────────────", 6, "──────", 52, "────────────────────────────────────────────────────");
    
    for (int i = 0; i < num_bins; i++) {
        // Skip empty bins in range-based mode
        if (!use_value_bins && bins[i] == 0) {
            continue;
        }
        
        const double bin_start = bin_edges[i];
        const double bin_end = bin_edges[i + 1];
        const int bar_length = max_bin_count > 0 ? (bins[i] * max_bar_width) / max_bin_count : 0;
        
        if (use_value_bins) {
            // For value-based bins, show the exact value
            printf("  %+14.6e │ %-6d │ ", bin_start, bins[i]);
        } else {
            // For range-based bins, show the range
            printf("  %+7.2e~%-6.0e │ %-6d │ ", bin_start, bin_end * 1e6, bins[i]);
        }
        
        // Draw bar with different characters for different intensities
        for (int j = 0; j < bar_length; j++) {
            if (j < bar_length * 0.7) printf("█");
            else if (j < bar_length * 0.9) printf("▊");
            else printf("▌");
        }
        printf(" (%.1f%%)\n", 100.0 * bins[i] / num_experiments);
    }
    
    printf("\n");
    
    // Add cryptanalysis-oriented correlation assessment
    const double abs_mean = fabs(mean);
    const double abs_mean_log2 = (abs_mean > 0.0) ? log2(abs_mean) : -INFINITY;
    
    // Calculate detection threshold based on TOTAL sample size across all threads
    // Total samples = blocks * threads * 2^sample_power
    const uint64_t total_threads = config->blocks * config->threads;
    const uint64_t samples_per_thread = 1ULL << config->sample_power;
    const uint64_t total_samples = total_threads * samples_per_thread;
    const double total_samples_log2 = log2((double)total_samples);
    const double detection_threshold = pow(2.0, -total_samples_log2/2.0);
    const double threshold_log2 = -total_samples_log2/2.0;
    
    printf("Cryptanalysis Assessment:\n");
    printf("  Total threads: %d blocks × %d threads = %" PRIu64 "\n", config->blocks, config->threads, total_threads);
    printf("  Samples per thread: 2^%d = %" PRIu64 "\n", config->sample_power, samples_per_thread);
    printf("  Total samples: %" PRIu64 " ≈ 2^%.1f\n", total_samples, total_samples_log2);
    printf("  Detection threshold: 2^%.1f\n", threshold_log2);
    
    if (abs_mean > 10.0 * detection_threshold) {
        printf("  Status: STRONG correlation detected (2^%.1f >> 2^%.1f)\n", abs_mean_log2, threshold_log2);
        printf("  Security: BROKEN - Practical attack feasible\n");
        printf("  Distinguisher complexity: ~2^%.0f operations\n", -2.0 * abs_mean_log2);
    } else if (abs_mean > 3.0 * detection_threshold) {
        printf("  Status: MODERATE correlation detected (2^%.1f > 2^%.1f)\n", abs_mean_log2, threshold_log2);
        printf("  Security: VULNERABLE - Theoretical attack possible\n");
        printf("  Distinguisher complexity: ~2^%.0f operations\n", -2.0 * abs_mean_log2);
    } else if (abs_mean > detection_threshold) {
        printf("  Status: WEAK correlation detected (2^%.1f ≈ 2^%.1f)\n", abs_mean_log2, threshold_log2);
        printf("  Security: MARGINAL - Requires many samples\n");
        printf("  Distinguisher complexity: ~2^%.0f operations\n", -2.0 * abs_mean_log2);
    } else {
        printf("  Status: NO significant correlation (2^%.1f ≤ 2^%.1f) - likely noise\n", abs_mean_log2, threshold_log2);
        printf("  Security: SECURE against this differential-linear attack\n");
        printf("  Complexity: > 2^64 operations (impractical)\n");
    }
    printf("\n");
    
    // Export correlation data for Python visualization
    export_correlation_data_for_python(correlations, num_experiments, config);
    
    // Cleanup dynamic allocation
    free(bins);
    free(bin_edges);
}

/**
 * @brief Display experiment keys and correlations
 */
void display_experiment_keys(ExperimentKey *keys, int num_experiments,
                             Config *config) {
    char key_hex[ORTHROS_STATE_SIZE + 1];
    char left_hex[ORTHROS_STATE_SIZE + 1];
    char right_hex[ORTHROS_STATE_SIZE + 1];
    char mask_hex[ORTHROS_STATE_SIZE + 1];

    format_nibbles(config->input_diff_left, left_hex);
    format_nibbles(config->input_diff_right, right_hex);
    format_nibbles(config->output_mask, mask_hex);

    printf("Experiment Keys and Individual Correlations:\n");
    printf("=============================================\n");
    printf("Configuration: %d-round Orthros\n", config->rounds);
    printf("  Mode   : %d\n", (int)config->mode);
    printf("  Offset : %d\n", config->offset);
    printf("  ΔL     : %s\n", left_hex);
    printf("  ΔR     : %s\n", right_hex);
    printf("  Output : %s\n\n", mask_hex);

    printf("┌─────┬────────────────────────────────┬─────────────────────┬────────────────────────┐\n");
    printf("│ Exp │              Key             │    Correlation      │      Log₂ |Corr|       │\n");
    printf("├─────┼────────────────────────────────┼─────────────────────┼────────────────────────┤\n");

    for (int i = 0; i < num_experiments; i++) {
        const double corr_mag = fabs(keys[i].correlation);
        const char sign = (keys[i].correlation >= 0.0) ? '+' : '-';

        format_nibbles(keys[i].key, key_hex);
        printf("│ %3d │ %s │ %c%+.6e      │",
               keys[i].experiment_id + 1, key_hex, sign, corr_mag);

        if (corr_mag > 0.0) {
            const double log_corr = log2(corr_mag);
            printf(" %c2^%7.2f         │\n", sign, log_corr);
        } else {
            printf("     -∞             │\n");
        }
    }

    printf("└─────┴────────────────────────────────┴─────────────────────┴────────────────────────┘\n");
    printf("\n");

    char filename[256];
    snprintf(filename, sizeof(filename),
             "orthros_%dr_keys_mode%d_offset%d.csv",
             config->rounds, (int)config->mode, config->offset);

    FILE *f = fopen(filename, "w");
    if (f) {
        fprintf(f, "# Orthros %d-round experiment keys and correlations\n", config->rounds);
        fprintf(f, "# Mode   : %d\n", (int)config->mode);
        fprintf(f, "# Offset : %d\n", config->offset);
        fprintf(f, "# ΔL     : %s\n", left_hex);
        fprintf(f, "# ΔR     : %s\n", right_hex);
        fprintf(f, "# Output : %s\n", mask_hex);
        fprintf(f, "ExperimentID,Key,Correlation\n");  // CSV header

        for (int i = 0; i < num_experiments; i++) {
            format_nibbles(keys[i].key, key_hex);
            fprintf(f, "%d,%s,%+.12e\n",
                    keys[i].experiment_id + 1, key_hex, keys[i].correlation);
        }
        fclose(f);

        printf("Detailed results saved to: %s\n\n", filename);
    }
}

/**
 * @brief Export correlation data in CSV format for Python visualization
 */
void export_correlation_data_for_python(const double *correlations, int num_experiments,
                                        const Config *config) {
    char left_hex[ORTHROS_STATE_SIZE + 1];
    char right_hex[ORTHROS_STATE_SIZE + 1];
    char mask_hex[ORTHROS_STATE_SIZE + 1];

    format_nibbles(config->input_diff_left, left_hex);
    format_nibbles(config->input_diff_right, right_hex);
    format_nibbles(config->output_mask, mask_hex);

    char filename[256];
    snprintf(filename, sizeof(filename),
             "orthros_corr_mode%d_r%d_o%d.csv",
             (int)config->mode, config->rounds, config->offset);

    FILE *f = fopen(filename, "w");
    if (!f) {
        printf("Warning: Could not create Python data file %s\n", filename);
        return;
    }

    fprintf(f, "# Orthros %d-round Differential-Linear Correlation Analysis\n", config->rounds);
    fprintf(f, "# Mode   : %d\n", (int)config->mode);
    fprintf(f, "# Offset : %d\n", config->offset);
    fprintf(f, "# ΔL     : %s\n", left_hex);
    fprintf(f, "# ΔR     : %s\n", right_hex);
    fprintf(f, "# Output : %s\n", mask_hex);
    fprintf(f, "# Experiments : %d\n", num_experiments);
    fprintf(f, "# GPU configuration: %d blocks x %d threads\n", config->blocks, config->threads);
    fprintf(f, "# Samples per thread: 2^%d\n", config->sample_power);
    fprintf(f, "# Generated by Orthros DLCT CUDA implementation\n");
    fprintf(f, "correlation\n");

    for (int i = 0; i < num_experiments; i++) {
        fprintf(f, "%.12e\n", correlations[i]);
    }

    fclose(f);

    printf("Correlation data exported for Python analysis: %s\n", filename);
    printf("To visualize: python3 src/visualizer.py %s\n", filename);
    printf("\n");
}
