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
*/

#ifndef ORTHROS_COMMON_TYPES_H
#define ORTHROS_COMMON_TYPES_H

#include <stdint.h>
#include <stdbool.h>

#define ORTHROS_STATE_SIZE 32
#define ORTHROS_MAX_ROUNDS 12
#define ORTHROS_BATCH_SIZE 64

// Mode identifiers (matching CPU difflin.c constants)
typedef enum {
    ORTHROS_MODE_LEFT  = 0,
    ORTHROS_MODE_RIGHT = 1,
    ORTHROS_MODE_PRF   = 2
} OrthrosMode;

// Configuration structure for Orthros DLCT analysis
typedef struct {
    int blocks;
    int threads;
    int rounds;
    int offset;
    int experiments;
    int sample_power;
    OrthrosMode mode;
    bool verbose;
    bool self_test_only;
    bool generate_dlct;
    bool use_wildcard_mask;
    char wildcard_mask_pattern[ORTHROS_STATE_SIZE * 2 + 1];  // String pattern with 'x'
    unsigned char input_diff_left[ORTHROS_STATE_SIZE];
    unsigned char input_diff_right[ORTHROS_STATE_SIZE];
    unsigned char output_mask[ORTHROS_STATE_SIZE];
} Config;

// Structure for experiment key logging
typedef struct {
    int experiment_id;
    unsigned char key[ORTHROS_STATE_SIZE];
    double correlation;
} ExperimentKey;

#endif // ORTHROS_COMMON_TYPES_H