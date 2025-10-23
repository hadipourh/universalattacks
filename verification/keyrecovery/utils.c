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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __APPLE__
#include <unistd.h>
#include <fcntl.h>
#else
#include <sys/random.h>
#endif
#include <regex.h>
#include<ctype.h>


void printState(const unsigned char S[32], const char* str) {
    if (str != NULL) {
        printf("%s> ", str); // Print the string
    }
    for (size_t i = 0; i < 32; ++i) {
        printf("%X ", S[i] & 0xF); // Print each nibble in hex
    }
    // printf("\n");
}

int hammingWeight(int* arr, size_t size) {
    int weight = 0;
    for (size_t i = 0; i < size; ++i) {
        weight += arr[i];
    }
    return weight;
}

int compareArrays(const unsigned char arr1[], const unsigned char arr2[], size_t n) {
    for (size_t i = 0; i < n; ++i) {
        if (arr1[i] != arr2[i]) {
            // Arrays are different
            printf("Different\n");
            return 0;
        }
    }
    // Arrays are the same
    printf("Equivalent\n");
    return 1;
}

void setNibbles(unsigned char S[32], int pos[][8], unsigned char nib[][2], size_t m) {
    for (size_t i = 0; i < m; ++i) {
        unsigned char nib0 = nib[i][0];
        unsigned char nib1 = nib[i][1];

        // printf("nib   %u %u\n", nib0,nib1);

        // Clear the bits at the specified positions
        for (size_t j = 0; j < 4; ++j) {
            S[pos[i][j] / 4] &= ~(1 << (3 - pos[i][j] % 4)); // Clear the bit at the correct position
        }

        for (size_t j = 4; j < 8; ++j) {
            S[pos[i][j] / 4] &= ~(1 << (3 - pos[i][j] % 4)); // Clear the bit at the correct position
        }

        // Set the bits for nib0
        for (size_t j = 0; j < 4; ++j) {
            S[pos[i][j] / 4] |= ((nib0 >> (3 - j)) & 1) << (3 - pos[i][j] % 4); // Set the bit at the correct position
        }

        // Set the bits for nib1
        for (size_t j = 4; j < 8; ++j) {
            S[pos[i][j] / 4] |= ((nib1 >> (3 - (j - 4))) & 1) << (3 - pos[i][j] % 4); // Set the bit at the correct position
        }
    }
    printState(S,"\n S ");
}

int isKeyContradict(int keyBitPos[][8], unsigned char nib[][2], size_t numActNib){
    // for(int i=0;i<numActNib;i++){
    //         printf("in nib> %d %d\n", nib[i][0], nib[i][1]);
    // }
    int tempKeyContradict=0;
    for(int compNib=0;compNib<numActNib;compNib++){
        for(int compIndex=0;compIndex<8;compIndex++){
            for(int targNib=0;targNib<numActNib;targNib++){
                for(int targIndex=0;targIndex<8;targIndex++){
                    if(keyBitPos[compNib][compIndex]==keyBitPos[targNib][targIndex]){
                        // printf("\n(%d,%d) (%d,%d)--- %d-%d>> %d - %d",compNib,compIndex,targNib,targIndex,keyBitPos[compNib][compIndex],keyBitPos[targNib][targIndex],(nib[compNib][compIndex/4]>>(3-compIndex%4)) & 0x01, (nib[targNib][targIndex/4]>>(3-targIndex%4)) & 0x01);                        
                        if(((nib[compNib][compIndex/4]>>(3-compIndex%4)) & 0x01)!= ((nib[targNib][targIndex/4]>>(3-targIndex%4)) & 0x01)){
                            tempKeyContradict=1;
                        }
                    }
                }
            }
        }
    }
    // printf("\n");
    return tempKeyContradict;
}


unsigned char stringToNibble(const char *str){
    int val = atoi(str); // Convert string to integer
    if (val < 0 || val > 15) {
        printf("Input must be between 0 and 15.\n");
        exit(1);
    }
    return (unsigned char)val; // Convert integer to unsigned char (nibble)
}

void trimString(char *str){
    // Trim leading spaces
    char *start = str;
    while (*start && isspace(*start)) {
        ++start;
    }

    // Trim trailing spaces
    char *end = str + strlen(str) - 1;
    while (end > start && isspace(*end)) {
        *end-- = '\0';
    }

    // Shift the trimmed string to the beginning of the original string
    if (start != str) {
        memmove(str, start, end - start + 2); // Add 2 to include the null terminator
    }
}

int hex_to_int(char c) {
    if (c >= '0' && c <= '9') {
        return c - '0';
    } else if (c >= 'a' && c <= 'f') {
        return c - 'a' + 10;
    } else if (c >= 'A' && c <= 'F') {
        return c - 'A' + 10;
    }
    return -1;  // Invalid hex character
}

int hex_string_to_bytes(char *hex_string, unsigned char *byte_array, size_t array_size) {
    size_t len = strlen(hex_string);
    if (len != array_size) {
        return -1;  // Invalid hex string length or array size
    }

    for (size_t i = 0; i < len; i ++) {
        int high = hex_to_int(hex_string[i]);
        if (high == -1) {
            return -1;  // Invalid hex character
        }
        byte_array[i] = (unsigned char)(high);
    }

    return 0;  // Success
}

//----------------------------------
// Dot product
//----------------------------------

int dot_prod(unsigned char *A, unsigned char *B, size_t array_size){
    int output=0;
    for (size_t i=0;i<array_size;i++){
        char C=(A[i] & B[i]);
        while (C>0){
            output ^= (C & 0x1);
            C >>= 1;
        }
    }
    return output;
}

//----------------------------------
// Init PRNG
//----------------------------------

unsigned int init_prng(unsigned int offset) {
    unsigned int initial_seed = 0;
    ssize_t temp;
    
#ifdef __APPLE__
    // macOS alternative using /dev/urandom
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd >= 0) {
        temp = read(fd, &initial_seed, sizeof(initial_seed));
        close(fd);
    } else {
        temp = -1;
    }
#else
    temp = getrandom(&initial_seed, sizeof(initial_seed), 0);
#endif
    
    //temp=5;
    if (temp == -1) perror("error!");
    initial_seed += offset;
    // srand(initial_seed);
    // srand(time(NULL));
	printf("[+] PRNG initialized to 0x%08X\n", initial_seed);
    return initial_seed;
}


int isZero(unsigned char *diff, int size){
    for(int i=0;i<size;i++){
        if(diff[i]!=0){
            return 0;
        }
    }
    return 1;
}

char generate_random_hexchar() {
    int random_num = rand() % 16; // Generate a random number between 0 and 15
    if (random_num < 10) {
        return '0' + random_num; // For numbers 0-9
    } else {
        return 'a' + (random_num - 10); // For letters a-f
    }
}


#define MAX_LINE_LENGTH 1024
#define INITIAL_BUFFER_SIZE 10

// Function to read a file line by line and return the lines
char **read_file_line_by_line(const char *filename, int *num_lines) {
    FILE *file;
    char line[MAX_LINE_LENGTH];
    char **lines = NULL;
    int buffer_size = INITIAL_BUFFER_SIZE;
    int count = 0;

    // Open the file for reading
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Allocate initial buffer for lines
    lines = (char **)malloc((size_t)buffer_size * sizeof(char *));
    if (lines == NULL) {
        perror("Error allocating memory");
        fclose(file);
        return NULL;
    }

    // Read the file line by line
    while (fgets(line, sizeof(line), file)) {
        // Remove the newline character
        line[strcspn(line, "\n")] = '\0';

        // Reallocate memory if necessary
        if (count >= buffer_size) {
            buffer_size *= 2;
            char **tmp = (char **)realloc(lines, (size_t)buffer_size * sizeof(char *));
            if (tmp == NULL) {
                perror("Error reallocating memory");
                fclose(file);
                for (int i = 0; i < count; ++i) {
                    free(lines[i]);
                }
                free(lines);
                return NULL;
            }
            lines = tmp;
        }

        // Allocate memory for the line and copy it
        lines[count] = (char *)malloc((strlen(line) + 1) * sizeof(char));
        if (lines[count] == NULL) {
            perror("Error allocating memory for line");
            fclose(file);
            for (int i = 0; i < count; ++i) {
                free(lines[i]);
            }
            free(lines);
            return NULL;
        }
        strcpy(lines[count], line);
        count++;
    }

    // Close the file
    fclose(file);

    // Set the number of lines read
    *num_lines = count;

    return lines;
}

#define MAX_MATCHES 10
#define MAX_SPLIT_PARTS 100

// Function to split a string based on a pattern
char **split_string(const char *string, const char *pattern, int *num_parts) {
    regex_t regex;
    regmatch_t matches[MAX_MATCHES];
    int ret, start, end;
    int buffer_size = MAX_SPLIT_PARTS;
    int count = 0;
    char **parts = (char **)malloc((size_t)buffer_size * sizeof(char *));
    const char *current = string;
    if (parts == NULL) {
        perror("Error allocating memory for split parts");
        return NULL;
    }

    // Compile the regular expression
    ret = regcomp(&regex, pattern, REG_EXTENDED);
    if (ret) {
        fprintf(stderr, "Could not compile regex\n");
        free(parts);
        return NULL;
    }

    // Execute the regular expression and split the string
    while ((ret = regexec(&regex, current, MAX_MATCHES, matches, 0)) == 0) {
        start = matches[0].rm_so;
        end = matches[0].rm_eo;

        // Allocate memory for the substring and copy it
        if (start > 0) {
            parts[count] = (char *)malloc((size_t)(start + 1) * sizeof(char));
            if (parts[count] == NULL) {
                perror("Error allocating memory for split part");
                for (int i = 0; i < count; ++i) {
                    free(parts[i]);
                }
                free(parts);
                regfree(&regex);
                return NULL;
            }
            strncpy(parts[count], current, (size_t)start);
            parts[count][start] = '\0';
            count++;

            if (count >= buffer_size) {
                buffer_size *= 2;
                char **tmp = (char **)realloc(parts, (size_t)buffer_size * sizeof(char *));
                if (tmp == NULL) {
                    perror("Error reallocating split parts");
                    for (int i = 0; i < count; ++i) {
                        free(parts[i]);
                    }
                    free(parts);
                    regfree(&regex);
                    return NULL;
                }
                parts = tmp;
            }
        }

        current += end;
    }

    // Copy the remaining part of the string
    if (*current != '\0') {
        parts[count] = strdup(current);
        if (parts[count] == NULL) {
            perror("Error duplicating remaining string");
            for (int i = 0; i < count; ++i) {
                free(parts[i]);
            }
            free(parts);
            regfree(&regex);
            return NULL;
        }
        count++;
    }

    regfree(&regex);
    *num_parts = count;
    return parts;
}

// Function to retrieve substrings between two specific characters
char** retrieve_items_between_chars(const char *string, char start_char, char end_char, int *count) {
    const char *start = string;
    const char *end;
    char **substrings = NULL;
    int substr_count = 0;

    while ((start = strchr(start, start_char)) != NULL) {
        start++; // Move past start_char
        end = strchr(start, end_char);
        
        if (end != NULL) {
            // Calculate the length of the substring
            int len = end - start;
            
            // Allocate memory for the new substring and copy it
            char *substring = (char *)malloc((size_t)len + 1);
            if (substring == NULL) {
                perror("Failed to allocate substring");
                for (int i = 0; i < substr_count; ++i) {
                    free(substrings[i]);
                }
                free(substrings);
                return NULL;
            }
            strncpy(substring, start, (size_t)len);
            substring[len] = '\0';
            
            // Reallocate memory for the array of substrings
            char **tmp = (char **)realloc(substrings, (size_t)(substr_count + 1) * sizeof(char*));
            if (tmp == NULL) {
                perror("Failed to reallocate memory");
                free(substring);
                for (int i = 0; i < substr_count; ++i) {
                    free(substrings[i]);
                }
                free(substrings);
                return NULL;
            }
            substrings = tmp;
            substrings[substr_count] = substring;
            substr_count++;
            start = end + 1; // Move past end_char
        } else {
            break; // No more end_char found
        }
    }
    
    *count = substr_count;
    return substrings;
}

int containsSubstring(const char *str, const char *sub){
    int len1 = strlen(str);
    int len2 = strlen(sub);
    
    for (int i = 0; i <= len1 - len2; i++)
    {
        int j;
        for (j = 0; j < len2; j++)
        {
            if (str[i + j] != sub[j])
                break;
        }
        if (j == len2)
            return 1; // Substring found
    }
    return 0; // Substring not found
}


// Function to write a variable string to a file using a file pointer
int write_string_to_file(FILE *file, const char *string_to_write) {
    if (file == NULL) {
        perror("Error: File pointer is NULL");
        return 1;
    }

    if (fprintf(file, "%s\n", string_to_write) < 0) {
        perror("Error writing to file");
        return 1;
    }

    return 0; // Success
}

// int main() {
//     const char *string = "one,two,three,four";
//     const char *pattern = ",";
//     int num_parts, i;

//     char **parts = split_string(string, pattern, &num_parts);
//     if (parts == NULL) {
//         printf("Failed to split the string.\n");
//         return 1;
//     }

//     for (i = 0; i < num_parts; i++) {
//         printf("%s\n", parts[i]);
//         free(parts[i]);
//     }

//     free(parts);
//     return 0;
// }

// int main() {
//     int num_lines, i, num_parts, count;
//     char **lines = read_file_line_by_line("input.txt", &num_lines);
//     const char *pattern = "(ad|da| )";
//     char start_char = 's';
//     char end_char = 'd';


//     if (lines == NULL) {
//         printf("Failed to read lines from file\n");
//         return 1;
//     }

//     // Print the lines
//     for (i = 0; i < num_lines; i++) {
//         char **parts = split_string(lines[i], pattern, &num_parts);
//         printf("'%s'\n", lines[i]);

//         char **substrings = retrieve_items_between_chars(lines[i], start_char, end_char, &count);
        
//         if (substrings != NULL) {
//             for (i = 0; i < count; i++) {
//                 printf("Substring: %s\n", substrings[i]);
//                 free(substrings[i]); // Free each individual substring
//             }
//             free(substrings); // Free the array of substring pointers
//         }

//         for (int j = 0; j < num_parts; j++) {
//             printf("%s\n", parts[j]);
//             free(parts[j]);
//         }
//         free(lines[i]); // Free each line
//     }

//     free(lines); // Free the array of pointers

//     return 0;
// }

// int main() {
//     const char *hex_string = "a947436710924ccd47f2d571deea8f05";
//     unsigned char byte_array[32];

//     if (hex_string_to_bytes(hex_string, byte_array, sizeof(byte_array)) == 0) {
//         // Print the result
//         for (size_t i = 0; i < sizeof(byte_array); i++) {
//             printf("%x", byte_array[i]);
//         }
//         printf("\n");
//     } else {
//         printf("Error: Invalid hex string or array size.\n");
//     }

//     return 0;
// }
