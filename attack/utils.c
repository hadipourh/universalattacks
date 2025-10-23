#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

// Cross-platform random support
#ifdef __linux__
    #include <sys/random.h>
    #ifndef GRND_NONBLOCK
    #define GRND_NONBLOCK 0x0001
    #endif
#endif

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

void init_prng(unsigned int offset) {
    unsigned int initial_seed = 0;
    
#ifdef __linux__
    // Linux: Use getrandom if available
    ssize_t temp = getrandom(&initial_seed, sizeof(initial_seed), GRND_NONBLOCK);
    if (temp == -1) {
        // Fallback to time-based seeding
        initial_seed = (unsigned int)time(NULL);
    }
#else
    // macOS and other systems: Use time-based seeding
    initial_seed = (unsigned int)time(NULL) + (unsigned int)getpid();
#endif
    
    initial_seed += offset;
    srand(initial_seed);
    // printf("[+] PRNG initialized to 0x%08X\n", initial_seed);
}

int isZero(unsigned char *diff, int size){
    for(int i=0;i<size;i++){
        if(diff[i]!=0){
            return 0;
        }
    }
    return 1;
}

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
