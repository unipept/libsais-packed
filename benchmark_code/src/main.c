#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "bitpacking.h"
#include "libsais16x64.h"
#include "libsais32x64.h"
#include "libsais64.h"


void print_help() {
    printf("Usage: ./libsais64 <input_file> <output_file>\n");
}

int main(int argc, char *argv[]) {

    if (argc != 3) {
        print_help();
        return 1;
    }

    // Open the input file for reading
    FILE *input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        perror("Failed to open file");
        print_help();
        return 1;
    }

    // Seek to the end to determine the length of the file
    fseek(input_file, 0, SEEK_END);
    size_t length = ftell(input_file);
    rewind(input_file);  // Return to the beginning of the file

    // Allocate memory for the text buffer
    uint8_t *text = (uint8_t *)malloc(length + 1);  // +1 for the null-terminator
    if (text == NULL) {
        perror("Failed to allocate memory");
        print_help();
        fclose(input_file);
        return 1;
    }

    // Read the contents of the file into the buffer
    fread(text, 1, length, input_file);
    text[length] = '\0';  // Null-terminate the string
    fclose(input_file);

    size_t packed_len = 0;
    int64_t sparseness_factor = 3;
    uint32_t* packed_text = bitpack_text_32(text, length, sparseness_factor, &packed_len, 0);

    // Allocate memory for the suffix array (sa)
    int64_t *sa = (int64_t *)malloc(packed_len * sizeof(int64_t));
    if (sa == NULL) {
        perror("Failed to allocate memory for suffix array");
        print_help();
        free(text);
        return 1;
    }

    libsais32x64(packed_text, sa, packed_len, 1 << 15, 0, NULL);

    // Open the output binary file for writing
    FILE *output_file = fopen(argv[2], "wb");
    if (output_file == NULL) {
        perror("Failed to open output file");
        print_help();
        free(sa);
        free(text);
        return 1;
    }

    // Write the suffix array to the binary file
    fwrite(sa, sizeof(uint64_t), packed_len, output_file);

    // Clean up
    fclose(output_file);
    free(sa);
    free(text);
    return 0;
}

