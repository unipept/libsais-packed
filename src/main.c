#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h> 

#include "bitpacking.h"
#include "libsais16x64.h"
#include "libsais32x64.h"
#include "libsais64.h"


void print_usage() {
    printf("Usage: ./build_ssa -s <sparseness> [-cdo] <input_file> <output_file>\n\n");
    printf("-s <sparseness>     : Defines the sparseness factor (an integer).\n");
    printf("-d                  : Flag to specify whether the input file is DNA data or protein data.\n");
    printf("-c                  : Flag to specify whether the output SA is compressed by bitpacking.\n");
    printf("-u                  : If enabled, the program will compute the SSA unoptimized, by computing the full SA and subsampling afterwards.\n");
    printf("<input_file>        : The path to the input file containing the DNA data.\n");
    printf("<output_file>       : The path where the output will be saved.\n");
}

uint8_t* read_text(char* input_fn, size_t* length) {

    // Open the input file for reading
    FILE *input_file = fopen(input_fn, "r");
    if (input_file == NULL) {
        perror("Failed to open file");
        print_usage();
        exit(1);
    }

    // Seek to the end to determine the length of the file
    fseek(input_file, 0, SEEK_END);
    *length = ftell(input_file);
    rewind(input_file);  // Return to the beginning of the file

    // Allocate memory for the text buffer
    uint8_t *text = (uint8_t *)malloc(*length + 1);  // +1 for the null-terminator
    if (text == NULL) {
        perror("Failed to allocate memory");
        print_usage();
        fclose(input_file);
        exit(1);
    }

    // Read the contents of the file into the buffer
    fread(text, 1, *length, input_file);
    text[*length] = '\0';  // Null-terminate the string
    fclose(input_file);

    return text;

}

int64_t* allocate_sa(size_t sa_length) {
    // Allocate memory for the suffix array (sa)
    int64_t *sa = (int64_t *)malloc(sa_length * sizeof(int64_t));
    if (sa == NULL) {
        perror("Failed to allocate memory for suffix array");
        print_usage();
        exit(1);
    }

    return sa;
}

int64_t* build_sa_optimized(uint8_t* text, size_t length, int64_t sparseness_factor, size_t sa_length, int dna) {
    int64_t* sa = allocate_sa(sa_length);

    uint8_t bits_per_char = dna > 0? BITS_PER_CHAR_DNA: BITS_PER_CHAR;

    int64_t required_bits = bits_per_char * sparseness_factor;
    
    if (sparseness_factor == 1) {

        libsais64(text, sa, sa_length, 0, NULL);
        free(text);
    
    } else if (required_bits <= 8) {
        
        uint8_t* packed_text = bitpack_text_8(text, length, sparseness_factor, sa_length, dna);
        free(text);
        libsais64(packed_text, sa, sa_length, 0, NULL);

    } else if (required_bits <= 16) {
        
        uint16_t* packed_text = bitpack_text_16(text, length, sparseness_factor, sa_length, dna);
        free(text);
        libsais16x64(packed_text, sa, sa_length, 0, NULL);

    } else if (required_bits <= 32) {
        
        uint32_t* packed_text = bitpack_text_32(text, length, sparseness_factor, sa_length, dna);
        free(text);
        libsais32x64(packed_text, sa, sa_length, 1 << required_bits, 0, NULL);

    } else {
        perror("Alphabet too big\n");
    }

    if (sparseness_factor > 1) {
        for (size_t i = 0; i < sa_length; i ++) {
            sa[i] *= sparseness_factor;
        }
    }

    return sa;
}

int64_t* build_sa(uint8_t* text, size_t length, int64_t sparseness_factor) {

    // Allocate memory for the suffix array (sa)
    int64_t *sa = (int64_t *)malloc(length * sizeof(int64_t));
    if (sa == NULL) {
        perror("Failed to allocate memory for suffix array");
        print_usage();
        free(text);
        exit(1);
    }

    libsais64(text, sa, length, 0, NULL);

    // Sample the suffix array
    if (sparseness_factor > 1) {
        size_t ssa_index = 0;
        for (size_t i = 0; i < length; i ++) {
            if (sa[i] % sparseness_factor == 0) {
                sa[ssa_index] = sa[i];
                ssa_index ++;
            }
        }
    }

    return sa;
}

void compress_sa(uint64_t* sa, size_t* sa_length, uint8_t bits_per_element) {

    int64_t element = 0;
    int8_t shift_element = 64 - bits_per_element;
    size_t compressed_i = 0;
    for (size_t i = 0; i < *sa_length; i ++) {
        if (shift_element < 0) { // new element does not fit in element
            element |= sa[i] >> (-1 * shift_element);
            sa[compressed_i] = element;
            compressed_i ++;
            element = 0;
            shift_element += 64;
        }
        element |= sa[i] << shift_element;
        shift_element -= bits_per_element;
    }

    sa[compressed_i] = element;

    *sa_length = compressed_i + 1;
}


uint64_t* decompress_sa(uint64_t* sa, size_t orig_sa_length, uint8_t bits_per_element) {
    uint64_t*  decompressed_sa = malloc(orig_sa_length * sizeof(int64_t));
    
    int8_t start_shift_element = 0;
    size_t compressed_i = 0;
    for (size_t i = 0; i < orig_sa_length; i ++) {
        decompressed_sa[i] |= (sa[compressed_i] << start_shift_element) >> (64 - bits_per_element);
        start_shift_element += bits_per_element;
        
        if (start_shift_element >= 64) {
            compressed_i ++;
            start_shift_element -= 64;
            if (start_shift_element > 0) { // new element does not fit in element
                decompressed_sa[i] |= sa[compressed_i] >> (64 - start_shift_element);
            }
        }
    }

    return decompressed_sa;
}

void write_sa(char* output_fn, uint8_t sparseness_factor, uint64_t* sa, size_t sa_length, int compressed) {
    // Open the output binary file for writing
    FILE *output_file = fopen(output_fn, "wb");
    if (output_file == NULL) {
        perror("Failed to open output file");
        print_usage();
        free(sa);
        exit(1);
    }

    // Write bits per element to the binary file
    uint8_t bits_per_element = 64;
    if (compressed > 0) {
        bits_per_element = (uint8_t) log2(sa_length * sparseness_factor) + 1;
    }
    fwrite(&bits_per_element, sizeof(uint8_t), 1, output_file);

    // Write sparseness factor to the binary file
    fwrite(&sparseness_factor, sizeof(uint8_t), 1, output_file);

    // Write length of sa to binary file
    uint64_t sa_length64 = (uint64_t) sa_length;
    fwrite(&sa_length64, sizeof(uint64_t), 1, output_file);

    // Write the suffix array to the binary file
    if (compressed > 0) {
        compress_sa(sa, &sa_length, bits_per_element);
    }
    fwrite(sa, sizeof(int64_t), sa_length, output_file);
    fclose(output_file);
}

void translate_L_to_I(uint8_t* text, size_t length) {
    for (size_t i = 0; i < length; i ++) {
        if (text[i] == 'L') {
            text[i] = 'I';
        }
    }
}

int main(int argc, char *argv[]) {
    printf("Command being executed: ");
    for (int i = 0; i < argc; i ++) {
        printf("%s ", argv[i]);
    }
    printf("\n");

    int opt;
    int compressed = 0, dna = 0, optimized = 1;
    char *sparseness = NULL;
    char *input_file = NULL;
    char *output_file = NULL;

    // Parse command-line options
    while ((opt = getopt(argc, argv, "s:cdu")) != -1) {
        switch (opt) {
            case 's': // Required argument
                sparseness = optarg;
                break;
            case 'c':
                compressed = 1;
                break;
            case 'd':
                dna = 1;
                break;
            case 'u':
                optimized = 0;
                break;
            default:
                print_usage();
                return EXIT_FAILURE;
        }
    }

    // Remaining arguments after options
    if (optind + 2 != argc) {
        fprintf(stderr, "Error: You must provide input and output files.\n");
        print_usage();
        return EXIT_FAILURE;
    }

    input_file = argv[optind];
    output_file = argv[optind + 1];

    // Validate required arguments
    if (!sparseness) {
        fprintf(stderr, "Error: Missing required option -s <sparseness>\n");
        print_usage();
        return EXIT_FAILURE;
    }

    int64_t sparseness_factor = atoi(sparseness);

    clock_t start_reading = clock();
    printf("Started reading input file from %s ...\n", input_file);
    size_t length;
    uint8_t* text = read_text(input_file, &length);
    if (dna <= 0) { // Translate L to I in protein alfabet
        translate_L_to_I(text, length);
    }
    printf("Done reading input file in %fs\n", ((double) clock() - start_reading) / CLOCKS_PER_SEC);

    clock_t start_sa = clock();
    printf("Started building SA...\n");
    size_t sparseness_factor_size = (size_t)sparseness_factor;
    size_t sa_length = (length + sparseness_factor_size - 1) / sparseness_factor_size;
    int64_t* sa;
    if (optimized > 0) {
        sa = build_sa_optimized(text, length, sparseness_factor, sa_length, dna);
    } else {
        sa = build_sa(text, length, sparseness_factor);
        free(text);
    }
    printf("Done building SA in %fs\n", ((double) clock() - start_sa) / CLOCKS_PER_SEC);

    clock_t start_writing = clock();
    printf("Started writing results...\n");
    write_sa(output_file, (uint8_t) sparseness_factor_size, (uint64_t*) sa, sa_length, compressed);
    printf("Done writing results to %s in %fs\n\n", output_file, ((double) clock() - start_writing) / CLOCKS_PER_SEC);

    return 0;
}

