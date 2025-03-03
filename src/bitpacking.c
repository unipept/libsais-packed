#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "bitpacking.h"

uint8_t* build_char_to_rank(const uint8_t* text, size_t text_len, uint8_t* alphabet_size) {

    uint8_t occuring[256] = {0};
    for (size_t i = 0; i < text_len; i++) {
        if (! occuring[text[i]]) {
            occuring[text[i]] = 1;
        }
    }

    uint8_t* char_to_rank = calloc(256, 1);
    uint8_t chars = 0;
    for (int i = 0; i < 256; i++) {
        if (occuring[i]) {
            char_to_rank[i] = chars;
            chars ++;
        }
    }

    *alphabet_size = chars;

    return char_to_rank;

}

// Function to get the rank of a character in a proteomic alphabet
uint8_t get_rank_aa(uint8_t c) {
    switch (c) {
        case '$': return 0;
        case '-': return 1;
        default: return 2 + (c - 'A');
    }
}

// Function to get the rank of a character in a nucleotide alphabet
uint8_t get_rank_dna(uint8_t c) {
    switch(c) {
        case '$': return 0;
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }
    printf("other char %c \n", c);
    return 0;
}


// Function to bit-pack the text based on the sparseness factor
uint8_t* bitpack_text_8(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char) {
    size_t sparseness_factor_size = (size_t)sparseness_factor;
        
    uint8_t *text_packed = (uint8_t *)calloc(packed_len, sizeof(uint8_t));
    if (text_len == 0 || text_packed == NULL) {
        return text_packed;
    }

    for (size_t i = 0; i < (packed_len - 1); i++) {
        size_t ti = i * sparseness_factor_size;
        uint8_t element = 0;
        for (size_t j = 0; j < sparseness_factor_size; j++) {
            uint8_t c = text[ti + j];
            uint8_t rank_c = char_to_rank[c];
            element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    uint8_t last_element = 0;
    size_t last_el_start = sparseness_factor_size * (packed_len - 1);
    for (size_t i = 0; i < (text_len - 1) % sparseness_factor_size + 1; i++) {
        uint8_t c = text[last_el_start + i];
        uint8_t rank_c = char_to_rank[c];
        last_element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - i));
    }
    text_packed[packed_len - 1] = last_element;

    return text_packed;
}

// Function to bit-pack the text based on the sparseness factor
uint16_t* bitpack_text_16(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char) {
    size_t sparseness_factor_size = (size_t)sparseness_factor;
        
    uint16_t *text_packed = (uint16_t *)calloc(packed_len, sizeof(uint16_t));
    if (text_len == 0 || text_packed == NULL) {
        return text_packed;
    }

    for (size_t i = 0; i < (packed_len - 1); i++) {
        size_t ti = i * sparseness_factor_size;
        uint16_t element = 0;
        for (size_t j = 0; j < sparseness_factor_size; j++) {
            uint8_t c = text[ti + j];
            uint16_t rank_c = (uint16_t) char_to_rank[c];
            element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    uint16_t last_element = 0;
    size_t last_el_start = sparseness_factor_size * (packed_len - 1);
    for (size_t i = 0; i < (text_len - 1) % sparseness_factor_size + 1; i++) {
        uint8_t c = text[last_el_start + i];
        uint16_t rank_c = (uint16_t) char_to_rank[c];
        last_element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - i));
    }
    text_packed[packed_len - 1] = last_element;

    return text_packed;
}


// Function to bit-pack the text based on the sparseness factor
uint32_t* bitpack_text_32(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char) {
    size_t sparseness_factor_size = (size_t)sparseness_factor;
    
    uint32_t *text_packed = (uint32_t *)calloc(packed_len, sizeof(uint32_t));
    if (text_len == 0 || text_packed == NULL) {
        return text_packed;
    }

    for (size_t i = 0; i < (packed_len - 1); i++) {
        size_t ti = i * sparseness_factor_size;
        uint32_t element = 0;
        for (size_t j = 0; j < sparseness_factor_size; j++) {
            uint8_t c = text[ti + j];
            uint32_t rank_c = (uint32_t) char_to_rank[c];
            element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    uint32_t last_element = 0;
    size_t last_el_start = sparseness_factor_size * (packed_len - 1);
    for (size_t i = 0; i < (text_len - 1) % sparseness_factor_size + 1; i++) {
        uint8_t c = text[last_el_start + i];
        uint32_t rank_c = (uint32_t) char_to_rank[c];
        last_element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - i));
    }
    text_packed[packed_len - 1] = last_element;

    return text_packed;
}

// Function to bit-pack the text based on the sparseness factor
int64_t* bitpack_text_64(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char) {
    size_t sparseness_factor_size = (size_t)sparseness_factor;
    
    int64_t *text_packed = (int64_t *)calloc(packed_len, sizeof(int64_t));
    if (text_len == 0 || text_packed == NULL) {
        return text_packed;
    }

    for (size_t i = 0; i < (packed_len - 1); i++) {
        size_t ti = i * sparseness_factor_size;
        int64_t element = 0;
        for (size_t j = 0; j < sparseness_factor_size; j++) {
            uint8_t c = text[ti + j];
            uint64_t rank_c = (uint64_t) char_to_rank[c];
            element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - j));
        }
        text_packed[i] = element;
    }

    // Handle the last element
    int64_t last_element = 0;
    size_t last_el_start = sparseness_factor_size * (packed_len - 1);
    for (size_t i = 0; i < (text_len - 1) % sparseness_factor_size + 1; i++) {
        uint8_t c = text[last_el_start + i];
        uint64_t rank_c = (uint64_t) char_to_rank[c];
        last_element |= rank_c << (bits_per_char * (sparseness_factor_size - 1 - i));
    }
    text_packed[packed_len - 1] = last_element;

    return text_packed;
}