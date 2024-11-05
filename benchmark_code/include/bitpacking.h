

#ifndef BITPACKING_H
#define BITPACKING_H

uint8_t* bitpack_text_8(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t *packed_len, int dna);

uint16_t* bitpack_text_16(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t *packed_len, int dna);

uint32_t* bitpack_text_32(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t *packed_len, int dna);

int64_t* bitpack_text_64(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t *packed_len, int dna);

#endif
