

#ifndef BITPACKING_H
#define BITPACKING_H

uint8_t* build_char_to_rank(const uint8_t* text, size_t text_len, uint8_t* alphabet_size);

uint8_t* bitpack_text_8(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char);

uint16_t* bitpack_text_16(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char);

uint32_t* bitpack_text_32(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char);

int64_t* bitpack_text_64(const uint8_t *text, size_t text_len, uint8_t sparseness_factor, size_t packed_len, uint8_t* char_to_rank, uint8_t bits_per_char);

#endif
