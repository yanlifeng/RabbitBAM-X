#include "BamTools.h"
#include "libdeflate.h"
#include <htslib/hts_endian.h>
#include <cstring>
#include <climits>
#include <cassert>

// athread.h is available on Sunway platform
#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

// BlockProcessData structure for passing data to slave cores
// Simplified structure: only pointers, lengths, CRC, and result pair
struct BlockProcessData {
    const uint8_t *compressed_data;    // Aligned compressed data (comp->data + 18)
    size_t compressed_len;              // Compressed data length
    uint8_t *uncompressed_data;         // Aligned uncompressed data buffer
    size_t uncompressed_len;             // Uncompressed data buffer size (input) / actual size (output)
    uint32_t expected_crc;              // Pre-computed CRC
    std::pair<int, int> *result_pair;   // Result: (split_pos, bam_number)
};

// Helper function: le_to_u32 (slave core version with safe unaligned access)
// This version uses byte-by-byte reading to avoid alignment issues on Sunway slave cores
static inline uint32_t le_to_u32_slave(const uint8_t *buf) {
    // Read 4 bytes in little-endian order, safe for unaligned access
    return ((uint32_t) buf[0] |
            ((uint32_t) buf[1] << 8) |
            ((uint32_t) buf[2] << 16) |
            ((uint32_t) buf[3] << 24));
}

// Helper function: Rabbit_memcpy (local implementation for slave core)
inline void Rabbit_memcpy_slave(void *target, unsigned char *source, unsigned int length) {
    memcpy((uint8_t *) target, source, length);
}

// Helper function: bgzf_uncompress (simplified version for slave core)
// Note: libdeflate functions are available on slave cores
// Sunway slave cores require 64-byte alignment for memory access
// Data should already be aligned by host core before calling this function
int bgzf_uncompress_slave(uint8_t *dst, size_t *dlen, const uint8_t *src, size_t slen, uint32_t expected_crc) {
    // Check alignment for Sunway slave cores (64-byte alignment required)
    const uintptr_t ALIGN_MASK = 63;  // 64-byte alignment mask
    if (((uintptr_t)src & ALIGN_MASK) != 0 || ((uintptr_t)dst & ALIGN_MASK) != 0) {
        // Alignment check failed - this should not happen if host core prepared data correctly
        return -1;
    }
    
    // Use libdeflate functions available on slave cores (from libdeflate.h)
    struct libdeflate_decompressor *z = libdeflate_alloc_decompressor();
    if (!z) {
        return -1;
    }
    
    size_t actual_dlen = *dlen;
    int ret = libdeflate_deflate_decompress(z, src, slen, dst, actual_dlen, &actual_dlen);
    libdeflate_free_decompressor(z);
    
    if (ret != 0) {
        return -1;
    }
    
    *dlen = actual_dlen;
    
    uint32_t crc = libdeflate_crc32(0, (unsigned char *) dst, *dlen);
    if (crc != expected_crc) {
        return -2;
    }
    
    return 0;
}

// block_decode_func_slave implementation using simplified data structure
// This version takes aligned pointers directly, CRC is pre-computed
int block_decode_func_slave_simple(const uint8_t *compressed_data, size_t compressed_len,
                                    uint8_t *uncompressed_data, size_t *uncompressed_len,
                                    uint32_t expected_crc) {
    // Decompress using pre-aligned data and pre-computed CRC
    // printf("compressed data pointer: %p\n", compressed_data);
    // printf("compressed data length: %zu\n", compressed_len);
    // printf("Uncompressed data pointer: %p\n", uncompressed_data);
    // printf("Uncompressed data length: %zu\n", *uncompressed_len);
    int ret = bgzf_uncompress_slave(uncompressed_data, uncompressed_len,
                                     compressed_data, compressed_len, expected_crc);
    return ret;
}

// find_divide_pos_and_get_read_number implementation for slave core
std::pair<int, int> find_divide_pos_and_get_read_number_slave(bam_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ans = 0;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        Rabbit_memcpy_slave(&ret, block->data + divide_pos, 4);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {
                break;
            }
            Rabbit_memcpy_slave(x, block->data + divide_pos + 4, 32);
            int pos = (int32_t) x[1];
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;//block_len + c->l_extranul
            if (new_l_data > INT_MAX || l_qseq < 0 || l_qname < 1) {
                divide_pos += 4 + 32;
                continue;
            }
            if (((uint64_t) n_cigar << 2) + l_qname + l_extranul
                + (((uint64_t) l_qseq + 1) >> 1) + l_qseq > (uint64_t) new_l_data) {
                divide_pos += 4 + 32;
                continue;
            }
            if (divide_pos + 4 + 32 + l_qname > block->length) {
                break;
            }
            char fg_char;
            Rabbit_memcpy_slave(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
            if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {
                if (divide_pos + 4 + 32 + l_qname > block->length) {
                    break;
                }
                divide_pos += 4 + 32 + l_qname;
                continue;
            }
            
            if (divide_pos + 4 + ret > block->length) {
                break;
            }
            divide_pos += 4 + ret;
            ans++;
        } else {
            if (divide_pos + 4 > block->length) {
                break;
            }
            divide_pos += 4;
        }
    }
    return std::pair<int, int>(divide_pos, ans);
}

// Slave core function for parallel block processing
// Note: This function is compiled with -mslave flag for slave cores
extern "C" void slave_block_process(void *arg) {
    // Convert arg to BlockProcessData array pointer, then access element using _PEN
    BlockProcessData* array = (BlockProcessData *)arg;
    BlockProcessData* data = &array[_PEN];
    
    
    // Create temporary bam_block structure for find_divide_pos_and_get_read_number_slave
    // (it expects a bam_block*, but we only have data pointer and length)
    bam_block temp_block;
    temp_block.data = data->uncompressed_data;
    temp_block.length = BGZF_MAX_BLOCK_SIZE;
    temp_block.pos = 0;
    
    // Parallel decode using simplified interface
    int ret = block_decode_func_slave_simple(data->compressed_data, data->compressed_len,
                                              data->uncompressed_data, &data->uncompressed_len,
                                              data->expected_crc);
    if (ret) {
        assert(false);
        return;
    }
    // Update temp_block length after decompression
    temp_block.length = data->uncompressed_len;
    
    // Parallel pre-parse
    *(data->result_pair) = find_divide_pos_and_get_read_number_slave(&temp_block, 0);
}

// WriteCompressData structure for passing write block compression data to slave cores
struct WriteCompressData {
    uint8_t *uncompressed_data;      // Uncompressed data buffer (aligned)
    size_t uncompressed_len;          // Uncompressed data length
    uint8_t *compressed_data;         // Compressed data buffer (aligned)
    size_t *compressed_len;           // Compressed data length (output)
    int compress_level;                // Compression level
};

// Slave core function for parallel write block compression
extern "C" void slave_write_compress_process(void *arg) {
    WriteCompressData* array = (WriteCompressData *)arg;
    WriteCompressData* data = &array[_PEN];
    // Check alignment for Sunway slave cores (64-byte alignment required)
    const uintptr_t ALIGN_MASK = 63;
    if (((uintptr_t)data->uncompressed_data & ALIGN_MASK) != 0 || 
        ((uintptr_t)data->compressed_data & ALIGN_MASK) != 0) {
        // Alignment check failed
        *(data->compressed_len) = 0;
        return;
    }
    
    // Skip if input size is 0
    if (data->uncompressed_len == 0) {
        *(data->compressed_len) = 0;
        return;
    }
    
    // Use libdeflate compressor with gzip format
    struct libdeflate_compressor *compressor = libdeflate_alloc_compressor(data->compress_level);
    if (!compressor) {
        *(data->compressed_len) = 0;
        return;
    }
    
    // Get compressed size bound for gzip format
    size_t bound = libdeflate_gzip_compress_bound(compressor, data->uncompressed_len);
    
    // Perform gzip compression
    size_t out_size = libdeflate_gzip_compress(compressor, data->uncompressed_data, data->uncompressed_len,
                                               data->compressed_data, bound);
    
    libdeflate_free_compressor(compressor);
    
    // Set output size (0 if compression failed)
    *(data->compressed_len) = out_size;
}
