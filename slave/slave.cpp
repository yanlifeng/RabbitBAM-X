#include "BamTools.h"
#include "libdeflate.h"
#include <htslib/hts_endian.h>
#include <cstring>
#include <climits>
#include <cassert>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

struct BlockProcessData {
    const uint8_t *compressed_data;
    size_t compressed_len;
    uint8_t *uncompressed_data;
    size_t uncompressed_len;
    uint32_t expected_crc;
    std::pair<int, int> *result_pair;
};

// BGZF decompression for slave core
int slave_bgzf_decompress(uint8_t *dst, size_t *dlen, const uint8_t *src, size_t slen, uint32_t expected_crc) {
    const uintptr_t ALIGN_MASK = 63;
    if (((uintptr_t)src & ALIGN_MASK) != 0 || ((uintptr_t)dst & ALIGN_MASK) != 0) {
        return -1;
    }
    
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

int slave_block_decompress(const uint8_t *compressed_data, size_t compressed_len,
                           uint8_t *uncompressed_data, size_t *uncompressed_len,
                           uint32_t expected_crc) {
    return slave_bgzf_decompress(uncompressed_data, uncompressed_len,
                                  compressed_data, compressed_len, expected_crc);
}

std::pair<int, int> slave_find_divide_pos(bam_block *block, int last_pos) {
    int divide_pos = last_pos;
    int ans = 0;
    int ret = 0;
    uint32_t x[8], new_l_data;
    while (divide_pos < block->length) {
        memcpy(&ret, block->data + divide_pos, 4);
        if (ret >= 32) {
            if (divide_pos + 4 + 32 > block->length) {
                break;
            }
            memcpy(x, block->data + divide_pos + 4, 32);
            int l_qname = x[2] & 0xff;
            int l_extranul = (l_qname % 4 != 0) ? (4 - l_qname % 4) : 0;
            int n_cigar = x[3] & 0xffff;
            int l_qseq = x[4];
            new_l_data = ret - 32 + l_extranul;
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
            memcpy(&fg_char, block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
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

extern "C" void read_process(void *arg) {
    BlockProcessData* array = (BlockProcessData *)arg;
    BlockProcessData* data = &array[_PEN];
    
    bam_block temp_block;
    temp_block.data = data->uncompressed_data;
    temp_block.length = BGZF_MAX_BLOCK_SIZE;
    temp_block.pos = 0;
    
    int ret = slave_block_decompress(data->compressed_data, data->compressed_len,
                                      data->uncompressed_data, &data->uncompressed_len,
                                      data->expected_crc);
    if (ret) {
        assert(false);
        return;
    }
    
    temp_block.length = data->uncompressed_len;
    *(data->result_pair) = slave_find_divide_pos(&temp_block, 0);
}

int slave_bgzf_compress(uint8_t *uncompressed_data, size_t uncompressed_len,
                        uint8_t *deflate_buffer, size_t *deflate_len,
                        int compress_level) {
    const uintptr_t ALIGN_MASK = 63;
    if (((uintptr_t)uncompressed_data & ALIGN_MASK) != 0 || 
        ((uintptr_t)deflate_buffer & ALIGN_MASK) != 0) {
        *deflate_len = 0;
        return -1;
    }
    
    if (uncompressed_len == 0) {
        *deflate_len = 0;
        return 0;
    }
    
    struct libdeflate_compressor *compressor = libdeflate_alloc_compressor(compress_level);
    if (!compressor) {
        *deflate_len = 0;
        return -1;
    }
    
    // Compress to DEFLATE format (raw DEFLATE, not gzip)
    // Same as main core: compress to buffer, host will add BGZF header/footer
    size_t deflate_bound = libdeflate_deflate_compress_bound(compressor, uncompressed_len);
    size_t deflate_out_size = libdeflate_deflate_compress(compressor, uncompressed_data, uncompressed_len,
                                                          deflate_buffer, deflate_bound);
    
    libdeflate_free_compressor(compressor);
    
    if (deflate_out_size == 0) {
        *deflate_len = 0;
        return -1;
    }
    
    *deflate_len = deflate_out_size;
    return 0;
}

int slave_block_compress(uint8_t *uncompressed_data, size_t uncompressed_len,
                         uint8_t *deflate_buffer, size_t *deflate_len,
                         int compress_level) {
    return slave_bgzf_compress(uncompressed_data, uncompressed_len,
                                deflate_buffer, deflate_len, compress_level);
}

struct WriteCompressData {
    uint8_t *uncompressed_data;
    size_t uncompressed_len;
    uint8_t *deflate_buffer;
    size_t *deflate_len;
    int compress_level;
};

extern "C" void slave_write_process(void *arg) {
    WriteCompressData* array = (WriteCompressData *)arg;
    WriteCompressData* data = &array[_PEN];
    
    slave_block_compress(data->uncompressed_data, data->uncompressed_len,
                         data->deflate_buffer, data->deflate_len,
                         data->compress_level);
}
