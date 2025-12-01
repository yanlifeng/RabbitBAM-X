#include "BamTools.h"
#include "libdeflate.h"
#include <htslib/hts_endian.h>
#include <cstring>
#include <climits>
#include <cassert>

#ifdef PLATFORM_SUNWAY
#include <slave.h>
#endif

#define use_dynamic_task

#define cpe_num_slave 64
#define global_pen (_PEN)

#ifdef use_dynamic_task
//__uncached __cross long work_counter;
__uncached long work_counter;

int acquire_task(int block_num) {
    int cur_id;
    asm volatile("faal %0, 0(%1)\n\t"
                 : "=r"(cur_id)
                 : "r"(&work_counter)
                 : "memory");
    return cur_id;
}
#endif

struct BlockProcessData {
    const uint8_t *compressed_data;
    size_t compressed_len;
    uint8_t *uncompressed_data;
    size_t uncompressed_len;
    uint32_t expected_crc;
    std::pair<int, int> *result_pair;
};

// Wrapper structure to pass batch count along with data array
struct BlockProcessDataWrapper {
    BlockProcessData *data_array;
    int batch_count;
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

extern "C" void slave_read_process(void *arg) {
    BlockProcessDataWrapper* wrapper = (BlockProcessDataWrapper *)arg;
    BlockProcessData* array = wrapper->data_array;
    int batch_count = wrapper->batch_count;
    
#ifdef use_dynamic_task
    if(global_pen == 0) work_counter = 0;
    athread_ssync_array();
    
    int block_num = batch_count;
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        BlockProcessData* data = &array[i];
        
        bam_block temp_block;
        temp_block.data = data->uncompressed_data;
        temp_block.length = BGZF_MAX_BLOCK_SIZE;
        temp_block.pos = 0;
        
        int ret = slave_block_decompress(data->compressed_data, data->compressed_len,
                                          data->uncompressed_data, &data->uncompressed_len,
                                          data->expected_crc);
        if (ret) {
            assert(false);
            continue;
        }
        
        temp_block.length = data->uncompressed_len;
        *(data->result_pair) = slave_find_divide_pos(&temp_block, 0);
    }
#else
    // Static allocation: each slave core processes 16 tasks
    int tasks_per_core = batch_count / cpe_num_slave;
    int start_task = global_pen * tasks_per_core;
    int end_task = (global_pen == cpe_num_slave - 1) ? batch_count : (global_pen + 1) * tasks_per_core;
    
    for(int i = start_task; i < end_task; i++) {
        BlockProcessData* data = &array[i];
        
        bam_block temp_block;
        temp_block.data = data->uncompressed_data;
        temp_block.length = BGZF_MAX_BLOCK_SIZE;
        temp_block.pos = 0;
        
        int ret = slave_block_decompress(data->compressed_data, data->compressed_len,
                                          data->uncompressed_data, &data->uncompressed_len,
                                          data->expected_crc);
        if (ret) {
            assert(false);
            continue;
        }
        
        temp_block.length = data->uncompressed_len;
        *(data->result_pair) = slave_find_divide_pos(&temp_block, 0);
    }
#endif
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

// Wrapper structure to pass batch count along with data array
struct WriteCompressDataWrapper {
    WriteCompressData *data_array;
    int batch_count;
};

extern "C" void slave_write_process(void *arg) {
    WriteCompressDataWrapper* wrapper = (WriteCompressDataWrapper *)arg;
    WriteCompressData* array = wrapper->data_array;
    int batch_count = wrapper->batch_count;
    
#ifdef use_dynamic_task
    if(global_pen == 0) work_counter = 0;
    athread_ssync_array();
    
    int block_num = batch_count;
    for(int i = acquire_task(block_num); i < block_num; i = acquire_task(block_num)) {
        WriteCompressData* data = &array[i];
        
        slave_block_compress(data->uncompressed_data, data->uncompressed_len,
                             data->deflate_buffer, data->deflate_len,
                             data->compress_level);
    }
#else
    // Static allocation: each slave core processes 16 tasks
    int tasks_per_core = batch_count / cpe_num_slave;
    int start_task = global_pen * tasks_per_core;
    int end_task = (global_pen == cpe_num_slave - 1) ? batch_count : (global_pen + 1) * tasks_per_core;
    
    for(int i = start_task; i < end_task; i++) {
        WriteCompressData* data = &array[i];
        
        slave_block_compress(data->uncompressed_data, data->uncompressed_len,
                             data->deflate_buffer, data->deflate_len,
                             data->compress_level);
    }
#endif
}
