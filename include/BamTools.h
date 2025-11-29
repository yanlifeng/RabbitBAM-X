#ifndef BAMSTATUS_BAMTOOLS_H
#define BAMSTATUS_BAMTOOLS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

#include <htslib/khash.h>
#include "header.h"
#include <stdint.h>
#include <limits.h>
#include "config.h"

// Define SIZE_MAX if not available (needed for some Sunway toolchains)
#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif
#include <cstdlib>
#include <mutex>
#include <libdeflate.h>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <vector>
#ifdef PLATFORM_X86
#include <omp.h>
#endif
#include <cmath>

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8
//#define BGZF_MAX_BLOCK_SIZE 0x10000
#define BGZF_MAX_BLOCK_COMPLETE_SIZE 0x40000
//#define BGZF_MAX_BLOCK_COMPLETE_SIZE 0x10000
#define THREAD_NUM_P 6

//#define use_parallel_read
//#define use_parallel_write

// Debug flags - uncomment to enable debug output
//#define DEBUG  // Enable all debug output

#include <sys/time.h>


inline double GetTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double) tv.tv_sec + (double) tv.tv_usec / 1000000;
}

inline uint64_t GetCycle() {
    unsigned int lo, hi;
    asm volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t) hi << 32) | lo;
}

// Aligned memory allocation helper function
// Allocates size bytes aligned to alignment (must be power of 2)
// Returns aligned pointer, stores original pointer before aligned address for deallocation
inline unsigned char* aligned_alloc_custom(size_t alignment, size_t size) {
    // Allocate extra space for alignment and storing original pointer
    size_t total_size = size + alignment + sizeof(void*);
    unsigned char* raw_ptr = new unsigned char[total_size];
    if (!raw_ptr) {
        return nullptr;
    }
    
    // Calculate aligned address
    uintptr_t raw_addr = (uintptr_t)raw_ptr;
    uintptr_t aligned_addr = (raw_addr + sizeof(void*) + alignment - 1) & ~(uintptr_t)(alignment - 1);
    unsigned char* aligned_ptr = (unsigned char*)aligned_addr;
    
    // Store original pointer before aligned address
    void** ptr_storage = (void**)(aligned_ptr - sizeof(void*));
    *ptr_storage = raw_ptr;
    
    return aligned_ptr;
}

// Free aligned memory allocated by aligned_alloc_custom
inline void aligned_free_custom(unsigned char* aligned_ptr) {
    if (aligned_ptr) {
        // Retrieve original pointer
        void** ptr_storage = (void**)(aligned_ptr - sizeof(void*));
        unsigned char* raw_ptr = (unsigned char*)(*ptr_storage);
        delete[] raw_ptr;
    }
}

typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t
)
static const uint8_t g_magic[19] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\0\0";

struct bam_block {
    unsigned int errcode;
    unsigned char *data;  // Changed from array to pointer for better alignment control
    unsigned int length;
    unsigned int pos;
    int64_t block_address;
    unsigned int split_pos;
    unsigned int bam_number;
};

struct bam_complete_block {
    unsigned int errcode;
    unsigned char *data = nullptr;
    unsigned int data_size;
    unsigned int length;
    unsigned int pos;
    int64_t block_address;
};

struct bam_write_block {
    int status = 0; // 0:uncompress 1:compress
    int block_num = -1; // -1 : not used correctly
    int block_offset;
    int block_length;
    uint8_t *uncompressed_data;
    uint8_t *compressed_data;
};

struct bgzf_cache_t {
    khash_t(cache) *h;
    khint_t last_pos;
};


typedef struct bam_block bam_block;
typedef struct bam_complete_block bam_complete_block;
typedef struct bam_write_block bam_write_block;


inline void packInt16(uint8_t *buffer, uint16_t value) {
    buffer[0] = value;
    buffer[1] = value >> 8;
}

inline void packInt32(uint8_t *buffer, uint32_t value) {
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

int sam_realloc_bam_data(bam1_t *b, size_t desired);

inline int realloc_bam_data(bam1_t *b, size_t desired);

inline int possibly_expand_bam_data(bam1_t *b, size_t bytes);

int bam_tag2cigar(bam1_t *b, int recal_bin,
                  int give_warning); // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen);

void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host);

inline int unpackInt16(const uint8_t *buffer);

int load_block_from_cache(BGZF *fp, int64_t block_address);

int check_header(const uint8_t *header);

const char *bgzf_zerr(int errnum, z_stream *zs);

int bgzf_uncompress(uint8_t *dst, size_t *dlen, const uint8_t *src, size_t slen, uint32_t expected_crc);

//int bgzf_uncompress(uint8_t *dst, unsigned int *dlen,const uint8_t *src, unsigned int slen,uint32_t expected_crc);
int fixup_missing_qname_nul(bam1_t *b);

int block_decode_func(struct bam_block *comp, struct bam_block *un_comp);

int read_block(BGZF *fp, struct bam_block *j);

int Rabbit_bgzf_read(struct bam_block *fq, void *data, unsigned int length);

void Rabbit_memcpy(void *target, unsigned char *source, unsigned int length);

int Rabbit_bgzf_read(struct bam_complete_block *fq, void *data, unsigned int length);

int read_bam(struct bam_block *fq, bam1_t *b, int is_be);

int read_bam(struct bam_complete_block *fq, bam1_t *b, int is_be);

int rabbit_bgzf_compress(void *_dst, size_t *dlen, const void *src, size_t slen, int level);

int rabbit_bgzf_gzip_compress(BGZF *fp, void *_dst, size_t *dlen, const void *src, size_t slen, int level);


int find_divide_pos(bam_block *block, int last_pos = 0);

int find_divide_pos(bam_complete_block *block, int last_pos = 0);


std::pair<int, int> find_divide_pos_and_get_read_number(bam_block *block, int last_pos = 0);

std::pair<int, int> find_divide_pos_and_get_read_number(bam_complete_block *block, int last_pos = 0);

int change_data_size(bam_complete_block *block);


#endif //BAMSTATUS_BAMTOOLS_H
