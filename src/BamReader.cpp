#include "BamReader.h"
#include <htslib/hts_endian.h>  // For le_to_u32

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#include <mutex>
// Forward declaration for slave function (C linkage)
extern "C" void slave_block_process(void *arg);
extern "C" void slave_write_compress_process(void *arg);
// Global mutex to coordinate reader and writer spawn
std::mutex g_athread_spawn_mutex;
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
#endif

void read_pack(BGZF *fp, BamRead *read) {
    bam_block *b;
    b = read->getEmpty();
    int count = 0;
    while (read_block(fp, b) == 0) {
        read->inputBlock(b);
        b = read->getEmpty();
    }
    read->ReadComplete();
}

#ifdef PLATFORM_SUNWAY
// Sunway version: parallel processing with athread
void compress_pack(BamRead *read, BamCompress *compress) {
    const int BATCH_SIZE = 64;
    pair<bam_block *, int> comp;
    bam_block *compressed_blocks[BATCH_SIZE];
    bam_block *uncompressed_blocks[BATCH_SIZE];
    int comp_values[BATCH_SIZE];
    std::pair<int, int> result_pairs[BATCH_SIZE];
    BlockProcessData process_data[BATCH_SIZE];
    int batch_count = 0;
    bool has_more_data = true;
    int batch_id = 0;
    
    // Aligned buffers for compressed data (64-byte aligned for Sunway)
    static uint8_t aligned_comp_buf[BATCH_SIZE][65536 + 64] __attribute__((aligned(64)));
    const uintptr_t ALIGN_MASK = 63;  // 64-byte alignment mask
    double tot_time = 0;
    double all_t_start = GetTime();
    
    while (has_more_data) {
        // Collect blocks until we have 64 or reach end
        batch_count = 0;
        while (batch_count < BATCH_SIZE) {
            comp = read->getReadBlock();
            if (comp.second < 0) {
                has_more_data = false;
                break;
            }
            compressed_blocks[batch_count] = comp.first;
            uncompressed_blocks[batch_count] = compress->getEmpty();
            comp_values[batch_count] = comp.second;
            
            // Pre-compute CRC on host core (before sending to slave)
            uint32_t crc = le_to_u32((uint8_t *) comp.first->data + comp.first->length - 8);
            
            // Get compressed data pointer (skip 18-byte header)
            const uint8_t *comp_data = comp.first->data + 18;
            size_t comp_len = comp.first->length - 18;
            
            // Ensure compressed data is 64-byte aligned for Sunway slave cores
            const uint8_t *aligned_comp_data = comp_data;
            bool need_align = ((uintptr_t)comp_data & ALIGN_MASK) != 0;
            if (need_align && comp_len <= sizeof(aligned_comp_buf[batch_count]) - 64) {
                memcpy(aligned_comp_buf[batch_count], comp_data, comp_len);
                aligned_comp_data = aligned_comp_buf[batch_count];
            }
            
            // Fill simplified BlockProcessData structure
            process_data[batch_count].compressed_data = aligned_comp_data;
            process_data[batch_count].compressed_len = comp_len;
            process_data[batch_count].uncompressed_data = uncompressed_blocks[batch_count]->data;
            process_data[batch_count].uncompressed_len = BGZF_MAX_BLOCK_SIZE;
            process_data[batch_count].expected_crc = crc;
            process_data[batch_count].result_pair = &result_pairs[batch_count];
            
            batch_count++;
        }
        
        if (batch_count == 0) {
            break;
        }
        
        double t_start = GetTime();
        
        // Spawn slave cores for parallel processing (with mutex to coordinate with writer)
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            __real_athread_spawn((void *)slave_block_process, process_data, 1);
        }
        
        // Wait for all slave cores to complete
        athread_join();
        
        double t_end = GetTime();
        tot_time += t_end - t_start;
        
        // Main core processes the results sequentially
        for (int i = 0; i < batch_count; i++) {
            read->backBlock(compressed_blocks[i]);
            
            // Update uncompressed block with results from slave cores
            uncompressed_blocks[i]->pos = 0;
            uncompressed_blocks[i]->length = process_data[i].uncompressed_len;
            uncompressed_blocks[i]->split_pos = result_pairs[i].first;
            uncompressed_blocks[i]->bam_number = result_pairs[i].second;
            
            compress->inputUnCompressData(uncompressed_blocks[i], comp_values[i]);
        }
        
        batch_id++;
    }
#ifdef DEBUG
    printf("DEBUG compress_pack total time: %lf\n", tot_time);
    printf("DEBUG compress_pack all time: %lf\n", GetTime() - all_t_start);
#endif
    
    compress->CompressThreadComplete();
}
#else
// x86 version: original sequential processing
void compress_pack(BamRead *read, BamCompress *compress) {
    pair < bam_block * , int > comp;
    bam_block *un_comp = compress->getEmpty();
    while (1) {
        comp = read->getReadBlock();
        if (comp.second < 0) {
            break;
        }
        block_decode_func(comp.first, un_comp);
        read->backBlock(comp.first);

        std::pair<int, int> tmp_pair = find_divide_pos_and_get_read_number(un_comp);
        un_comp->split_pos = tmp_pair.first, un_comp->bam_number = tmp_pair.second;
        compress->inputUnCompressData(un_comp, comp.second);
        un_comp = compress->getEmpty();
    }
    compress->CompressThreadComplete();
}
#endif

void assign_pack(BamCompress *compress, BamCompleteBlock *completeBlock) {
    bam_block *un_comp = nullptr;
    bam_complete_block *assign_block = completeBlock->getEmpty();
    int need_block_len = 0, input_length = 0;
    int last_use_block_length = 0;
    bool isclean = true;
    int ret = -1;
    static int block_count = 0;
    static int complete_block_count = 0;
    
#ifdef DEBUG
    printf("DEBUG assign_pack Starting thread\n");
#endif
    
    while (1) {
        if (isclean && un_comp != nullptr) {
            compress->backEmpty(un_comp);
        }
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
#ifdef DEBUG
            printf("DEBUG assign_pack No more un_comp blocks, exiting\n");
#endif
            break;
        }
        block_count++;
        ret = un_comp->split_pos;
        need_block_len = ret;
        if (assign_block->length + need_block_len > assign_block->data_size) {
            complete_block_count++;
#ifdef DEBUG
            printf("DEBUG assign_pack Block full, sending to completeBlock[%d] length=%u\n",
                   complete_block_count, assign_block->length);
#endif
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
#ifdef DEBUG
            if (assign_block == nullptr) {
                printf("DEBUG assign_pack WARNING: Got nullptr assign_block!\n");
            }
#endif
        }
        if (ret != un_comp->length) {

            memcpy(assign_block->data + assign_block->length, un_comp->data, ret * sizeof(char));
            assign_block->length += ret;
            complete_block_count++;
#ifdef DEBUG
            printf("DEBUG assign_pack Split block: copied %d bytes, sending to completeBlock[%d]\n",
                   ret, complete_block_count);
#endif
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();

            memcpy(assign_block->data + assign_block->length, un_comp->data + ret,
                   (un_comp->length - ret) * sizeof(char));
            assign_block->length += (un_comp->length - ret);
            compress->backEmpty(un_comp);
            un_comp = compress->getUnCompressData();
            if (un_comp == nullptr) {
#ifdef DEBUG
                printf("DEBUG assign_pack WARNING: Got nullptr un_comp after split!\n");
#endif
                break;
            }
            block_count++;
            memcpy(assign_block->data + assign_block->length, un_comp->data, un_comp->length * sizeof(char));
            assign_block->length += un_comp->length;


            int divide_pos = 0;
            int ret = 0;
            uint32_t x[8], new_l_data;
            while (divide_pos < assign_block->length) {
                Rabbit_memcpy(&ret, assign_block->data + divide_pos, 4);
                if (ret >= 32) {
                    if (divide_pos + 4 + 32 > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    Rabbit_memcpy(x, assign_block->data + divide_pos + 4, 32);
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
                    while (divide_pos + 4 + 32 + l_qname > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    char fg_char;
                    Rabbit_memcpy(&fg_char, assign_block->data + divide_pos + 4 + 32 + l_qname - 1, 1);
                    if (fg_char != '\0') {
                    }
                    if (fg_char != '\0' && l_extranul <= 0 && new_l_data > INT_MAX - 4) {

                        while (divide_pos + 4 + 32 + l_qname > assign_block->length) {
                            compress->backEmpty(un_comp);
                            un_comp = compress->getUnCompressData();
                            if (assign_block->length + un_comp->length > assign_block->data_size) {
                                change_data_size(assign_block);
                            }
                            memcpy(assign_block->data + assign_block->length, un_comp->data,
                                   un_comp->length * sizeof(char));
                            assign_block->length += un_comp->length;
                        }
                        divide_pos += 4 + 32 + l_qname;
                        continue;
                    }

                    while (divide_pos + 4 + ret > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos += 4 + ret;
                } else {
                    if (divide_pos + 4 > assign_block->length) {
                        compress->backEmpty(un_comp);
                        un_comp = compress->getUnCompressData();
                        if (assign_block->length + un_comp->length > assign_block->data_size) {
                            change_data_size(assign_block);
                        }
                        memcpy(assign_block->data + assign_block->length, un_comp->data,
                               un_comp->length * sizeof(char));
                        assign_block->length += un_comp->length;
                    }
                    divide_pos += 4;
                }

            }

            complete_block_count++;
#ifdef DEBUG
            printf("DEBUG assign_pack Finished parsing, sending assign_block[%d] length=%u\n",
                   complete_block_count, assign_block->length);
#endif
            completeBlock->inputCompleteBlock(assign_block);
            assign_block = completeBlock->getEmpty();
        } else {
            if (ret != un_comp->length) {
#ifdef DEBUG
                printf("DEBUG assign_pack ERROR: ret(%d) != un_comp->length(%u), breaking\n", ret, un_comp->length);
#endif
                break;
            }
            memcpy(assign_block->data + assign_block->length, un_comp->data, ret * sizeof(char));
            assign_block->length += ret;
            last_use_block_length = 0;
            isclean = true;
        }
    }
    if (assign_block != nullptr && assign_block->length != 0) {
        complete_block_count++;
#ifdef DEBUG
        printf("DEBUG assign_pack Final assign_block[%d] length=%u, sending to completeBlock\n",
               complete_block_count, assign_block->length);
#endif
        completeBlock->inputCompleteBlock(assign_block);
    }
#ifdef DEBUG
    printf("DEBUG assign_pack Finished: processed %d blocks, created %d complete blocks\n",
           block_count, complete_block_count);
#endif
    completeBlock->is_over();

}

void compress_test_pack(BamCompress *compress) {
    bam_block *un_comp = nullptr;
    while (1) {
        un_comp = compress->getUnCompressData();
        if (un_comp == nullptr) {
            break;
        }
        compress->backEmpty(un_comp);

    }
}


BamReader::BamReader(std::string file_name, int n_thread, bool is_tgs) {


    if ((sin = sam_open(file_name.c_str(), "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open file: %s\n", file_name.c_str());
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }


    read = new BamRead(1024);
    this->n_thread = n_thread;
    compress = new BamCompress(1024, n_thread);
    completeBlock = new BamCompleteBlock(256);


    read_thread = new thread(&read_pack, sin->fp.bgzf, read);
    compress_thread = new thread *[n_thread];

    for (int i = 0; i < n_thread; i++) {
        compress_thread[i] = new thread(&compress_pack, read, compress);
    }

    thread *assign_thread = new thread(&assign_pack, compress, completeBlock);

    if(is_tgs) {
        un_comp = completeBlock->getCompleteBlock();
#ifdef DEBUG
        printf("DEBUG BamReader get un_comp %p\n", (void*)un_comp);
#endif
    }
//
//#ifdef use_parallel_read
//#else
//    un_comp = completeBlock->getCompleteBlock();
//#endif


}


BamReader::BamReader(std::string file_name, int read_block, int compress_block, int compress_complete_block,
                     int n_thread, bool is_tgs) {

    if ((sin = sam_open(file_name.c_str(), "r")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open file: %s\n", file_name.c_str());
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
    }

    read = new BamRead(read_block);
    this->n_thread = n_thread;
    compress = new BamCompress(compress_block, n_thread);
    completeBlock = new BamCompleteBlock(compress_block);

    read_thread = new thread(&read_pack, sin->fp.bgzf, read);
    read_thread->join();
#ifdef DEBUG
    printf("DEBUG BamReader read_thread joined\n");
#endif
    compress_thread = new thread *[n_thread];

    for (int i = 0; i < n_thread; i++) {
        compress_thread[i] = new thread(&compress_pack, read, compress);
        compress_thread[i]->join();
#ifdef DEBUG
        printf("DEBUG BamReader compress_thread[%d] joined\n", i);
#endif
    }

    thread *assign_thread = new thread(&assign_pack, compress, completeBlock);

    if(is_tgs) {
        un_comp = completeBlock->getCompleteBlock();
    }

//#ifdef use_parallel_read
//#else
//    un_comp = completeBlock->getCompleteBlock();
//#endif
}


sam_hdr_t *BamReader::getHeader() {
    return hdr;
}

std::vector<bam1_t *> BamReader::getBam1_t_parallel(std::vector<bam1_t *> b_vec[THREAD_NUM_P]) {
#define THREAD_NUM_PR 4

    std::vector < bam1_t * > res_vec;

    bam_complete_block *blocks[THREAD_NUM_PR];

    int pre_vec_pos[THREAD_NUM_PR] = {0};
    int bam_num[THREAD_NUM_PR];

    for (int k = 0; k < 16; k++) {
        for (int i = 0; i < THREAD_NUM_PR; i++) {
            blocks[i] = completeBlock->getCompleteBlock();
        }


#ifdef PLATFORM_X86
#pragma omp parallel for num_threads(THREAD_NUM_PR) schedule(static)
#endif
        for (int i = 0; i < THREAD_NUM_PR; i++) {
            bam_num[i] = 0;
            if (blocks[i] == nullptr) continue;
            int now_num = pre_vec_pos[i];
            while (read_bam(blocks[i], b_vec[i][now_num], 0) >= 0) {
                now_num++;
                if (now_num >= b_vec[i].size()) {
                    fprintf(stderr, "ERROR: b_vec size not big enough: %d >= %zu\n", now_num, b_vec[i].size());
                    exit(1);
                }
            }
            bam_num[i] = now_num - pre_vec_pos[i];
        }

        for (int i = 0; i < THREAD_NUM_PR; i++) {
            if (blocks[i] != nullptr) completeBlock->backEmpty(blocks[i]);
        }

        for (int i = 0; i < THREAD_NUM_PR; i++) {
            res_vec.insert(res_vec.end(), b_vec[i].begin() + pre_vec_pos[i],
                           b_vec[i].begin() + pre_vec_pos[i] + bam_num[i]);
            pre_vec_pos[i] += bam_num[i];
        }
    }

    return res_vec;
}


bool BamReader::getBam1_t(bam1_t *b) {
    int ret;
    while (un_comp != nullptr) {
        if ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            return true;
        } else {
            completeBlock->backEmpty(un_comp);
            un_comp = completeBlock->getCompleteBlock();
        }
    }
    return false;
}


bam_complete_block *BamReader::getBamCompleteClock() {
    return completeBlock->getCompleteBlock();
}

void BamReader::backBamCompleteBlock(bam_complete_block *un_comp) {
    completeBlock->backEmpty(un_comp);
}
