#include "BamWriter.h"

#ifdef PLATFORM_SUNWAY
#include <athread.h>
#include <mutex>
extern "C" void slave_write_process(void *arg);
// External mutex from BamReader.cpp
extern std::mutex g_athread_spawn_mutex;
// WriteCompressData structure (must match slave/slave.cpp)
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
#endif

int rabbit_write_deflate_block(BGZF *fp, bam_write_block *write_block) {
    size_t comp_size = BGZF_MAX_BLOCK_SIZE;
    int ret;
    if (!fp->is_gzip)
        ret = rabbit_bgzf_compress(write_block->compressed_data, &comp_size, write_block->uncompressed_data,
                                   write_block->block_offset, fp->compress_level);
    else
        ret = rabbit_bgzf_gzip_compress(fp, write_block->compressed_data, &comp_size, write_block->uncompressed_data,
                                        write_block->block_offset, fp->compress_level);

    if (ret != 0) {
        hts_log_debug("Compression error %d", ret);
        fp->errcode |= BGZF_ERR_ZLIB;
        return -1;
    }
    return comp_size;
}

int rabbit_bgzf_flush(BGZF *fp, bam_write_block *write_block) {
    while (write_block->block_offset > 0) {
        int block_length;
#ifdef DEBUG
        printf("DEBUG rabbit_bgzf_flush Block offset: %d\n", write_block->block_offset);
#endif
        block_length = rabbit_write_deflate_block(fp, write_block);
        if (block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, write_block->compressed_data, block_length) != block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }


        write_block->block_offset = 0;
        fp->block_address += block_length;
    }
    write_block->block_offset = 0;
    return 0;
}

int rabbit_bgzf_mul_flush(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block) {

    bam_write_compress->inputUnCompressData(write_block);
    write_block = bam_write_compress->getEmpty();
    return 0;
}

int rabbit_bgzf_write(BGZF *fp, bam_write_block *&write_block, const void *data, size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            if (rabbit_bgzf_flush(fp, write_block) != 0) return -1;
        }
    }
    return length - remaining;
}

int rabbit_bgzf_mul_write_fast(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block,
                               const void *data, size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        //if (write_block->block_offset == BGZF_BLOCK_SIZE) {
        //    if (rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block) != 0) return -1;
        //}
    }
    return length - remaining;
}


int
rabbit_bgzf_mul_write(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, const void *data,
                      size_t length) {
    const uint8_t *input = (const uint8_t *) data;
    ssize_t remaining = length;
    while (remaining > 0) {
        uint8_t *buffer = (uint8_t *) write_block->uncompressed_data;
        int copy_length = BGZF_BLOCK_SIZE - write_block->block_offset;
        if (copy_length > remaining) copy_length = remaining;
        memcpy(buffer + write_block->block_offset, input, copy_length);
        write_block->block_offset += copy_length;
        input += copy_length;
        remaining -= copy_length;
        if (write_block->block_offset == BGZF_BLOCK_SIZE) {
            if (rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block) != 0) return -1;
        }
    }
    return length - remaining;
}

int rabbit_bgzf_flush_try(BGZF *fp, bam_write_block *write_block, ssize_t size) {
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_flush(fp, write_block);
    }
    return 0;
}

int
rabbit_bgzf_mul_flush_try(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, ssize_t size) {
    if (write_block->block_offset + size > BGZF_BLOCK_SIZE) {
        return rabbit_bgzf_mul_flush(fp, bam_write_compress, write_block);
    }
    return 0;
}

int bam_write_pack(BGZF *fp, BamWriteCompress *bam_write_compress) {
    bam_write_block *block;
    int write_count = 0;
#ifdef DEBUG
    printf("DEBUG bam_write_pack Output thread started\n");
#endif
    while (1) {
        block = bam_write_compress->getCompressData();
        if (block == nullptr) {
#ifdef DEBUG
            printf("DEBUG bam_write_pack No more blocks, exiting. Total written: %d\n", write_count);
#endif
            return 0;
        }
        if (block->block_length < 0) {
            hts_log_debug("Deflate block operation failed: %s", bgzf_zerr(block->block_length, NULL));
            return -1;
        }
        if (hwrite(fp->fp, block->compressed_data, block->block_length) != block->block_length) {
            hts_log_error("File write failed (wrong size)");
            fp->errcode |= BGZF_ERR_IO; // possibly truncated file
            return -1;
        }
        write_count++;
#ifdef DEBUG
        printf("DEBUG bam_write_pack Wrote block %d length=%d total=%d\n",
               block->block_num, block->block_length, write_count);
#endif
        block->block_offset = 0;
        fp->block_address += block->block_length;
        bam_write_compress->backEmpty(block);
    }

}

#ifdef PLATFORM_SUNWAY
// Sunway version: batch compression with slave cores
void bam_write_compress_pack_sunway(BGZF *fp, BamWriteCompress *bam_write_compress) {
    const int BATCH_SIZE = 64 * 16;
    bam_write_block *blocks[BATCH_SIZE];
    WriteCompressData process_data[BATCH_SIZE];
    size_t compressed_lens[BATCH_SIZE];  // Temporary array for size_t values
    int batch_count = 0;
    bool has_more_data = true;
    int batch_id = 0;
    
#ifdef DEBUG
    printf("DEBUG bam_write_compress_pack_sunway Compression thread started\n");
#endif
    while (has_more_data) {
        // Collect blocks until we have 64*16 (1024) or reach end
        batch_count = 0;
        while (batch_count < BATCH_SIZE) {
            blocks[batch_count] = bam_write_compress->getUnCompressData();
            if (blocks[batch_count] == nullptr) {
                has_more_data = false;
                break;
            }
            
            // Fill WriteCompressData structure
            // From slave core: compress DEFLATE data to aligned buffer (compressed_data is 64-byte aligned)
            // Host core will copy DEFLATE data to compressed_data + BLOCK_HEADER_LENGTH and add header/footer
            process_data[batch_count].uncompressed_data = blocks[batch_count]->uncompressed_data;
            process_data[batch_count].uncompressed_len = blocks[batch_count]->block_offset;
            process_data[batch_count].deflate_buffer = blocks[batch_count]->compressed_data;  // Use aligned start of buffer
            process_data[batch_count].deflate_len = &compressed_lens[batch_count];
            process_data[batch_count].compress_level = fp->compress_level;
            
            batch_count++;
        }
        
        if (batch_count == 0) {
#ifdef DEBUG
            printf("DEBUG bam_write_compress_pack_sunway No more blocks, exiting\n");
#endif
            break;
        }
        
        // Spawn slave cores for parallel compression (with mutex to coordinate with reader)
        {
            std::lock_guard<std::mutex> lock(g_athread_spawn_mutex);
            WriteCompressDataWrapper wrapper;
            wrapper.data_array = process_data;
            wrapper.batch_count = batch_count;
            __real_athread_spawn((void *)slave_write_process, &wrapper, 1);
        }
        
        // Wait for all slave cores to complete
        athread_join();
        
        // Process results - assemble complete BGZF block from slave core DEFLATE data
        for (int i = 0; i < batch_count; i++) {
            size_t deflate_len = compressed_lens[i];
            if (deflate_len == 0) {
                blocks[i]->block_length = -1;
                continue;
            }
            
            uint8_t *dst = blocks[i]->compressed_data;
            size_t uncompressed_len = blocks[i]->block_offset;
            
            // Slave core compressed DEFLATE data to the start of compressed_data (64-byte aligned)
            // Now we need to: 1) Copy DEFLATE data to correct position, 2) Add header, 3) Add footer
            // Copy DEFLATE data first (from offset 0 to offset 18) - use memmove for overlapping regions
            uint8_t *deflate_src = dst;  // Slave core compressed here
            uint8_t *deflate_dst = dst + BLOCK_HEADER_LENGTH;  // Final position
            // memmove handles overlapping correctly (copying from back to front when dst > src)
            if (deflate_dst != deflate_src) {
                memmove(deflate_dst, deflate_src, deflate_len);
            }
            
            // Add BGZF header (this won't overwrite DEFLATE data since header is only 18 bytes)
            memcpy(dst, g_magic, BLOCK_HEADER_LENGTH);
            
            // Calculate total block length (header + deflate data + footer)
            size_t total_len = deflate_len + BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
            packInt16(&dst[16], total_len - 1);
            
            // Add BGZF footer
            uint32_t crc = libdeflate_crc32(0, blocks[i]->uncompressed_data, uncompressed_len);
            packInt32(&dst[total_len - 8], crc);
            packInt32(&dst[total_len - 4], uncompressed_len);
            
            blocks[i]->block_length = (int)total_len;
            bam_write_compress->inputCompressData(blocks[i]);
        }
        batch_id++;
    }
#ifdef DEBUG
    printf("DEBUG bam_write_compress_pack_sunway Compression thread completed, total batches=%d\n", batch_id);
#endif
    bam_write_compress->CompressThreadComplete();
}
#endif

void bam_write_compress_pack(BGZF *fp, BamWriteCompress *bam_write_compress) {
#ifdef PLATFORM_SUNWAY
    // Use Sunway version with slave cores
    bam_write_compress_pack_sunway(fp, bam_write_compress);
#else
    // Original x86 version
    bam_write_block *block;
    while (1) {
        block = bam_write_compress->getUnCompressData();
        if (block == nullptr) {
            break;
        }
        block->block_length = rabbit_write_deflate_block(fp, block);
        bam_write_compress->inputCompressData(block);
    }
    bam_write_compress->CompressThreadComplete();
#endif
}

int rabbit_bam_write_test(BGZF *fp, bam_write_block *write_block, bam1_t *b) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR
    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t) c->bin << 16 | c->qual << 8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t) c->flag << 16 | 2;
    else x[3] = (uint32_t) c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (rabbit_bgzf_flush_try(fp, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_write(fp, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_write(fp, write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok) ok = (rabbit_bgzf_write(fp, write_block, b->data + c->l_qname, b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"
            PRIhts_pos
            " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                    bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t *) bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t) c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t) cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, buf, 8) >= 0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_en], b->l_data - cigar_en) >=
                  0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, buf, 4) >= 0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_write(fp, write_block, &b->data[cigar_st], c->n_cigar * 4) >= 0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}

int
rabbit_bam_write_mul_test(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, bam1_t *b) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR


    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t) c->bin << 16 | c->qual << 8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t) c->flag << 16 | 2;
    else x[3] = (uint32_t) c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (rabbit_bgzf_mul_flush_try(fp, bam_write_compress, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, x, 32) >= 0);
    }
    if (ok) ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, b->data, c->l_qname - c->l_extranul) >= 0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok)
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, b->data + c->l_qname,
                                        b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"
            PRIhts_pos
            " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                    bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t *) bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t) c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t) cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, buf, 8) >=
                  0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &b->data[cigar_en],
                                        b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, buf, 4) >=
                  0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write(fp, bam_write_compress, write_block, &b->data[cigar_st], c->n_cigar * 4) >=
                  0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}


void benchmark_write_pack(BamCompleteBlock *completeBlock, samFile *output, sam_hdr_t *hdr, int level) {

    uint8_t *compress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    uint8_t *uncompress_block_test = new uint8_t[BGZF_BLOCK_SIZE];
    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->uncompressed_block = uncompress_block_test;
    output->fp.bgzf->compressed_block = compress_block_test;
    output->fp.bgzf->compress_level = level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0) {
        fprintf(stderr, "ERROR: Header write failed\n");
        return;
    }
    bam_complete_block *un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block = new bam_write_block();
    write_block->block_offset = 0;
    write_block->uncompressed_data = new uint8_t[BGZF_BLOCK_SIZE];
    write_block->compressed_data = new uint8_t[BGZF_BLOCK_SIZE];
    write_block->status = 0;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            rabbit_bam_write_test(output->fp.bgzf, write_block, b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    rabbit_bgzf_flush(output->fp.bgzf, write_block);
#ifdef DEBUG
    printf("DEBUG BamWriter Bam1_t Number: %lld\n", ans);
    printf("DEBUG BamWriter Block Number: %lld\n", res);
#endif
}

void benchmark_write_mul_pack(BamCompleteBlock *completeBlock, BamWriteCompress *bam_write_compress, samFile *output,
                              sam_hdr_t *hdr, int level) {

    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;
    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }
    if (sam_hdr_write(output, hdr) != 0) {
        fprintf(stderr, "ERROR: Header write failed\n");
        return;
    }
    bam_complete_block *un_comp;
    long long ans = 0;
    long long res = 0;
    bam_write_block *write_block = bam_write_compress->getEmpty();
    int bam_num = 1;
    while (1) {
        un_comp = completeBlock->getCompleteBlock();
        if (un_comp == nullptr) {
            break;
        }
        int ret;
        while ((ret = (read_bam(un_comp, b, 0))) >= 0) {
            rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);
            ans++;
        }
        res++;
        completeBlock->backEmpty(un_comp);
    }
    if (write_block->block_offset > 0) {
        bam_write_compress->inputUnCompressData(write_block);
    }
    bam_write_compress->WriteComplete();
}

int
rabbit_bam_write_mul_parallel(BGZF *fp, BamWriteCompress *bam_write_compress, bam_write_block *&write_block, bam1_t *b,
                              std::vector<bam_write_block *> &block_vec) {
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok = 1;
    if (c->l_qname - c->l_extranul > 255) {
        hts_log_error("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        errno = EOVERFLOW;
        return -1;
    }
    if (c->n_cigar > 0xffff) block_len += 16; // "16" for "CGBI", 4-byte tag length and 8-byte fake CIGAR


    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        hts_log_error("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t) c->bin << 16 | c->qual << 8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t) c->flag << 16 | 2;
    else x[3] = (uint32_t) c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    if (write_block->block_offset + 4 + block_len >= BGZF_BLOCK_SIZE) {
        block_vec.push_back(write_block);
        write_block = bam_write_compress->getEmpty();
    }
    //ok = (rabbit_bgzf_mul_flush_try(fp, bam_write_compress, write_block, 4 + block_len) >= 0);
    if (fp->is_be) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, ed_swap_4p(&y), 4) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &block_len, 4) >= 0);
        }
    }
    if (ok) {
        ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, x, 32) >= 0);
    }
    if (ok)
        ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, b->data, c->l_qname - c->l_extranul) >=
              0);
    if (c->n_cigar <= 0xffff) { // no long CIGAR; write normally
        if (ok)
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, b->data + c->l_qname,
                                             b->l_data - c->l_qname) >= 0);
    } else { // with long CIGAR, insert a fake CIGAR record and move the real CIGAR to the CG:B,I tag
        uint8_t buf[8];
        uint32_t cigar_st, cigar_en, cigar[2];
        hts_pos_t cigreflen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(b));
        if (cigreflen >= (1 << 28)) {
            // Length of reference covered is greater than the biggest
            // CIGAR operation currently allowed.
            hts_log_error("Record %s with %d CIGAR ops and ref length %"
            PRIhts_pos
            " cannot be written in BAM.  Try writing SAM or CRAM instead.\n",
                    bam_get_qname(b), c->n_cigar, cigreflen);
            return -1;
        }
        cigar_st = (uint8_t *) bam_get_cigar(b) - b->data;
        cigar_en = cigar_st + c->n_cigar * 4;
        cigar[0] = (uint32_t) c->l_qseq << 4 | BAM_CSOFT_CLIP;
        cigar[1] = (uint32_t) cigreflen << 4 | BAM_CREF_SKIP;
        u32_to_le(cigar[0], buf);
        u32_to_le(cigar[1], buf + 4);
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, buf, 8) >=
                  0); // write cigar: <read_length>S<ref_length>N
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &b->data[cigar_en],
                                             b->l_data - cigar_en) >= 0); // write data after CIGAR
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, "CGBI", 4) >= 0); // write CG:B,I
        }
        u32_to_le(c->n_cigar, buf);
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, buf, 4) >=
                  0); // write the true CIGAR length
        }
        if (ok) {
            ok = (rabbit_bgzf_mul_write_fast(fp, bam_write_compress, write_block, &b->data[cigar_st], c->n_cigar * 4) >=
                  0); // write the real CIGAR
        }
    }
    if (fp->is_be) swap_data(c, b->l_data, b->data, 0);
    return ok ? 4 + block_len : -1;
}

void BamWriter::write_parallel(std::vector<bam1_t *> b_vec) {

#define THREAD_NUM_PW 4
    int vec_size = b_vec.size();
    std::vector < bam_write_block * > todo_push_blocks[THREAD_NUM_PW];
    for (int i = 0; i < THREAD_NUM_PW; i++) {
        blocks[i] = bam_write_compress->getEmpty();
    }

#ifdef PLATFORM_X86
#pragma omp parallel for num_threads(THREAD_NUM_PW) schedule(static)
    for (int i = 0; i < vec_size; i++) {
        int tid = omp_get_thread_num();
        rabbit_bam_write_mul_parallel(output->fp.bgzf, bam_write_compress, blocks[tid], b_vec[i],
                                      todo_push_blocks[tid]);
    }
#else
    // Sunway version: sequential processing
    for (int i = 0; i < vec_size; i++) {
        int tid = 0;  // Single thread for Sunway
        rabbit_bam_write_mul_parallel(output->fp.bgzf, bam_write_compress, blocks[tid], b_vec[i],
                                      todo_push_blocks[tid]);
    }
#endif

    for (int i = 0; i < THREAD_NUM_PW; i++) {
        for (auto item: todo_push_blocks[i]) {
            bam_write_compress->inputUnCompressData(item);
        }
        bam_write_compress->inputUnCompressData(blocks[i]);
    }


}


void BamWriter::write(bam1_t *b) {
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);

}

void BamWriter::bam_write(bam1_t *b) {
    rabbit_bam_write_mul_test(output->fp.bgzf, bam_write_compress, write_block, b);
}


BamWriter::BamWriter(int threadNumber, int level, int BufferSize) {

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


}


BamWriter::BamWriter(std::string filename, int threadNumber, int level, int BufferSize, bool is_tgs) {

    if ((output = sam_open(filename.c_str(), "wb")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open file\n");
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++)
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);


    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;

    if(is_tgs) {
        write_block=bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif
    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}

}

BamWriter::BamWriter(std::string filename, sam_hdr_t *hdr, int threadNumber, int level, int BufferSize, bool is_tgs) {

    if ((output = sam_open(filename.c_str(), "wb")) == NULL) {
        fprintf(stderr, "ERROR: Cannot open file\n");
    }
    if (sam_hdr_write(output, hdr) != 0) {
        fprintf(stderr, "ERROR: Header write failed\n");
    }

    n_thread_write = threadNumber;
    bam_write_compress = new BamWriteCompress(BufferSize, n_thread_write);


    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++) {
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);
    }

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);

    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = level;


    if(is_tgs) {
        write_block = bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif
    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}


}

void BamWriter::hdr_write(sam_hdr_t *hdr) {
    if (sam_hdr_write(output, hdr) != 0) {
        fprintf(stderr, "ERROR: Header write failed\n");
    }
}


void BamWriter::set_output(samFile *output, bool is_tgs) {
    this->output = output;
    write_compress_thread = new std::thread *[n_thread_write];
    for (int i = 0; i < n_thread_write; i++)
        write_compress_thread[i] = new std::thread(&bam_write_compress_pack, output->fp.bgzf, bam_write_compress);

    write_output_thread = new std::thread(&bam_write_pack, output->fp.bgzf, bam_write_compress);


    output->fp.bgzf->block_offset = 0;
    output->fp.bgzf->compress_level = 6;

    if(is_tgs) {
        write_block=bam_write_compress->getEmpty();
    }

//#ifdef use_parallel_write
//#else
//    write_block=bam_write_compress->getEmpty();
//#endif

    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    blocks[i] = bam_write_compress->getEmpty();
    //}


}

void BamWriter::over_parallel() {

    //for(int i = 0; i < THREAD_NUM_PW; i++) {
    //    if (blocks[i]->block_offset > 0) {
    //        bam_write_compress->inputUnCompressData(blocks[i]);
    //    }
    //}
    bam_write_compress->WriteComplete();
    for (int i = 0; i < n_thread_write; i++) write_compress_thread[i]->join();
    write_output_thread->join();
    int ret = sam_close(output);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        //exit_code = EXIT_FAILURE;
    }

}


void BamWriter::over() {

    if (write_block->block_offset > 0) {
        bam_write_compress->inputUnCompressData(write_block);
        write_block = bam_write_compress->getEmpty();
    }
    bam_write_compress->WriteComplete();
    for (int i = 0; i < n_thread_write; i++) write_compress_thread[i]->join();
    write_output_thread->join();
    int ret = sam_close(output);
    if (ret < 0) {
        fprintf(stderr, "Error closing output.\n");
        //exit_code = EXIT_FAILURE;
    }
}

