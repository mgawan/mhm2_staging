const double Q2Perror[] = { 1.0,
                            0.7943,	0.6309,	0.5012,	0.3981,	0.3162,
                            0.2512,	0.1995,	0.1585,	0.1259,	0.1,
                            0.07943,	0.06310,	0.05012,	0.03981,	0.03162,
                            0.02512,	0.01995,	0.01585,	0.01259,	0.01,
                            0.007943,	0.006310,	0.005012,	0.003981,	0.003162,
                            0.002512,	0.001995,	0.001585,	0.001259,	0.001,
                            0.0007943,	0.0006310,	0.0005012,	0.0003981,	0.0003162,
                            0.0002512,	0.0001995,	0.0001585,	0.0001259,	0.0001,
                            7.943e-05,	6.310e-05,	5.012e-05,	3.981e-05,	3.162e-05,
                            2.512e-05,	1.995e-05,	1.585e-05,	1.259e-05,	1e-05,
                            7.943e-06,	6.310e-06,	5.012e-06,	3.981e-06,	3.162e-06,
                            2.512e-06,	1.995e-06,	1.585e-06,	1.259e-06,	1e-06,
                            7.943e-07,	6.310e-07,	5.012e-07,	3.981e-07,	3.1622e-07,
                            2.512e-07,	1.995e-07,	1.585e-07,	1.259e-07,	1e-07,
                            7.943e-08,	6.310e-08,	5.012e-08,	3.981e-08,	3.1622e-08,
                            2.512e-08,	1.995e-08,	1.585e-08,	1.259e-08,	1e-08};

// This is a very simple algorithm (eg compared to bbmerge). We just take the longest perfect overlap
// between the two reads. If there is no perfect overlap, then there is no merge.
// the reads_fname will be replaced with the new name of the checkpoint
// modified to allow some % errors in the overlap based on quality score differences
// and pick the base from the read with the higher quality score
static void merge_reads_file(fq_reader_t fqr, fq_reader_t fqr2, char *merged_reads_checkpoint, char * libName, int files_per_pair, int hexify_libi, int qual_offset) {
    if (files_per_pair == 0) DIE("Invalid program logic, can not merge reads that are not paired!\n");
    START_TIMER(T_MERGE_READS);
    Buffer id_buf = initBuffer(MAX_READ_NAME_LEN), seq_buf = initBuffer(1024), quals_buf = initBuffer(1024);
    Buffer prev_id_buf = initBuffer(MAX_READ_NAME_LEN), prev_seq_buf = initBuffer(1024), prev_quals_buf = initBuffer(1024);
    Buffer tmpseq_buf = initBuffer(1024), tmpquals_buf = initBuffer(1024);
    int64_t num_read_pairs = 0, num_reads = 0, num_merged = 0, num_ambiguous = 0;
    uint64_t read_id_1 = MYTHREAD, read_id_2 = MYTHREAD;
    int64_t bytes_written = 0;
    int64_t mergedLen = 0, overlapLen = 0, maxReadLen = 0;
    const int16_t MIN_OVERLAP = 12;
    const int16_t EXTRA_TEST_OVERLAP = 2;
    const int16_t MAX_MISMATCHES = 3; // allow up to 3 mismatches, with MAX_PERROR
    const int Q2PerrorSize = sizeof(Q2Perror) / sizeof(*Q2Perror);
    assert(qual_offset == 33 || qual_offset == 64);

    // illumina reads generally accumulate errors at the end, so allow more mismatches in the overlap as long as differential quality indicates a clear winner
    const double MAX_PERROR = 0.025; // allow up to 2.5% accumulated mismatch probablity of error within the overlap by differential quality score
    const int16_t EXTRA_MISMATCHES_PER_1000 = (int) 150; // allow additional mismatches per 1000 bases of overlap before aborting an overlap test
    const uint8_t MAX_MATCH_QUAL = 41 + qual_offset;
    GZIP_FILE f = openCheckpoint(merged_reads_checkpoint, "w");
    char pair = '0', prev_pair = '0';
    while (get_next_fq_record(fqr, id_buf, seq_buf, quals_buf)) {
        pair = num_reads % 2 == 0 ? '1' : '2';
        if (hexify_libi >= 0) {
            hexifyId(getStartBuffer(id_buf), hexify_libi, (pair == '1' ? &read_id_1 : NULL), (pair == '2' ? &read_id_2 : NULL), THREADS);
        }
        DBG2("Merging %s\n", getStartBuffer(id_buf));
        num_reads++;
        if (files_per_pair == 2) {
            // read read2
            swapBuffer(prev_id_buf, id_buf);
            swapBuffer(prev_seq_buf, seq_buf);
            swapBuffer(prev_quals_buf, quals_buf);
            int readLen = getLengthBuffer(prev_seq_buf);
            if (maxReadLen < readLen) maxReadLen = readLen;
    
            if (!get_next_fq_record(fqr2, id_buf, seq_buf, quals_buf)) {
                DIE("Did not find matching paired read for lib=%s, num_reads=%lld hexid=%s\n", libName, (lld) num_reads, getStartBuffer(id_buf));
            }
            DBG2("with %s\n", getStartBuffer(id_buf));
            prev_pair = pair;
            pair = '2';
            if (hexify_libi >= 0) {
                hexifyId(getStartBuffer(id_buf), hexify_libi, NULL, &read_id_2, THREADS);
            }
            num_reads++;
        }
        char *id = getStartBuffer(id_buf);
        char *seq = getStartBuffer(seq_buf);
        char *quals = getStartBuffer(quals_buf);

        // ensure prev buffers are large enough for a merged read
        growBuffer(prev_seq_buf, getLengthBuffer(seq_buf));
        growBuffer(prev_quals_buf, getLengthBuffer(quals_buf));

        char *prev_id = getStartBuffer(prev_id_buf);
        char *prev_seq = getStartBuffer(prev_seq_buf);
        char *prev_quals = getStartBuffer(prev_quals_buf);
        int needs_pair_suffix = 1;
        if (strlen(prev_id)) {
            if (files_per_pair == 2 || (prev_pair == '1' && pair == '2')) {
                uint8_t namelen = strlen(id);
                uint8_t prev_namelen = strlen(prev_id);
                if (namelen != prev_namelen || strncmp(prev_id, id, namelen - 1) != 0) {
                    DIE("Error in paired reads file %s line %ld: names don't match, %s != %s\n",
                        libName, num_reads * 4, prev_id, id);
                }
                if (namelen >= 3 && id[namelen-2] == '/') needs_pair_suffix = 0;

                int8_t is_merged = 0;
                int8_t abortMerge = 0;

                // now make a copy and revcomp the second mate pair and reverse the second quals
                copyBuffer(tmpseq_buf, seq_buf);
                char *revcomp_seq = getStartBuffer(tmpseq_buf);
                switch_code(reverse(revcomp_seq));

                copyBuffer(tmpquals_buf, quals_buf);
                char *rev_quals = getStartBuffer(tmpquals_buf);
                reverse(rev_quals);

                // calcualte the sequenc length
                int prev_len = strlen(prev_seq), revcomp_len = strlen(revcomp_seq);
                if (prev_len > 32767 || revcomp_len > 37676) DIE("Can not merge ultra-long reads %lld or %lld size\n", prev_len, revcomp_len);

                // use start_i to offset inequal lengths which can be very different but still overlap near the end.  250 vs 178..
                int16_t len = (revcomp_len < prev_len) ? revcomp_len : prev_len;
                int16_t start_i = ((len == prev_len) ? 0 : prev_len - len);
                int16_t found_i = -1;
                int16_t best_i = -1;
                int16_t best_mm = len;
                double best_perror = -1.0;

                // slide along prev_seq
                for (int16_t i = 0; i < len - MIN_OVERLAP + EXTRA_TEST_OVERLAP; i++) { // test less overlap than MIN_OVERLAP
                    if (abortMerge) break;
                    int16_t overlap = len-i;
                    int16_t this_max_mismatch = MAX_MISMATCHES + (EXTRA_MISMATCHES_PER_1000 * overlap / 1000);
                    int16_t error_max_mismatch = this_max_mismatch * 4 / 3 + 1; // 33% higher
                    int fastMismatch = fastCountMismatches(prev_seq + start_i + i, revcomp_seq, overlap, error_max_mismatch);
                    if (fastMismatch > error_max_mismatch) {
                        continue;
                    }
                    int16_t matches = 0, mismatches = 0, bothNs = 0, Ncount = 0;
                    int16_t overlapChecked = 0;
                    double perror = 0.0;
                    for (int16_t j = 0; j < overlap; j++) {
                        overlapChecked++;
                        char ps = prev_seq[start_i + i + j];
                        char rs = revcomp_seq[j];
                        if (ps == rs) {
                            matches++;
                            if (ps == 'N') {
                                Ncount += 2;
                                if (bothNs++) {
                                    abortMerge++;
                                    num_ambiguous++;
                                    DBG2("Both have Ns at the same spot\n");
                                    break; // do not match multiple Ns in the same position -- 1 is okay
                                }
                            }
                        } else {
                            mismatches++;
                            if (ps == 'N') {
                                mismatches++; // N still counts as a mismatch
                                Ncount++;
                                prev_quals[start_i + i + j] = qual_offset;
                                assert(rev_quals[j] - qual_offset < Q2PerrorSize); 
                                assert(rev_quals[j] - qual_offset >= 0); 
                                perror += Q2Perror[rev_quals[j] - qual_offset];
                            } else if (rs == 'N') {
                                Ncount++;
                                mismatches++; // N still counts as a mismatch
                                rev_quals[j] = qual_offset;
                                assert(prev_quals[start_i + i + j] - qual_offset < Q2PerrorSize); 
                                assert(prev_quals[start_i + i + j] - qual_offset >= 0);
                                perror += Q2Perror[prev_quals[start_i + i + j] - qual_offset];
                            }
                            if (MAX_PERROR > 0.0) {

                                assert(prev_quals[start_i + i + j] >= qual_offset);
                                assert(rev_quals[j] >= qual_offset);
                                uint8_t q1 = prev_quals[start_i + i + j] - qual_offset;
                                uint8_t q2 = rev_quals[j] - qual_offset;
                                if (q1 < 0 || q2 < 0 || q1 >= Q2PerrorSize || q2 >= Q2PerrorSize) {
                                    DIE("Invalid quality score for read %s '%c' or %s '%c' assuming common qual_offset of %d. Please check your data and make sure it follows a single consistent quality scoring model (phred+64 vs. phred+33)\n", prev_id,  prev_quals[start_i + i + j], id, rev_quals[j], qual_offset);
                                }
                                // sum perror as the difference in q score perrors
                                uint8_t diffq = (q1 > q2) ? q1-q2 : q2-q1;
                                if (diffq <= 2) {
                                    perror += 0.5; // cap at flipping a coin when both quality scores are close
                                } else {
                                    assert(diffq < Q2PerrorSize);
                                    perror += Q2Perror[diffq];
                                }
                            }
                        }
                        if (Ncount > 3) {
                            abortMerge++;
                            num_ambiguous++;
                            break; // do not match reads with many Ns
                        }
                        if (mismatches > error_max_mismatch) break;
                    }
                    int16_t match_thres = overlap - this_max_mismatch;
                    if (match_thres < MIN_OVERLAP) match_thres = MIN_OVERLAP;
                    if (matches >= match_thres && overlapChecked == overlap && mismatches <= this_max_mismatch && perror/overlap <= MAX_PERROR) {
                        if (best_i < 0 && found_i < 0) {
                            DBG2("Found potential match at i=%d MM=%d PE=%f\n", i, mismatches, perror);
                            best_i = i;
                            best_mm = mismatches;
                            best_perror = perror;
                        } else {
                            // another good or ambiguous overlap detected
                            DBG2("Found ambiguous match at i=%d MM=%d PE=%f\n", i, mismatches, perror);
                            num_ambiguous++;
                            best_i = -1;
                            best_mm = len;
                            best_perror = -1.0;
                            break;
                        }
                    } else if (overlapChecked == overlap && mismatches <= error_max_mismatch && perror/overlap <= MAX_PERROR * 4 / 3) {
                        // lower threshold for detection of an ambigious overlap
                        found_i = i;
                        if (best_i >= 0) {
                            // ambiguous mapping found after a good one was
                            DBG2("Found ambiguous (poor) match at i=%d MM=%d PE=%f\n", i, mismatches, perror);
                            num_ambiguous++;
                            best_i = -1;
                            best_mm = len;
                            best_perror = -1.0;
                            break;
                        }
                    }
                }

                if (best_i >= 0 && abortMerge == 0) {
                    int16_t i = best_i;
                    int16_t overlap = len-i;
                    // pick the base with the highest quality score for the overlapped region
                    for(int16_t j = 0 ; j < overlap ; j++) {
                        if (prev_seq[start_i+i+j] == revcomp_seq[j]) {
                            // match boost quality up to the limit
                            uint16_t newQual = prev_quals[start_i+i+j] + rev_quals[j] - qual_offset;
                            prev_quals[start_i+i+j] = ((newQual > MAX_MATCH_QUAL) ? MAX_MATCH_QUAL : newQual);
                            assert( prev_quals[start_i+i+j] >= prev_quals[start_i+i+j] );
                            assert( prev_quals[start_i+i+j] >= rev_quals[j] );
                        } else {
                            uint8_t newQual;
                            if (prev_quals[start_i+i+j] < rev_quals[j]) {
                                // use rev base and discount quality
                                newQual = rev_quals[j] - prev_quals[start_i+i+j] + qual_offset;
                                prev_seq[start_i+i+j] = revcomp_seq[j];
                            } else {
                                // keep prev base, but still discount quality
                                newQual = prev_quals[start_i+i+j] - rev_quals[j] + qual_offset;
                            }
                            prev_quals[start_i+i+j] = ((newQual > (2+qual_offset)) ? newQual : (2+qual_offset)); // a bit better than random chance here
                        }
                        assert(prev_quals[start_i+i+j] >= qual_offset);
                    }

                    // include the remainder of the revcomp_seq and quals
                    strcpy(prev_seq + start_i + i + overlap, revcomp_seq + overlap);
                    strcpy(prev_quals + start_i + i + overlap, rev_quals + overlap);
                       
                    is_merged = 1;
                    num_merged++;
                    DBG2("Merged: %s OL=%d MM=%d PE=%0.2f\n", prev_id, overlap, best_mm, best_perror);
                    bytes_written += GZIP_PRINTF(f, "%s%s\n%s\n+\n%s\n", prev_id, (needs_pair_suffix ? "/1" : ""), prev_seq, prev_quals);
                    bytes_written += GZIP_PRINTF(f, "%s%s\nN\n+\n%c\n", id, (needs_pair_suffix ? "/2" : ""), qual_offset); // output an empty record, keeping a dummy entry for Read2
                    int readLen = strlen(prev_seq); // cacluate new merged length
                    if (maxReadLen < readLen) maxReadLen = readLen;
                    mergedLen += readLen;
                    overlapLen += overlap;
                }
                if (!is_merged) {
                    DBG2("NOT MERGED: %s abortMerge=%d\n", prev_id, abortMerge);
                    bytes_written += GZIP_PRINTF(f, "%s%s\n%s\n+\n%s\n", prev_id, (needs_pair_suffix ? "/1" : ""), prev_seq, prev_quals);
                    bytes_written += GZIP_PRINTF(f, "%s%s\n%s\n+\n%s\n", id, (needs_pair_suffix ? "/2" : ""), seq, quals);
                }
                num_read_pairs++;
                // ensure this read pair will not be used again
                resetBuffer(id_buf);
                resetBuffer(prev_id_buf);
                pair = '0';
                prev_pair = '0';
            } else {
                WARN("Mismatched pairs! '%s' '%s')\n", id, prev_id);
            }
        }
        if (maxReadLen < getLengthBuffer(seq_buf)) maxReadLen = getLengthBuffer(seq_buf);
        if (files_per_pair == 1) {
            swapBuffer(prev_id_buf, id_buf);
            swapBuffer(prev_seq_buf, seq_buf);
            swapBuffer(prev_quals_buf, quals_buf);
            prev_pair = pair;
        }
    }
    if (num_reads != num_read_pairs * 2) {
        DIE("Mismatch in number of reads (%ld) and num read pairs (%ld)\n", num_reads, num_read_pairs);
    }
    closeCheckpoint(f);
    // store the uncompressed size in a secondary file
    char uncompressed_fname[MAX_FILE_PATH+30];
    sprintf(uncompressed_fname, "%s.uncompressedSize", merged_reads_checkpoint);
    FILE *fsz = openCheckpoint0(uncompressed_fname, "w");
    fwrite_chk(&bytes_written, sizeof(int64_t), 1, fsz);
    closeCheckpoint0(fsz);

    // cleanup
    freeBuffer(id_buf);
    freeBuffer(seq_buf);
    freeBuffer(quals_buf);
    freeBuffer(prev_id_buf);
    freeBuffer(prev_seq_buf);
    freeBuffer(prev_quals_buf);
    freeBuffer(tmpseq_buf);
    freeBuffer(tmpquals_buf);

    // update the maximum read length
    int max_read_len = reduce_long(maxReadLen, UPC_MAX, ALL_DEST);

    if (!MYTHREAD) {
        serial_printf("Writing maxReadLen for %s: %d\n", merged_reads_checkpoint, max_read_len);
        char buf[MAX_FILE_PATH+30];
        snprintf(buf, MAX_FILE_PATH+30, "%s.maxReadLen.txt", merged_reads_checkpoint);
        put_num_in_file(NULL, buf, max_read_len);
        snprintf(buf, MAX_FILE_PATH+30, "%s-readlen.txt", libName);
        FILE * f = fopen_chk(buf, "a");
        fprintf(f, "%d", max_read_len);
        fclose_track(f);
    }

    UPC_TIMED_BARRIER;
    long all_num_read_pairs = reduce_long(num_read_pairs, UPC_ADD, SINGLE_DEST);
    long all_num_merged = reduce_long(num_merged, UPC_ADD, SINGLE_DEST);
    long all_num_ambiguous = reduce_long(num_ambiguous, UPC_ADD, SINGLE_DEST);
    double avg_merged_len = reduce_long(mergedLen, UPC_ADD, SINGLE_DEST) / (double) ((all_num_merged > 0) ? all_num_merged : 1);
    double avg_overlap_len = reduce_long(overlapLen, UPC_ADD, SINGLE_DEST) / (double) ((all_num_merged > 0) ? all_num_merged : 1);
    stop_timer(T_MERGE_READS);
    serial_printf("Merged %ld out of %ld read pairs (%.2f%%, %0.2f%% ambiguous) with %0.1f avgLength (%d max), %0.1f avgOverlap in %.2f s\n", all_num_merged, all_num_read_pairs,
                  100.0 * all_num_merged / all_num_read_pairs, 100.0* all_num_ambiguous / all_num_read_pairs, avg_merged_len, max_read_len, avg_overlap_len,
                  get_elapsed_time(T_MERGE_READS));
}

static void merge_reads_lib(char * libName, int files_per_pair, int hexify_libi, int cached_io, int cached_reads, int qual_offset)
{
    char *base_dir = BASE_DIR(cached_io);
    if (files_per_pair == 0) {
        DIE("merge_reads is called on a single library!\n");
        return;
    }

    serial_printf("Merging reads from %s\n", libName);
    UPC_LOGGED_BARRIER;

    // read the fofn, merge each, then replace fofn
    Buffer newFofnBuffer = initBuffer(MAX_FILE_PATH);
    char fofn[MAX_FILE_PATH+10];
    sprintf(fofn, "%s.fofn", libName);
    Buffer fofnBuffer = broadcast_file(fofn);
    char file1[MAX_FILE_PATH], file2[MAX_FILE_PATH];
    while ( getsBuffer(fofnBuffer, file1, MAX_FILE_PATH) ) {
        file1[strlen(file1)-1] = '\0';
        if (files_per_pair == 2) {
            if ( !getsBuffer(fofnBuffer, file2, MAX_FILE_PATH) ) { DIE("Missing second file for paired libs: %s\n", file1); }
            file2[strlen(file2)-1] = '\0';
        }

        char *base = strdup(file1);
        char *reads_fname_base = basename(base);
        char *isgz = strstr(reads_fname_base, ".gz");
        if (isgz) isgz[0] = '\0';
        char merged_reads_checkpoint[MAX_FILE_PATH];
        sprintf(merged_reads_checkpoint, "%s-merged.fastq" GZIP_EXT, reads_fname_base);
        free(base);

        printfBuffer(newFofnBuffer, "%s\n", merged_reads_checkpoint);

        if (doesCheckpointExist(merged_reads_checkpoint)) {
            serial_printf("Skipping merge of %s as they already exist\n", merged_reads_checkpoint);
            if (!doesLocalCheckpointExist(merged_reads_checkpoint)) {
                restoreLocalCheckpoint(merged_reads_checkpoint);
                strcat(merged_reads_checkpoint, ".uncompressedSize");
                if ( !doesLocalCheckpointExist(merged_reads_checkpoint) && doesGlobalCheckpointExist(merged_reads_checkpoint) ) {
                    restoreCheckpoint(merged_reads_checkpoint);
                }
            }
            continue;
        }
        serial_printf("Merging reads from %s %s into %s\n", file1, (files_per_pair==2 ? file2 : ""), merged_reads_checkpoint);
    
        fq_reader_t fqr = create_fq_reader(), fqr2;
        open_fq(fqr, file1, cached_reads, base_dir, cached_reads ? -1 : broadcast_file_size(file1));
        if (files_per_pair == 2) {
            fqr2 = create_fq_reader();
            open_fq(fqr2, file2, cached_reads, base_dir, cached_reads ? -1 : broadcast_file_size(file2));
        }

        merge_reads_file(fqr, fqr2, merged_reads_checkpoint, libName, files_per_pair, hexify_libi, qual_offset);
     
        destroy_fq_reader(fqr);
        if (files_per_pair == 2) {
            destroy_fq_reader(fqr2);
        }
    }
    UPC_LOGGED_BARRIER;
    if (!MYTHREAD) {
        serial_printf("Replacing fofn: %s\n", fofn);
        char new[MAX_FILE_PATH+30];
        sprintf(new, "%s-PREMERGE", fofn);
        rename(fofn, new);
        FILE * f = fopen_chk(fofn, "w");
        fwrite_chk(getStartBuffer(newFofnBuffer), 1, getLengthBuffer(newFofnBuffer), f);
        fclose_track(f);
    }
    UPC_LOGGED_BARRIER;
}

// the stage by lib
int merge_reads_by_lib(int argc, char **argv) {
    if (argc != 13) DIE("Invalid call: %d\n", argc);
    if (strcmp(argv[1], "-l") != 0) DIE("expected -l\n");
    char * libName = argv[2];
    if (strcmp(argv[3], "-p") != 0) DIE("expected -p\n");
    int files_per_pair = atoi(argv[4]);
    if (strcmp(argv[5], "-h") != 0) DIE("expected -h\n");
    int hexify_libi = atoi(argv[6]);
    if (strcmp(argv[7], "-c") != 0) DIE("expected -c\n");
    int cached_io = atoi(argv[8]);
    if (strcmp(argv[9], "-X") != 0) DIE("expected -X\n");
    int cached_reads = atoi(argv[10]);
    if (strcmp(argv[11], "-q") != 0) DIE("expected -q\n");
    int qual_offset = atoi(argv[12]);

    merge_reads_lib(libName, files_per_pair, hexify_libi, cached_io, cached_reads, qual_offset);
    return 0;
}

// merges all eligible libs, and updates lib values
void merge_reads(cfg_t *cfg, const char * all_inputs_fofn_name)
{
    if (cfg->merge_reads == 0) DIE("merge_reads should not be called when the merge_reads option == 0\n");
    serial_printf("Starting to merge overlapping reads\n");
    int didRun = 0, numToRun = 0;
    for (int libi = 0; libi < cfg->num_libs; libi++) {
        lib_t *lib = cfg->libs + libi;

        int canMergeLib = (lib->files_per_pair > 0  // paired or interleaved libs
                           && (lib->ins_avg == 0 || lib->ins_avg < ((lib->read_len?lib->read_len:250) - 12) * 2 + 3*lib->std_dev || lib->innie));
        
        if (canMergeLib) {
            serial_printf("Expecting many reads can be merged\n");
        }
        // but always merge a paired library so it becomes interleaved
        if ( canMergeLib || lib->files_per_pair == 2 ) {
            // merge reads into per_thread files
            numToRun++;
            char stageLabel[MAX_FILE_PATH+20];
            sprintf(stageLabel, "merge_reads-%s", lib->name);
            didRun += exec_stage(cfg->stages, NOT_SCAFF, merge_reads_by_lib, stageLabel,
                                 "-l %s", lib->name,
                                 "-p %d", lib->files_per_pair,
                                 "-h %d", libi,
                                 "-c %d", cfg->cached_io,
                                 "-X %d", (cfg->merge_reads == 1 ? 0 : cfg->cached_reads),
                                 "-q %d", cfg->qual_offset,
                                 NULL);

            lib->files_per_pair = 1; // paired libs are now interleaved
            lib->is_merged = 1;
            // get the new maxReadLen
            if (!MYTHREAD) {
                char buf[MAX_FILE_PATH+20];
                sprintf(buf, "%s.fofn", lib->name);
                FILE *f = fopen_chk(buf, "r");
                while( fgets(buf, MAX_FILE_PATH, f) ) {
                    buf[strlen(buf)-1] = '\0';
                    if (strchr(buf, '/')) DIE("Invalid entry in %s.fofn, which should not contain a slash after merging! %s\n", lib->name, buf);
                    char buf2[MAX_FILE_PATH+40];
                    sprintf(buf2, "%s.maxReadLen.txt", buf);
                    int maxlen = get_num_from_file(NULL, buf2); 
                    if (lib->read_len < maxlen) lib->read_len = maxlen;
                }
                fclose_track(f);
            }
            lib->read_len = broadcast_int(lib->read_len, 0);
            serial_printf("Updated max_read len for %s: %d\n", lib->name, lib->read_len);
        } else {
            serial_printf("Skipping library %s fpp=%d insert=%d stddev=%d read_len=%d that will not have any overlaps\n", lib->name, lib->files_per_pair, lib->ins_avg, lib->std_dev, lib->read_len);
            serial_printf("Loading %s to local checkpoints\n", lib->name);
            loadfq_library(lib->name, BASE_DIR(cfg->cached_io));
        }
    }
    if ((!MYTHREAD) && (numToRun == 0 || didRun > 0)) { // if any ran, then this would be the last restart to succeed
        rebuild_all_inputs_fofn(all_inputs_fofn_name, cfg, "MERGED_READS", numToRun > 0);
    }
    cfg->cached_reads = 1; // now all input files are in cached io or are checkpoints
}
