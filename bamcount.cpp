
#define __STDC_FORMAT_MACROS
#include <iostream>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <cerrno>
#include <cstring>
#include <unordered_map>
#include <vector>
#include <thread>
#include <math.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <sys/stat.h>

typedef std::vector<std::string> strvec;

static const void print_version() {
    //fprintf(stderr, "bamcount %s\n", string(BAMCOUNT_VERSION).c_str());
    std::cout << "bamcount " << std::string(BAMCOUNT_VERSION) << std::endl;
}

static const char USAGE[] = "BAM Multi-Threaded Reading Utility\n"
    "\n"
    "Usage:\n"
    "  bamcount <bam|bw|-> [options]\n"
    "\n"
    "Options:\n"
    "  -h --help            Show this screen.\n"
    "  --version            Show version.\n"
    "  --threads            # of threads to do: BAM decompression OR compute sums over multiple BigWigs in parallel\n"
    "                       if the 2nd is intended then a TXT file listing the paths to the BigWigs to process in parallel\n"
    "                       should be passed in as the main input file instead of a single BigWig file (EXPERIMENTAL).\n"
    "\n";

static const char* get_positional_n(const char ** begin, const char ** end, size_t n) {
    size_t i = 0;
    for(const char **itr = begin; itr != end; itr++) {
        if((*itr)[0] != '-' || strlen(*itr) == 1) {
            if(i++ == n) {
                return *itr;
            }
        }
    }
    return nullptr;
}

static bool has_option(const char** begin, const char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

/**
 * Return the argument after the given one, (or further downstream when shift > 0).
 */
static const char** get_option(
        const char** begin,
        const char** end,
        const std::string& option,
        unsigned shift = 0)
{
    const char** itr = std::find(begin, end, option);
    return itr + shift + 1;
}

void process_bam_worker(const int worker_idx, const char* bam_arg, char** regions, uint64_t region_offset, uint64_t num_regions) {
    htsFile *bam_fh = sam_open(bam_arg, "r");
    if(!bam_fh) {
        std::cerr << "ERROR: Could not open " << bam_arg << ": "
                  << std::strerror(errno) << std::endl;
        return;
    }
    bam_hdr_t* hdr = sam_hdr_read(bam_fh);
    hts_idx_t* idx = nullptr;
    idx = sam_index_load(bam_fh, bam_arg);
    if (idx == 0) { // index is unavailable
        fprintf(stderr, "[main_samview] random alignment retrieval only works for indexed BAM or CRAM files.\n");
        return;
    }
    bam1_t *rec = bam_init1();
    if(!rec) {
        std::cerr << "ERROR: Could not initialize BAM object: "
                  << std::strerror(errno) << std::endl;
        return;
    }
    uint64_t recs = 0;
    int result;
    uint64_t max_region_offset = region_offset + num_regions;
    for(int z = region_offset; z < max_region_offset; z++) {
        hts_itr_t *iter = sam_itr_querys(idx, hdr, regions[z]); // parse a region in the format like `chr2:100-200'
        while ((result = sam_itr_next(bam_fh, iter, rec)) >= 0) {
            recs++;
        }
        hts_itr_destroy(iter);
    }
    bam_destroy1(rec);
    hts_idx_destroy(idx); // destroy the BAM index
    fprintf(stdout,"Worker %d read %lu records\n",worker_idx,recs);
}

int go(const char* bam_arg, int argc, const char** argv, htsFile *bam_fh, bool is_bam) {

    int nthreads = 8;
    if(has_option(argv, argv+argc, "--threads")) {
        const char** nthreads_ = get_option(argv, argv+argc, "--threads");
        nthreads = atoi(*nthreads_);
    }
    bam_hdr_t *hdr = sam_hdr_read(bam_fh);
    if(!hdr) {
        std::cerr << "ERROR: Could not read header for " << bam_arg
                  << ": " << std::strerror(errno) << std::endl;
        return 1;
    }
    //find total bases to split up among threads
    int64_t total_bases = 0;
    std::vector<int> chr_idxs;
    for(int i = 0; i < hdr->n_targets; i++) {
        total_bases += hdr->target_len[i];
        chr_idxs.push_back(i);
    }

    char** regions = new char*[hdr->n_targets];
    //ds to hold 1 or more regions per thread, so an array of vectors
    //strvec* regions = new strvec[nthreads];
    //int regions_per_thread = hdr->n_targets / nthreads;
    uint64_t bases_per_thread = total_bases / nthreads;
    printf("bases per_thread %lu\n",bases_per_thread);
    //now sort list of chromosomes by bases ascending
    auto sortClosure = [&](int x, int y)
    {
            return hdr->target_len[x] < hdr->target_len[y];
    };
    std::sort(chr_idxs.begin(), chr_idxs.end(), sortClosure);
    /*for(int i = 0; i < hdr->n_targets; i++)
        printf("chr idx %d\n",chr_idxs[i]);*/
    for(int i = 0; i < hdr->n_targets; i++) {
        regions[i] = new char[100];
        sprintf(regions[i],"%s:1-%d",hdr->target_name[i],hdr->target_len[i]);
    }
        
    /*regions[0] = new char[50];
    regions[0] = "chr1:1-260000000";
    regions[3] = "chr4:1-260000000";*/
    std::vector<std::thread> threads;
    uint64_t offset = 0;
    uint64_t num_regions = 0;
    uint64_t cur_bytes = 0;
    int j = 0;
    int used_threads = 0;
    /*for(int i=0; i < nthreads && j < hdr->n_targets; i++) {
        while((cur_bytes <= bases_per_thread || i == (nthreads - 1)) && j < hdr->n_targets) {
            cur_bytes += hdr->target_len[chr_idxs[j++]];
            num_regions++;
        }
        if(cur_bytes == 0)
            break;
        threads.push_back(std::thread(process_bam_worker, i, std::ref(bam_arg), std::ref(regions), offset, num_regions));
        used_threads++;
        cur_bytes = 0;
        num_regions = 0;
        offset += num_regions;
    }
    for(int i=0; i < used_threads; i++)
        threads[i].join();*/

    while(j < hdr->n_targets) {
        for(int i=0; i < nthreads && j+i < hdr->n_targets; i++)
            //threads.push_back(std::thread(process_bam_worker, i, std::ref(bam_arg), std::ref(regions[j+i])));
            threads.push_back(std::thread(process_bam_worker, i, std::ref(bam_arg), std::ref(regions), j+i, 1));
        int z = threads.size();
        for(int i=0; i < z; i++)
            threads[i].join();
        threads.clear();
        j += z;
    }
    //uint64_t recs = 0;
    //fprintf(stdout,"Parent read %lu records\n",recs);
    fprintf(stdout,"Parent done\n");
    return 0;
}

int main(int argc, const char** argv) {
    argv++; argc--;  // skip binary name
    if(argc == 0 || has_option(argv, argv + argc, "--help") || has_option(argv, argv + argc, "--usage")) {
        print_version();
        std::cout << std::endl << USAGE << std::endl;
        return 0;
    }
    if(has_option(argv, argv + argc, "--version")) {
        print_version();
        return 0;
    }
    const char *bam_arg = get_positional_n(argv, argv+argc, 0);
    if(!bam_arg) {
        std::cerr << "ERROR: Could not find <bam|bw> positional arg" << std::endl;
        return 1;
    }
    std::cerr << "Processing \"" << bam_arg << "\"" << std::endl;

    htsFile *bam_fh = sam_open(bam_arg, "r");
    if(!bam_fh) {
        std::cerr << "ERROR: Could not open " << bam_arg << ": "
                  << std::strerror(errno) << std::endl;
        return 1;
    }
    const htsFormat* format = hts_get_format(bam_fh);
    const char* hts_format_ex = hts_format_file_extension(format);
    bool is_bam = true;
    std::ios::sync_with_stdio(false);
    return go(bam_arg, argc, argv, bam_fh, is_bam);
}
