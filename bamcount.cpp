
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

void process_bam_worker(const int worker_idx, const char* bam_arg, bam_hdr_t *hdr_) {
    htsFile *bam_fh = sam_open(bam_arg, "r");
    if(!bam_fh) {
        std::cerr << "ERROR: Could not open " << bam_arg << ": "
                  << std::strerror(errno) << std::endl;
        return;
    }
    //hts_set_threads(bam_fh, 1);
    bam1_t *rec = bam_init1();
    if(!rec) {
        std::cerr << "ERROR: Could not initialize BAM object: "
                  << std::strerror(errno) << std::endl;
        return;
    }
    bam_hdr_t *hdr = sam_hdr_read(bam_fh);
    uint64_t recs = 0;
    while(sam_read1(bam_fh, hdr, rec) >= 0) {
        recs++;
        bam1_core_t *c = &rec->core;
        //read name
        char* qname = bam_get_qname(rec);
    }
    fprintf(stdout,"Worker %d read %lu records\n",worker_idx,recs);
}

int go(const char* bam_arg, int argc, const char** argv, htsFile *bam_fh, bool is_bam) {
    int nthreads = 4;
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
    std::vector<std::thread> threads;
    for(int i=0; i < nthreads; i++) 
            threads.push_back(std::thread(process_bam_worker, i, std::ref(bam_arg), std::ref(hdr)));
    for(int i=0; i < threads.size(); i++)
        threads[i].join();
    uint64_t recs = 0;
    fprintf(stdout,"Parent read %lu records\n",recs);
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
