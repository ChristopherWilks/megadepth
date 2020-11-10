/* The MIT License

    Copyright (c) 2018-  by Christopher Wilks <cwilks3@jhu.edu> and Ben Langmead
                         <langmea@cs.jhu.edu>

    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

*/


#define __STDC_FORMAT_MACROS
#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <thread>
#include <iterator>
#include <numeric>

#include <zlib.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <sys/stat.h>
#include "bigWig.h"
#include "countlut.hpp"
#ifdef WINDOWS_MINGW
    #include <unordered_map>
    #include <unordered_set>
    #include "getline.h"
    #include "mingw-std-threads/mingw.thread.h"
    template<class K, class V>
    using hashmap = std::unordered_map<K, V>;
    template<class V2>
    using hashset = std::unordered_set<V2>;
#else
    #include "robin_hood.h"
    template<class K, class V>
    using hashmap = robin_hood::unordered_map<K, V>;
    template<class V2>
    using hashset = robin_hood::unordered_set<V2>;
#endif
#if defined(__AVX2__) || defined(__SSE2__)
#include <x86intrin.h>
#endif
#if defined(__GNUC__) || defined(__clang__)
#  ifndef unlikely
#    define unlikely(x) __builtin_expect(!!(x), 0)
#  endif
#  ifndef likely
#    define likely(x) __builtin_expect(!!(x), 1)
#  endif
#endif

int UNKNOWN_FORMAT=-1;
int BAM_FORMAT = 1;
int BW_FORMAT = 2;
int CRAM_FORMAT = 3;

//taken from HTSlib bgzip
int BGZF_WRITE_WINDOW_SZ = 64 * 1024;

//critical to use a high value here for remote BigWigs
//accesses, has much less (maybe no) effect on local processing
const uint32_t default_BW_READ_BUFFER = 1<<30;
uint32_t BW_READ_BUFFER = default_BW_READ_BUFFER;

bool SUMS_ONLY = false;

typedef std::vector<std::string> strvec;
typedef hashmap<std::string, uint64_t> mate2len;
typedef hashmap<std::string, double*> str2dblist;

uint64_t MAX_INT = (2^63);
//how many intervals to start with for a chromosome in a BigWig file
//uint64_t STARTING_NUM_INTERVALS = 1000;
uint64_t STARTING_NUM_INTERVALS = 1000000;
//used for --annotation where we read a 3+ column BED file
static const int CHRM_COL=0;
static const int START_COL=1;
static const int END_COL=2;
//1MB per line should be more than enough for CIO
static const int LINE_BUFFER_LENGTH=1048576;
static const int BIGWIG_INIT_VAL = 17;
static double SOFTCLIP_POLYA_TOTAL_COUNT_MIN=3;
static double SOFTCLIP_POLYA_RATIO_MIN=0.8;

//used for buffering up text/gz output
static const int OUT_BUFF_SZ=4000000;
static const int COORD_STR_LEN=34;

enum Op { csum, cmean, cmin, cmax };

static const void print_version() {
    std::cout << "megadepth " << std::string(MEGADEPTH_VERSION) << std::endl;
}

static const char USAGE[] = "BAM and BigWig utility.\n"
    "\n"
    "Usage:\n"
    "  megadepth <bam|bw|-> [options]\n"
    "\n"
    "Options:\n"
    "  -h --help                Show this screen.\n"
    "  --version                Show version.\n"
    "  --threads                # of threads to do: BAM decompression OR compute sums over multiple BigWigs in parallel\n"
    "                            if the 2nd is intended then a TXT file listing the paths to the BigWigs to process in parallel\n"
    "                            should be passed in as the main input file instead of a single BigWig file (EXPERIMENTAL).\n"
    "  --prefix                 String to use to prefix all output files.\n"
    "  --no-auc-stdout          Force all AUC(s) to be written to <prefix>.auc.tsv rather than STDOUT\n"
    "  --no-annotation-stdout   Force summarized annotation regions to be written to <prefix>.annotation.tsv rather than STDOUT\n"
    "  --no-coverage-stdout     Force covered regions to be written to <prefix>.coverage.tsv rather than STDOUT\n"
    "  --keep-order             Output annotation coverage in the order chromosomes appear in the BAM/BigWig file\n"
    "                           The default is to output annotation coverage in the order chromosomes appear in the annotation BED file.\n"
    "                           This is only applicable if --annotation is used for either BAM or BigWig input.\n"
    "\n"
    "BigWig Input:\n"
    "Extract regions and their counts from a BigWig outputting BED format if a BigWig file is detected as input (exclusive of the other BAM modes):\n"
    "                                          Extracts all reads from the passed in BigWig and output as BED format.\n"
    "                                           This will also report the AUC over the annotated regions to STDOUT.\n"
    "                                           If only the name of the BigWig file is passed in with no other args, it will *only* report total AUC to STDOUT.\n"
    "  --annotation <bed>                      Only output the regions in this BED applying the argument to --op to them.\n"
    "  --op <sum[default], mean, min, max>     Statistic to run on the intervals provided by --annotation\n"
    "  --sums-only                             Discard coordinates from output of summarized regions\n"
    "  --bwbuffer <1GB[default]>               Size of buffer for reading BigWig files, critical to use a large value (~1GB) for remote BigWigs.\n"
    "                                           Default setting should be fine for most uses, but raise if very slow on a remote BigWig.\n"
    "\n"
    "\n"
    "BAM Input:\n"
    "Extract basic junction information from the BAM, including co-occurrence\n"
    "If only the name of the BAM file is passed in with no other args, it will *only* report total AUC to STDOUT.\n"
    "  --fasta	            Path to the reference FASTA file if a CRAM file is passed as the input file (ignored otherwise)\n"
    "                       If not passed, references will be downloaded using the CRAM header.\n"
    "  --junctions          Extract jx coordinates, strand, and anchor length, per read\n"
    "                       writes to a TSV file <prefix>.jxs.tsv\n"
    "  --longreads          Modifies certain buffer sizes to accommodate longer reads such as PB/Oxford.\n"
    "  --filter-in          Integer bitmask, any bits of which alignments need to have to be kept (similar to samtools view -f).\n"
    "  --filter-out         Integer bitmask, any bits of which alignments need to have to be skipped (similar to samtools view -F).\n"
    "\n"
    "Non-reference summaries:\n"
    "  --alts                       Print differing from ref per-base coverages\n"
    "                               Writes to a CSV file <prefix>.alts.tsv\n"
    "  --include-softclip           Print a record to the alts CSV for soft-clipped bases\n"
    "                               Writes total counts to a separate TSV file <prefix>.softclip.tsv\n"
    "  --only-polya                 If --include-softclip, only print softclips which are mostly A's or T's\n"
    "  --include-n                  Print mismatch records when mismatched read base is N\n"
    "  --print-qual                 Print quality values for mismatched bases\n"
    "  --delta                      Print POS field as +/- delta from previous\n"
    "  --require-mdz                Quit with error unless MD:Z field exists everywhere it's\n"
    "                               expected\n"
    "  --head                       Print sequence names and lengths in SAM/BAM header\n"
    "\n"
    "Coverage and quantification:\n"
    "  --coverage           Print per-base coverage (slow but totally worth it)\n"
    "  --auc                Print per-base area-under-coverage, will generate it for the genome\n"
    "                       and for the annotation if --annotation is also passed in\n"
    "                       Defaults to STDOUT, unless other params are passed in as well, then\n"
    "                       if writes to a TSV file <prefix>.auc.tsv\n"
    "  --bigwig             Output coverage as BigWig file(s).  Writes to <prefix>.bw\n"
    "                       (also <prefix>.unique.bw when --min-unique-qual is specified).\n"
    "                       Requires libBigWig.\n"
    "  --annotation <bed>   Path to BED file containing list of regions to sum coverage over\n"
    "                       (tab-delimited: chrm,start,end)\n"
    "  --op <sum[default], mean>     Statistic to run on the intervals provided by --annotation\n"
    "  --no-index           If using --annotation, skip the use of the BAM index (BAI) for pulling out regions.\n"
    "                       Setting this can be faster if doing windows across the whole genome.\n"
    "  --min-unique-qual <int>\n"
    "                       Output second bigWig consisting built only from alignments\n"
    "                       with at least this mapping quality.  --bigwig must be specified.\n"
    "                       Also produces second set of annotation sums based on this coverage\n"
    "                       if --annotation is enabled\n"
    "  --double-count       Allow overlapping ends of PE read to count twice toward\n"
    "                       coverage\n"
    "  --num-bases          Report total sum of bases in alignments processed (that pass filters)\n"
    "  --gzip               Turns on gzipping of coverage output (no effect if --bigwig is passsed),\n"
    "                       this will also enable --no-coverage-stdout.\n"
    "\n"
    "Other outputs:\n"
    "  --read-ends          Print counts of read starts/ends, if --min-unique-qual is set\n"
    "                       then only the alignments that pass that filter will be counted here\n"
    "                       Writes to 2 TSV files: <prefix>.starts.tsv, <prefix>.ends.tsv\n"
    "  --frag-dist          Print fragment length distribution across the genome\n"
    "                       Writes to a TSV file <prefix>.frags.tsv\n"
    "  --echo-sam           Print a SAM record for each aligned read\n"
    "  --ends               Report end coordinate for each read (useful for debugging)\n"
    "  --test-polya         Lower Poly-A filter minimums for testing (only useful for debugging/testing)\n"
    "\n";

int my_write(void* fh, char* buf, uint32_t buf_len) {
#if USE_POSIX
    return ::write(::fileno(fh), buf, bu_len);
#else
    return std::fwrite(buf, 1, buf_len, (FILE *)fh);
#endif
}

int my_gzwrite(void* fh, char* buf, uint32_t buf_len) {
    return bgzf_write((BGZF*)fh, buf, buf_len);
    //return gzwrite(*((gzFile*) fh), buf, buf_len);
}

template <typename T>
int print_local(char* buf,const char* c, long start, long end, T val, double* local_vals, long z);

template <typename T>
int print_local_sums_only(char* buf,const char* c, long start, long end, T val, double* local_vals, long z);

template <typename T>
int print_shared(char* buf,const char* c, long start, long end, T val, double* local_vals, long z);

template <typename T>
int print_shared_sums_only(char* buf,const char* c, long start, long end, T val, double* local_vals, long z);

template <>
int print_local<long>(char* buf,const char* c, long start, long end, long val, double* local_vals, long z) {
        return sprintf(buf, "%s\t%lu\t%lu\t%lu\n", c, start, end, (long) local_vals[z]);
}

template <>
int print_local_sums_only<long>(char* buf,const char* c, long start, long end, long val, double* local_vals, long z) {
        return sprintf(buf, "%lu\n", (long) local_vals[z]);
}

template <>
int print_shared<long>(char* buf,const char* c, long start, long end, long val, double* local_vals, long z) {
        return sprintf(buf, "%s\t%lu\t%lu\t%lu\n", c, start, end, val);
}

template <>
int print_shared_sums_only<long>(char* buf,const char* c, long start, long end, long val, double* local_vals, long z) {
        return sprintf(buf, "%lu\n", val);
}

template <>
int print_shared<double>(char* buf, const char* c, long start, long end, double val, double* local_vals, long z) {
        return sprintf(buf, "%s\t%lu\t%lu\t%.2f\n", c, (long) start, (long) end, val);
}

template <>
int print_shared_sums_only<double>(char* buf, const char* c, long start, long end, double val, double* local_vals, long z) {
        return sprintf(buf, "%.2f\n", val);
}

template <>
int print_local<double>(char* buf, const char* c, long start, long end, double val, double* local_vals, long z) {
        return sprintf(buf, "%s\t%lu\t%lu\t%.2f\n", c, (long) start, (long) end, local_vals[z]);
}

template <>
int print_local_sums_only<double>(char* buf, const char* c, long start, long end, double val, double* local_vals, long z) {
        return sprintf(buf, "%.2f\n", local_vals[z]);
}

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

/**
 * Holds an MDZ "operation"
 * op can be
 */
struct MdzOp {
    char op;
    int run;
    char str[1024];
};

//from https://github.com/samtools/htslib/blob/7c04ea5c328547e9e8a9af4b932b87a3cb1939e6/hts.c#L82
int A_idx = 1;
int T_idx = 8;
static inline int polya_check(const uint8_t *str, size_t off, size_t run, char *c) {
    char seq_nt16_str_counts[16] = {0};
    for(size_t i = off; i < off + run; i++)
        seq_nt16_str_counts[bam_seqi(str, i)]++;
    int count = -1;
    if((seq_nt16_str_counts[A_idx] / (double) run) >= SOFTCLIP_POLYA_RATIO_MIN) {
        *c = 'A';
        count = seq_nt16_str_counts[A_idx];
    }
    else if((seq_nt16_str_counts[T_idx] / (double) run) >= SOFTCLIP_POLYA_RATIO_MIN) {
        *c = 'T';
        count = seq_nt16_str_counts[T_idx];
    }
    return count;
}

static const char seq_rev_nt16_str[] = "=TGMCRSVAWYHKDBN";
static inline std::ostream& seq_substring(std::ostream& os, const uint8_t *str, size_t off, size_t run, bool reverse=false) {
    if(reverse) {
        int i=(off+run)-1;
        while(((int) off) <= i) {
            int io = bam_seqi(str, i);
            os << seq_rev_nt16_str[io];
            i--;
        }
        return os;
    }
    for(size_t i = off; i < off + run; i++) {
        os << seq_nt16_str[bam_seqi(str, i)];
    }
    return os;
}

static inline std::ostream& kstring_out(std::ostream& os, const kstring_t *str) {
    for(size_t i = 0; i < str->l; i++) {
        os << str->s[i];
    }
    return os;
}

static inline std::ostream& cstr_substring(std::ostream& os, const uint8_t *str, size_t off, size_t run) {
    for(size_t i = off; i < off + run; i++) {
        os << (char)str[i];
    }
    return os;
}

static inline std::ostream& qstr_substring(std::ostream& os, const uint8_t *str, size_t off, size_t run, bool reverse=false) {
    if(reverse) {
        int i=(off+run)-1;
        while(((int) off) <= i) {
            os << (char)(str[i]+33);
            i--;
        }
        return os;
    }
    for(size_t i = off; i < off + run; i++) {
        os << (char)(str[i]+33);
    }
    return os;
}

/**
 * Parse given MD:Z extra field into a vector of MD:Z operations.
 */
static void parse_mdz(
        const uint8_t *mdz,
        std::vector<MdzOp>& ops)
{
    int i = 0;
    size_t mdz_len = strlen((char *)mdz);
    while(i < mdz_len) {
        if(isdigit(mdz[i])) {
            int run = 0;
            while(i < mdz_len && isdigit(mdz[i])) {
                run *= 10;
                run += (int)(mdz[i] - '0');
                i++;
            }
            if(run > 0) {
                ops.emplace_back(MdzOp{'=', run, ""});
                ops.back().str[0] = '\0';
            }
        } else if(isalpha(mdz[i])) {
            int st = i;
            while(i < mdz_len && isalpha(mdz[i])) i++;
            assert(i > st);
            ops.emplace_back(MdzOp{'X', i - st, ""});
            for(int j = 0; j < i ; j++) {
                ops.back().str[j] = mdz[st + j];
            }
            std::memcpy(ops.back().str, mdz + st, (size_t)(i - st));
            ops.back().str[i - st] = '\0';
        } else if(mdz[i] == '^') {
            i++;
            int st = i;
            while (i < mdz_len && isalpha(mdz[i])) i++;
            assert(i > st);
            ops.emplace_back(MdzOp{'^', i - st, ""});
            std::memcpy(ops.back().str, mdz + st, (size_t)(i - st));
            ops.back().str[i - st] = '\0';
        } else {
            std::stringstream ss;
            ss << "Unknown MD:Z operation: \"" << mdz[i] << "\"";
            throw std::runtime_error(ss.str());
        }
    }
}

static void output_from_cigar_mdz(
        const bam1_t *rec,
        std::vector<MdzOp>& mdz,
        std::fstream& fout,
        uint64_t* total_softclip_count,
        bool print_qual = false,
        bool include_sc = false,
        bool only_polya_sc = false,
        bool include_n_mms = false,
        bool delta = false)
{
    uint8_t *seq = bam_get_seq(rec);
    uint8_t *qual = bam_get_qual(rec);
    // If QUAL field is *. this array is just a bunch of 255s
    uint32_t *cigar = bam_get_cigar(rec);
    size_t mdzi = 0, seq_off = 0;
    int32_t ref_off = rec->core.pos;
    for(unsigned int k = 0; k < rec->core.n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        int run = bam_cigar_oplen(cigar[k]);
        if((strchr("DNMX=", BAM_CIGAR_STR[op]) != nullptr) && mdzi >= mdz.size()) {
            std::stringstream ss;
            ss << "Found read-consuming CIGAR op after MD:Z had been exhausted" << std::endl;
            throw std::runtime_error(ss.str());
        }
        if(op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) {
            // Look for block matches and mismatches in MD:Z string
            int runleft = run;
            while(runleft > 0 && mdzi < mdz.size()) {
                int run_comb = std::min(runleft, mdz[mdzi].run);
                runleft -= run_comb;
                assert(mdz[mdzi].op == 'X' || mdz[mdzi].op == '=');
                if(mdz[mdzi].op == '=') {
                    // nop
                } else {
                    assert(mdz[mdzi].op == 'X');
                    assert(strlen(mdz[mdzi].str) == run_comb);
                    int cread = bam_seqi(seq, seq_off);
                    if(!include_n_mms && run_comb == 1 && seq_nt16_str[cread] == 'N') {
                        // skip
                    } else {
                        fout << rec->core.tid << ',' << ref_off << ",X,";
                        seq_substring(fout, seq, seq_off, (size_t)run_comb);
                        if(print_qual) {
                            fout << ',';
                            cstr_substring(fout, qual, seq_off, (size_t)run_comb);
                        }
                        fout << '\n';
                    }
                }
                seq_off += run_comb;
                ref_off += run_comb;
                if(run_comb < mdz[mdzi].run) {
                    assert(mdz[mdzi].op == '=');
                    mdz[mdzi].run -= run_comb;
                } else {
                    mdzi++;
                }
            }
        } else if(op == BAM_CINS) {
            fout << rec->core.tid << ',' << ref_off << ",I,";
            seq_substring(fout, seq, seq_off, (size_t)run) << '\n';
            seq_off += run;
        } else if(op == BAM_CSOFT_CLIP) {
            if(include_sc) {
                char direction = '+';
                if(seq_off == 0)
                    direction = '-';
                (*total_softclip_count)+=run;
                if(only_polya_sc) {
                    char c;
                    int count_polya = polya_check(seq, seq_off, (size_t)run, &c);
                    if(count_polya != -1 && run >= SOFTCLIP_POLYA_TOTAL_COUNT_MIN) {
                        fout << rec->core.tid << ',' << ref_off << ",S,";
                        fout << run << ',' << direction << ',' << c << ',' << count_polya << '\n';
                    }
                }
                else {
                    fout << rec->core.tid << ',' << ref_off << ",S,";
                    seq_substring(fout, seq, seq_off, (size_t)run) << '\n';
                }
            }
            seq_off += run;
        } else if (op == BAM_CDEL) {
            assert(mdz[mdzi].op == '^');
            assert(run == mdz[mdzi].run);
            assert(strlen(mdz[mdzi].str) == run);
            mdzi++;
            fout << rec->core.tid << ',' << ref_off << ",D," << run << '\n';
            ref_off += run;
        } else if (op == BAM_CREF_SKIP) {
            ref_off += run;
        } else if (op == BAM_CHARD_CLIP) {
        } else if (op == BAM_CPAD) {
        } else {
            std::stringstream ss;
            ss << "No such CIGAR operation as \"" << op << "\"";
            throw std::runtime_error(ss.str());
        }
    }
    assert(mdzi == mdz.size());
}

static void output_from_cigar(const bam1_t *rec, std::fstream& fout, uint64_t* total_softclip_count, const bool include_sc, const bool only_polya_sc) {
    uint8_t *seq = bam_get_seq(rec);
    uint32_t *cigar = bam_get_cigar(rec);
    uint32_t n_cigar = rec->core.n_cigar;
    if(n_cigar == 1) {
        return;
    }
    int32_t refpos = rec->core.pos;
    int32_t seqpos = 0;
    for(uint32_t k = 0; k < n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        int run = bam_cigar_oplen(cigar[k]);
        switch(op) {
            case BAM_CDEL: {
                fout << rec->core.tid << ',' << refpos << ",D," << run << '\n';
                refpos += run;
                break;
            }
            case BAM_CSOFT_CLIP: {
                if(include_sc) {
                    char direction = '+';
                    if(seqpos == 0)
                        direction = '-';
                    (*total_softclip_count) += run;
                    if(only_polya_sc) {
                        char c;
                        int count_polya = polya_check(seq, (size_t)seqpos, (size_t)run, &c);
                        if(count_polya != -1 && run >= SOFTCLIP_POLYA_TOTAL_COUNT_MIN) {
                            fout << rec->core.tid << ',' << refpos << ',' << BAM_CIGAR_STR[op] << ',';
                            fout << run << ',' << direction << ',' << c << ',' << count_polya << '\n';
                        }
                    }
                    else {
                        fout << rec->core.tid << ',' << refpos << ',' << BAM_CIGAR_STR[op] << ',';
                        seq_substring(fout, seq, (size_t)seqpos, (size_t)run) << '\n';
                    }
                }
                seqpos += run;
                break;
            }
            case BAM_CINS: {
                fout << rec->core.tid << ',' << refpos << ',' << BAM_CIGAR_STR[op] << ',';
                seq_substring(fout, seq, (size_t)seqpos, (size_t)run) << '\n';
                seqpos += run;
                break;
            }
            case BAM_CREF_SKIP: {
                refpos += run;
                break;
            }
            case BAM_CMATCH:
            case BAM_CDIFF:
            case BAM_CEQUAL: {
                seqpos += run;
                refpos += run;
                break;
            }
            case 'H':
            case 'P': { break; }
            default: {
                std::stringstream ss;
                ss << "No such CIGAR operation as \"" << op << "\"";
                throw std::runtime_error(ss.str());
            }
        }
    }
}

static void print_header(const bam_hdr_t * hdr) {
    for(int32_t i = 0; i < hdr->n_targets; i++) {
        std::cout << '@' << i << ','
                  << hdr->target_name[i] << ','
                  << hdr->target_len[i] << std::endl;
    }
}

static const long get_longest_target_size(const bam_hdr_t * hdr) {
    long max = 0;
    for(int32_t i = 0; i < hdr->n_targets; i++) {
        if(hdr->target_len[i] > max)
            max = hdr->target_len[i];
    }
    return max;
}

static void reset_array(uint32_t* arr, const long arr_sz) {
#if USE_SIMD_ZERO
    #if __AVX2__
        __m256i zero = _mm256_setzero_si256();
        static constexpr size_t nper = sizeof(__m256i) / sizeof(uint32_t);
        const size_t nsimd = arr_sz / nper;
        const size_t nsimd4 = (nsimd / 4) * 4;
        size_t i = 0;
        for(; i < nsimd4; i += 4) {
            _mm256_storeu_si256((__m256i *)(arr + nper * i), zero);
            _mm256_storeu_si256((__m256i *)(arr + nper * (i + 1)), zero);
            _mm256_storeu_si256((__m256i *)(arr + nper * (i + 2)), zero);
            _mm256_storeu_si256((__m256i *)(arr + nper * (i + 3)), zero);
        }
        for(;i < nsimd; ++i) {
            _mm256_storeu_si256((__m256i *)(arr + nper * i), zero);
        }
        for(i *= sizeof(__m256i) / sizeof(uint32_t); i < arr_sz; ++i) {
            arr[i] = 0;
        }
    #elif __SSE2__
        __m128i zero = _mm_setzero_si128();
        const size_t nsimd = arr_sz / 4;
        const size_t nsimd4 = (nsimd / 4) * 4;
        size_t i = 0;
        for(; i < nsimd4; i += 4) {
            _mm_storeu_si128((__m128i *)(arr + 4 * i), zero);
            _mm_storeu_si128((__m128i *)(arr + 4 * (i + 1)), zero);
            _mm_storeu_si128((__m128i *)(arr + 4 * (i + 2)), zero);
            _mm_storeu_si128((__m128i *)(arr + 4 * (i + 3)), zero);
        }
        for(;i < nsimd; ++i) {
            _mm_storeu_si128((__m128i *)(arr + 4 * i), zero);
        }
        for(i *= 4; i < arr_sz; ++i) {
            arr[i] = 0;
        }
    #endif
#else
    std::memset(arr, 0, sizeof(uint32_t) * arr_sz);
#endif
}

template <typename T2>
static uint64_t print_array(const char* prefix,
                        char* chrm,
                        int32_t tid,
                        const T2* arr,
                        const long arr_sz,
                        const bool skip_zeros,
                        bigWigFile_t* bwfp,
                        FILE* cov_fh,
                        const bool dont_output_coverage = false,
                        bool no_region=true,
                        BGZF* gcov_fh = nullptr,
                        hts_idx_t* cidx = nullptr,
                        int* chrms_in_cidx = nullptr,
                        FILE* wcov_fh=nullptr,
                        BGZF* gwcov_fh=nullptr,
                        int window_size=0,
                        Op op = csum) {

    bool first = true;
    bool first_print = true;
    uint32_t running_value = 0;
    uint32_t last_pos = 0;
    uint64_t auc = 0;
    //from https://stackoverflow.com/questions/27401388/efficient-gzip-writing-with-gzprintf
    int chrnamelen = strlen(chrm);
    int total_line_len = chrnamelen + COORD_STR_LEN;
    int num_lines_per_buf = round(OUT_BUFF_SZ / total_line_len) - 3;
    int buf_written = 0;
    char* buf = nullptr;
    char* bufptr = nullptr;
    int (*printPtr) (void* fh, char* buf, uint32_t buf_len) = &my_write;
    void* cfh = nullptr;
    if(!bwfp) {
      buf = new char[OUT_BUFF_SZ];
      bufptr = buf;
      cfh = cov_fh;
      //writing gzip
      if(!cov_fh) {
        printPtr = &my_gzwrite;
        cfh = gcov_fh;
      }
    }

    //might only want to print windowed coverage
    bool print_windowed_coverage = window_size > 0 && (gwcov_fh || wcov_fh);
    void* wcfh = nullptr;
    if(print_windowed_coverage) {
      wcfh = wcov_fh; 
      //this assumes we're never going to have coverage and windowed coverage be different in terms of --gzip
      if(!wcov_fh) {
        printPtr = &my_gzwrite;
        wcfh = gwcov_fh; 
      }
    }


    uint32_t buf_len = 0;
    int bytes_written = 0;
    char* startp = new char[32];
    char* endp = new char[32];
    char* valuep = new char[32];
    float running_value_ = 0.0;
    uint32_t wcounter = 0;
    int64_t wsum = 0;
    char* wbuf = new char[1024];
    int window_bytes_written = -1;
    uint32_t window_start = 0;
    //make sure we track this chromosome in whatever index we're building
    //if we may it this far, means the chromosome had some alignments
    if(chrms_in_cidx && chrms_in_cidx[tid+1] == 0)
        chrms_in_cidx[tid+1] = ++chrms_in_cidx[0];

    for(uint32_t i = 0; i < arr_sz; i++) {
        if(first || (!no_region && running_value != arr[i]) || (no_region && arr[i] != 0)) {
            if(!first) {
                if(running_value > 0 || !skip_zeros) {
                    //based on wiggletools' AUC calculation
                    auc += (i - last_pos) * ((long) running_value);
                    if(not dont_output_coverage) {
                        if(bwfp && first_print) {
                            running_value_ = static_cast<float>(running_value);
                            bwAddIntervals(bwfp, &chrm, &last_pos, &i, &running_value_, 1);
                        }
                        else if(bwfp) {
                            running_value_ = static_cast<float>(running_value);
                            bwAppendIntervals(bwfp, &last_pos, &i, &running_value_, 1);
                        }
                        else {
                            memcpy(bufptr, chrm, chrnamelen);
                            char *oldbufptr = bufptr;
                            bufptr += chrnamelen;

                            *bufptr++='\t';
                            //idea from https://github.com/brentp/mosdepth/releases/tag/v0.2.9
                            uint32_t digits = u32toa_countlut(last_pos, bufptr, '\t');
                            bufptr+=digits+1;

                            digits = u32toa_countlut(i, bufptr, '\t');
                            bufptr+=digits+1;

                            digits = u32toa_countlut(running_value, bufptr, '\n');
                            bufptr+=digits+1;
                            buf_len += (bufptr - oldbufptr); // Track bytes written using the distance bufptr has traveled
                            bufptr[0]='\0';
                            (*printPtr)(cfh, buf, buf_len);
                            if(cidx) {
                                if(hts_idx_push(cidx, chrms_in_cidx[tid+1]-1, last_pos, i, bgzf_tell((BGZF*) cfh), 1) < 0) {
                                    fprintf(stderr,"error writing line in index at coordinates: %s:%u-%u, tid: %d idx tid: %d exiting\n",chrm,last_pos,i, tid, chrms_in_cidx[tid+1]-1);
                                    exit(-1);
                                }
                            }
                            buf_written++;
                            bufptr = buf;
                            buf_written = 0;
                            buf_len = 0;
                        } 
                        first_print = false;
                    }
                }
            }
            first = false;
            if(no_region)
                running_value += arr[i];
            else
                running_value = arr[i];
            last_pos = i;
        }
        if(print_windowed_coverage) {
            if(wcounter == window_size) {
                if(op == csum)
                    window_bytes_written = sprintf(wbuf, "%s\t%u\t%u\t%ld\n", chrm, window_start, i, wsum); 
                else if(op == cmean) {
                    double wmean = (double)wsum / (double)window_size;
                    window_bytes_written = sprintf(wbuf, "%s\t%u\t%u\t%.2f\n", chrm, window_start, i, wmean); 
                }

                (*printPtr)(wcfh, wbuf, window_bytes_written);
                wsum = 0;
                wcounter = 0;
                window_start = i;
            }
            wsum += running_value;
            wcounter++;
        }
    }
    char last_line[1024];
    if(!first) {
        if(running_value > 0 || !skip_zeros) {
            auc += (arr_sz - last_pos) * ((long) running_value);
            if(not dont_output_coverage) {
                if(bwfp) {
                    running_value_ = static_cast<float>(running_value);
                    if(first_print) {
                        bwAddIntervals(bwfp, &chrm, &last_pos, (uint32_t*) &arr_sz, &running_value_, 1);
                    } else {
                        bwAppendIntervals(bwfp, &last_pos, (uint32_t*) &arr_sz, &running_value_, 1);
                    }
                } else {
                    if(buf_written > 0) 
                        (*printPtr)(cfh, buf, buf_len);
                    // This printing step could also be u32toa_countlut-ified
                    buf_len = sprintf(last_line, "%s\t%u\t%lu\t%u\n", chrm, last_pos, arr_sz, running_value);
                    (*printPtr)(cfh, last_line, buf_len);
                    if(cidx)
                        if(hts_idx_push(cidx, chrms_in_cidx[tid+1]-1, last_pos, arr_sz, bgzf_tell((BGZF*) cfh), 1) < 0)
                            fprintf(stderr,"error writing last line of chromosome in index at coordinates: %s:%u-%ld, exiting\n",chrm,last_pos,arr_sz);
                }
            }
        }
        if(print_windowed_coverage) {
            if(op == csum)
                window_bytes_written = sprintf(wbuf, "%s\t%u\t%lu\t%ld\n", chrm, window_start, arr_sz, wsum); 
            else if(op == cmean) {
                window_size = arr_sz - window_start;
                double wmean = (double)wsum / (double)window_size;
                window_bytes_written = sprintf(wbuf, "%s\t%u\t%lu\t%.2f\n", chrm, window_start, arr_sz, wmean); 
            }
            (*printPtr)(wcfh, wbuf, window_bytes_written);
        }
    }
    return auc;
}

//generic function to loop through cigar
//and for each operation/lenth, call a list of functions to process
typedef std::vector<void*> args_list;
typedef void (*callback)(const int, const int, args_list*);
typedef std::vector<callback> callback_list;
static void process_cigar(int n_cigar, const uint32_t *cigar, char** cigar_str, callback_list* callbacks, args_list* outlist) {
    int cx = 0;
    for (int k = 0; k < n_cigar; ++k) {
        const int cigar_op = bam_cigar_op(cigar[k]);
        const int len = bam_cigar_oplen(cigar[k]);
        char op_char[2];
        op_char[0] = (char) bam_cigar_opchr(cigar[k]);
        op_char[1] = '\0';
        cx += sprintf((*cigar_str)+cx, "%d%s", len, op_char);
        int i = 0;
        //now call each callback function
        for(auto const& func : *callbacks) {
            (*func)(cigar_op, len, (args_list*) (*outlist)[i]);
            i++;
        }
    }
}

//mostly cribbed from htslib/sam.c
//calculates the mapped length of an alignment
static void maplength(const int op, const int len, args_list* out) {
    int type = bam_cigar_type(op);
    if ((type & 1) && (type & 2)) *((uint64_t*) (*out)[0]) += len;
}

static void end_genomic_coord(const int op, const int len, args_list* out) {
    int type = bam_cigar_type(op);
    if (type & 2) *((uint64_t*) (*out)[0]) += len;
}

static const int32_t align_length(const bam1_t *rec) {
    //bam_endpos processes the whole cigar string
    return bam_endpos(rec) - rec->core.pos;
}

typedef hashmap<std::string, char*> str2cstr;
typedef hashmap<std::string, int> str2int;
typedef std::vector<uint32_t> coords;
static void extract_junction(const int op, const int len, args_list* out) {
    uint32_t* base = (uint32_t*) (*out)[0];
    //not an intron
    if(op != BAM_CREF_SKIP) {
        //but track the length if consuming the ref
        if(bam_cigar_type(op) & 2)
            (*base) += len;
        return;
    }
    coords* jxs = (coords*) (*out)[1];
    jxs->push_back(*base);
    (*base) += len;
    jxs->push_back(*base);
}



static inline void decrement_coverages(uint32_t *coverages, uint32_t *unique_coverages, int start, int ninc, bool no_region=true) {
    coverages += start;
    unique_coverages += start;
    
    if(no_region) {
        int32_t* coverages_ = (int32_t*) coverages;
        int32_t* unique_coverages_ = (int32_t*) unique_coverages;

        coverages_[0]--;
        coverages_[ninc]++;

        unique_coverages_[0]--;
        unique_coverages_[ninc]++;
        return;
    }

#if __AVX512F__
    const size_t nper = sizeof(__m512) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __AVX2__
    const size_t nper = sizeof(__m256) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __SSE2__
    const size_t nper = sizeof(__m128) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#endif

    int i = 0;
#if __AVX512F__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm512_set1_epi32(-1);
        _mm512_storeu_si512((__m512i *)(coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(coverages + i * nper))));
        _mm512_storeu_si512((__m512i *)(unique_coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#elif __AVX2__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm256_set1_epi32(-1);
        _mm256_storeu_si256((__m256i *)(coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(coverages + i * nper))));
        _mm256_storeu_si256((__m256i *)(unique_coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#elif __SSE2__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm_set1_epi32(-1);
        _mm_storeu_si128((__m128i *)(coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(coverages + i * nper))));
        _mm_storeu_si128((__m128i *)(unique_coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#endif
    for(; i < ninc; ++i) {
        --coverages[i]; --unique_coverages[i];
    }
}

static inline void decrement_coverages(uint32_t *coverages, int ninc, bool no_region=true) {
    if(no_region) {
        int32_t* coverages_ = (int32_t*) coverages;
        coverages_[0]--;
        coverages_[ninc]++;
        return;
    }
#if __AVX512F__
    const size_t nper = sizeof(__m512) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __AVX2__
    const size_t nper = sizeof(__m256) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __SSE2__
    const size_t nper = sizeof(__m128) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#endif

    int i = 0;
#if __AVX512F__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm512_set1_epi32(-1);
        _mm512_storeu_si512((__m512i *)(coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(coverages + i * nper))));
    }
    i *= nper;
#elif __AVX2__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm256_set1_epi32(-1);
        _mm256_storeu_si256((__m256i *)(coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(coverages + i * nper))));
    }
    i *= nper;
#elif __SSE2__
    #pragma GCC unroll 4
    for(;i < nsimd; ++i) {
        auto s1 = _mm_set1_epi32(-1);
        _mm_storeu_si128((__m128i *)(coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(coverages + i * nper))));
    }
    i *= nper;
#endif
    for(; i < ninc; --coverages[i++]);
}

static inline void increment_coverages(uint32_t *coverages, int ninc, bool no_region=true) {
    if(no_region) {
        int32_t* coverages_ = (int32_t*) coverages;
        coverages_[0]++;
        coverages_[ninc]--;
        return;
    }
#if __AVX512F__
    const size_t nper = sizeof(__m512) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __AVX2__
    const size_t nper = sizeof(__m256) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __SSE2__
    const size_t nper = sizeof(__m128) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#endif

    int i = 0;
#if __AVX512F__
    for(;i < nsimd; ++i) {
        auto s1 = _mm512_set1_epi32(1);
        _mm512_storeu_si512((__m512i *)(coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(coverages + i * nper))));
    }
    i *= nper;
#elif __AVX2__
    for(;i < nsimd; ++i) {
        auto s1 = _mm256_set1_epi32(1);
        _mm256_storeu_si256((__m256i *)(coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(coverages + i * nper))));
    }
    i *= nper;
#elif __SSE2__
    for(;i < nsimd; ++i) {
        auto s1 = _mm_set1_epi32(1);
        _mm_storeu_si128((__m128i *)(coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(coverages + i * nper))));
    }
    i *= nper;
#endif
    for(; i < ninc; ++coverages[i++]);
}

static inline void increment_coverages(uint32_t *coverages, uint32_t *unique_coverages, int start, int ninc, bool no_region=true) {
    coverages += start;
    unique_coverages += start;
    if(no_region) {
        int32_t* coverages_ = (int32_t*) coverages;
        int32_t* unique_coverages_ = (int32_t*) unique_coverages;
        coverages_[0]++;
        coverages_[ninc]--;
        
        unique_coverages_[0]++;
        unique_coverages_[ninc]--;
        return;
    }
#if __AVX512F__
    const size_t nper = sizeof(__m512) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __AVX2__
    const size_t nper = sizeof(__m256) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#elif __SSE2__
    const size_t nper = sizeof(__m128) / sizeof(int32_t);
    size_t nsimd = ninc / nper;
#endif

    int i = 0;
#if __AVX512F__
    for(;i < nsimd; ++i) {
        auto s1 = _mm512_set1_epi32(1);
        _mm512_storeu_si512((__m512i *)(coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(coverages + i * nper))));
        _mm512_storeu_si512((__m512i *)(unique_coverages + i * nper), _mm512_add_epi32(s1, _mm512_loadu_si512((__m512i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#elif __AVX2__
    for(;i < nsimd; ++i) {
        auto s1 = _mm256_set1_epi32(1);
        _mm256_storeu_si256((__m256i *)(coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(coverages + i * nper))));
        _mm256_storeu_si256((__m256i *)(unique_coverages + i * nper), _mm256_add_epi32(s1, _mm256_loadu_si256((__m256i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#elif __SSE2__
    for(;i < nsimd; ++i) {
        auto s1 = _mm_set1_epi32(1);
        _mm_storeu_si128((__m128i *)(coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(coverages + i * nper))));
        _mm_storeu_si128((__m128i *)(unique_coverages + i * nper), _mm_add_epi32(s1, _mm_loadu_si128((__m128i *)(unique_coverages + i * nper))));
    }
    i *= nper;
#endif
    for(; i < ninc; ++i) {
        ++coverages[i]; ++unique_coverages[i];
    }
}

static uint64_t num_overlapping_pairs = 0;
//static uint32_t num_opairs[10024];

struct MateInfo {
    bool passing_qual;
    std::string qname;
    //char* qname;
    int32_t mrefpos;
    uint32_t n_cigar;
    uint32_t* cigar;
    bool erased;
};

//typedef hashmap<std::string, uint32_t*> read2len;
//typedef hashmap<uint32_t, uint32_t*> read2len;
typedef hashmap<uint32_t, std::vector<MateInfo*>*> read2len;
static const int32_t calculate_coverage(const bam1_t *rec, uint32_t* coverages,
                                        uint32_t* unique_coverages, const bool double_count,
                                        const int min_qual, read2len* overlapping_mates,
                                        int32_t* total_intron_length, bool no_region=true) {
    int32_t refpos = rec->core.pos;
    int32_t mrefpos = rec->core.mpos;
    int32_t refpos_to_hash = mrefpos;
    //lifted from htslib's bam_cigar2rlen(...) & bam_endpos(...)
    int32_t algn_end_pos = refpos;
    const uint32_t* cigar = bam_get_cigar(rec);
    int k, z;
    //check for overlapping mate and corect double counting if exists
    char* qname = bam_get_qname(rec);
    bool unique = min_qual > 0;
    bool passing_qual = rec->core.qual >= min_qual;
    //fix paired mate overlap double counting
    //fix overlapping mate pair, only if 1) 2nd mate and
    //2) we're either not unique, or we're higher than the required quality
    int32_t mendpos = 0;
    int n_mspans = 0;
    std::unique_ptr<int32_t[]> mspans;
    int mspans_idx = 0;
    const std::string tn(qname);
    int32_t end_pos = bam_endpos(rec);
    uint32_t mate_passes_quality = 0;
    //-----First Mate Check
    //if we're the first mate and
    //we're avoiding double counting and we're a proper pair
    //and we overlap with our mate, then store our cigar + length
    //for the later mate to adjust its coverage appropriately
    if(coverages && !double_count && (rec->core.flag & BAM_FPROPER_PAIR) == 2) {
        bool possible_overlap = rec->core.tid == rec->core.mtid && end_pos > mrefpos;
        bool first_mate_w_overlap = possible_overlap && refpos <= mrefpos;
        bool second_mate = possible_overlap && refpos >= mrefpos;
        if(second_mate)
            refpos_to_hash = refpos;
        //1) we're on the same chrm as our mate AND
        //2) we're either the first mate overlapping with the 2nd, or we're the 2nd mate
        //so we could have mate overlap
        if(first_mate_w_overlap || second_mate) {
            std::vector<MateInfo*>* mate_vec = nullptr;
            MateInfo* mate_info = nullptr;
            
            auto mit = overlapping_mates->find(refpos_to_hash);
            bool potential_mate_found = mit != overlapping_mates->end();

            //if we found a potential mate in the hash based on pos
            int mvi = 0;
            if(potential_mate_found) {
                mate_vec = mit->second;
                for(auto mate : *mate_vec) {
                    //fprintf(stderr,"name check for refpos %u mrefpos %u: %s vs. %s\n",refpos, mrefpos, tn.c_str(), mate->qname);
                    if(!mate->erased && mate->qname == tn) {
                        mate_info = mate;
                        break;
                    }
                    mvi++;
                }
            }

            //first mate in the pair
            if(first_mate_w_overlap && !mate_info) {
                const uint32_t* mcigar = bam_get_cigar(rec);
                uint32_t n_cigar = rec->core.n_cigar;
                mate_info = new MateInfo;
                mate_info->passing_qual = unique && passing_qual;
                mate_info->qname = tn;
                mate_info->mrefpos = refpos;
                mate_info->n_cigar = n_cigar;
                mate_info->cigar = new uint32_t[n_cigar];
                std::memcpy(mate_info->cigar, mcigar, 4*n_cigar);
                mate_info->erased = false;
                //if we didn't find a previous vector, create one
                if(!potential_mate_found) {
                    mate_vec = new std::vector<MateInfo*>;
                    overlapping_mates->emplace(mrefpos, mate_vec);
                }
                mate_vec->push_back(mate_info);
                num_overlapping_pairs++;
            }
            //-------Second Mate Check
            else if(second_mate && mate_info) {
                uint32_t mn_cigar = mate_info->n_cigar;
                mate_passes_quality = mate_info->passing_qual;
                uint32_t* mcigar = mate_info->cigar;
                int32_t real_mate_pos = mate_info->mrefpos;
                int32_t malgn_end_pos = real_mate_pos;
                //bash cigar to get spans of overlap
                mspans.reset(new int32_t[mn_cigar * 2]);
                for (k = 0; k < mn_cigar; ++k) {
                    const int cigar_op = bam_cigar_op(mcigar[k]);
                    if(bam_cigar_type(cigar_op)&2) {
                        const int32_t len = bam_cigar_oplen(mcigar[k]);
                        if(bam_cigar_type(cigar_op)&1) {
                            mspans[mspans_idx * 2] = malgn_end_pos;
                            mspans[mspans_idx * 2 + 1] = malgn_end_pos + len;
                            mspans_idx++;
                        }
                        malgn_end_pos += len;
                    }
                }
                delete[] mcigar;
                mate_info->erased = true;
                //overlapping_mates->erase(mit);
                delete mate_info;
                mate_vec->erase(mate_vec->begin()+mvi);
                if(mate_vec->size() == 0) {
                    //mate_vec->shrink_to_fit();
                    //std::vector<MateInfo*>().swap(*mate_vec);
                    //fprintf(stderr, "erasing vector\n");
                    //delete mate_vec;
                    delete mate_vec;
                    overlapping_mates->erase(mit);
                }
                n_mspans = mspans_idx;
                mendpos = malgn_end_pos;
            }
        }
    }
    mspans_idx = 0;
    if(unique && passing_qual) {
        int32_t lastref = 0;
        for (k = 0; k < rec->core.n_cigar; ++k) {
            const int cigar_op = bam_cigar_op(cigar[k]);
            //do we consume ref?
            if(bam_cigar_type(cigar_op)&2) {
                const int32_t len = bam_cigar_oplen(cigar[k]);
                if(cigar_op == BAM_CREF_SKIP)
                    (*total_intron_length) = (*total_intron_length) + len;
                //are we calc coverages && do we consume query?
                if(coverages && bam_cigar_type(cigar_op)&1) {
                    increment_coverages(coverages, unique_coverages, algn_end_pos, len, no_region);
                    //now fixup overlapping segment but only if mate passed quality
                    if(n_mspans > 0 && algn_end_pos < mendpos) {
                        //loop until we find the next overlapping span
                        //if are current segment is too early we just keep the span index where it is
                        while(mspans_idx < n_mspans && algn_end_pos >= mspans[mspans_idx * 2 + 1])
                            mspans_idx++;
                        int32_t cur_end = algn_end_pos + len;
                        int32_t left_end = algn_end_pos;
                        if(left_end < mspans[mspans_idx * 2])
                            left_end = mspans[mspans_idx * 2];
                        //check 1) we've still got mate spans 2) current segment overlaps the current mate span
                        while(mspans_idx < n_mspans && left_end < mspans[mspans_idx * 2 + 1]
                                                    && cur_end > mspans[mspans_idx * 2]) {
                            //set right end of segment to decrement
                            int32_t right_end = cur_end;
                            int32_t next_left_end = left_end;
                            if(right_end >= mspans[mspans_idx * 2 + 1]) {
                                right_end = mspans[mspans_idx * 2 + 1];
                                //if our segment is greater than the previous mate's
                                //also increment the mate spans index
                                mspans_idx++;
                                if(mspans_idx < n_mspans)
                                    next_left_end = mspans[mspans_idx * 2];
                            }
                            else {
                                next_left_end = mspans[mspans_idx * 2 + 1];
                            }
                            decrement_coverages(coverages + left_end, right_end - left_end, no_region);
                            if(mate_passes_quality)
                                decrement_coverages(unique_coverages + left_end, right_end - left_end, no_region);
                            left_end = next_left_end;
                        }
                    }
                }
                algn_end_pos += len;
            }
        }
    } else {
        for (k = 0; k < rec->core.n_cigar; ++k) {
            const int cigar_op = bam_cigar_op(cigar[k]);
            //do we consume ref?
            if(bam_cigar_type(cigar_op)&2) {
                const int32_t len = bam_cigar_oplen(cigar[k]);
                if(cigar_op == BAM_CREF_SKIP)
                    (*total_intron_length) = (*total_intron_length) + len;
                //are we calc coverages && do we consume query?
                if(coverages && bam_cigar_type(cigar_op)&1) {
                    increment_coverages(&coverages[algn_end_pos], len, no_region);
                    //now fixup overlapping segment
                    if(n_mspans > 0 && algn_end_pos < mendpos) {
                        //loop until we find the next overlapping span
                        //if are current segment is too early we just keep the span index where it is
                        while(mspans_idx < n_mspans && algn_end_pos >= mspans[mspans_idx * 2 + 1])
                            mspans_idx++;
                        int32_t cur_end = algn_end_pos + len;
                        int32_t left_end = algn_end_pos;
                        if(left_end < mspans[mspans_idx * 2])
                            left_end = mspans[mspans_idx * 2];
                        //check 1) we've still got mate spans 2) current segment overlaps the current mate span
                        while(mspans_idx < n_mspans && left_end < mspans[mspans_idx * 2 + 1]
                                                    && cur_end > mspans[mspans_idx * 2]) {
                            //set right end of segment to decrement
                            int32_t right_end = cur_end;
                            int32_t next_left_end = left_end;
                            if(right_end >= mspans[mspans_idx * 2 + 1]) {
                                right_end = mspans[mspans_idx * 2 + 1];
                                //if our segment is greater than the previous mate's
                                //also increment the mate spans index
                                //delete[] mspans[mspans_idx];
                                mspans_idx++;
                                if(mspans_idx < n_mspans)
                                    next_left_end = mspans[mspans_idx * 2];
                            }
                            else {
                                next_left_end = mspans[mspans_idx * 2 + 1];
                            }
                            decrement_coverages(&coverages[left_end], right_end - left_end, no_region);
                            left_end = next_left_end;
                        }
                    }
                }
                algn_end_pos += len;
            }
        }
    }
    return algn_end_pos;
}

template <typename T>
using annotation_map_t = hashmap<std::string, std::vector<T*>>;
typedef std::vector<char*> strlist;
//about 3x faster than the sstring/string::getline version
template <typename T>
static const int process_region_line(char* line, const char* delim, annotation_map_t<T>* amap, strlist* chrm_order, bool keep_order) {
    char* tok = strtok(line, delim);
    int i = 0;
    char* chrm = nullptr;
    long start = -1;
    long end = -1;
    int ret = 0;
    int last_col = END_COL;
    while(tok != nullptr) {
        if(i > last_col)
            break;
        if(i == CHRM_COL)
            chrm = strdup(tok);
        else if(i == START_COL)
            start = atol(tok);
        else if(i == END_COL)
            end = atol(tok);
        i++;
        tok = strtok(nullptr, delim);
    }
    //if we need to keep the order, then we'll store values here
    const int alen = keep_order?4:2;
    T* coords = new T[alen];
    coords[0] = start;
    coords[1] = end;
    std::fill(coords + 2, coords + alen, 0);
    auto it = amap->find(chrm);
    if(it == amap->end()) {
        chrm_order->push_back(chrm);
        it = amap->emplace(chrm, std::vector<T*>()).first;
    }
    it->second.push_back(coords);
    return ret;
}

template <typename T>
static const int read_annotation(FILE* fin, annotation_map_t<T>* amap, strlist* chrm_order, bool keep_order, uint64_t* num_annotations) {
    char *line = (char *)std::malloc(LINE_BUFFER_LENGTH);
    size_t length = LINE_BUFFER_LENGTH;
    assert(fin);
    ssize_t bytes_read = getline(&line, &length, fin);
    int err = 0;
    while(bytes_read != -1) {
        err = process_region_line(line, "\t", amap, chrm_order, keep_order);
        if(err) {
            std::cerr << "Error: " << err << " in process_region_line.\n";
            break;
        }
        assert(err==0);
        (*num_annotations)++;
        bytes_read = getline(&line, &length, fin);
    }
    std::free(line);
    std::cerr << "building whole annotation region map done\n";
    return err;
}

typedef hashmap<std::string, int> str2op;

template <typename T>
static void sum_annotations(const uint32_t* coverages, const std::vector<T*>& annotations, const long chr_size, const char* chrm, FILE* ofp, uint64_t* annotated_auc, Op op, bool just_auc = false, int keep_order_idx = -1) {
    unsigned long z, j;
    int (*printPtr) (char* buf, const char*, long, long, T, double*, long) = &print_shared;
    int (*outputFunc)(void* fh, char* buf, uint32_t buf_len) = &my_write;
    if(SUMS_ONLY)
        printPtr = &print_shared_sums_only;
    char* buf = new char[1024];
    for(z = 0; z < annotations.size(); z++) {
        T sum = 0;
        T start = annotations[z][0];
        T end = annotations[z][1];
        T local_sum = 0;
        for(j = start; j < end; j++) {
            assert(j < chr_size);
            local_sum += coverages[j];
        }
        sum += local_sum;
        (*annotated_auc) = (*annotated_auc) + sum;
        if(!just_auc) {
            if(op == cmean) 
                sum = (double)local_sum / ((double)(end-start));
            if(keep_order_idx == -1) {
                int buf_len = (*printPtr)(buf, chrm, (long) start, (long) end, sum, nullptr, 0);
                (*outputFunc)(ofp, buf, buf_len);
            }
            else
                annotations[z][keep_order_idx] = sum;
        }
    }
}


static bigWigFile_t* create_bigwig_file(const bam_hdr_t *hdr, const char* out_fn, const char *suffix) {
    if(bwInit(BW_READ_BUFFER) != 0) {
        fprintf(stderr, "Failed when calling bwInit with %d init val\n", BIGWIG_INIT_VAL);
        exit(-1);
    }
    char fn[1024] = "";
    sprintf(fn, "%s.%s", out_fn, suffix);
    bigWigFile_t* bwfp = bwOpen(fn, nullptr, "w");
    if(!bwfp) {
        fprintf(stderr, "Failed when attempting to open BigWig file %s for writing\n", fn);
        exit(-1);
    }
    //create with up to 10 zoom levels (though probably less in practice)
    bwCreateHdr(bwfp, 10);
    bwfp->cl = bwCreateChromList(hdr->target_name, hdr->target_len, hdr->n_targets);
    bwWriteHdr(bwfp);
    return bwfp;
}

int KALLISTO_MAX_FRAG_LENGTH = 1000;
typedef hashmap<int32_t, uint32_t> fraglen2count;
static void print_frag_distribution(const fraglen2count* frag_dist, FILE* outfn)
{
    double mean = 0.0;
    uint64_t count = 0;
    //track a Kallisto-comparable version separately
    double kmean = 0.0;
    uint64_t kcount = 0;
    uint64_t mode = 0;
    uint64_t mode_count = 0;
    for(auto kv: *frag_dist) {
        fprintf(outfn, "%d\t%u\n", kv.first, kv.second);
        count += kv.second;
        mean += (kv.first*kv.second);
        if(kv.first < KALLISTO_MAX_FRAG_LENGTH) {
            kcount += kv.second;
            kmean += (kv.first*kv.second);
        }
        if(kv.second > mode_count) {
            mode_count = kv.second;
            mode = kv.first;
        }
    }
    mean /= count;
    kmean /= kcount;
    fprintf(outfn, "STAT\tCOUNT\t%" PRIu64 "\n", count);
    fprintf(outfn, "STAT\tMEAN_LENGTH\t%.3f\n", mean);
    fprintf(outfn, "STAT\tMODE_LENGTH\t%" PRIu64 "\n", mode);
    fprintf(outfn, "STAT\tMODE_LENGTH_COUNT\t%" PRIu64 "\n", mode_count);
    fprintf(outfn, "STAT\tKALLISTO_COUNT\t%" PRIu64 "\n", kcount);
    fprintf(outfn, "STAT\tKALLISTO_MEAN_LENGTH\t%.3f\n", kmean);
}

void output_read_sequence_and_qualities(char* qname, int midx, uint8_t* seq, uint8_t* qual, size_t l_qseq, bool reversed, std::ostream* outfh, bool one_file) {
    (*outfh) << "@" << qname;
    if(!one_file)
        (*outfh) << "/" << midx;
    (*outfh) << "\n";
    seq_substring(*outfh, seq, 0, l_qseq, reversed);
    (*outfh) << "\n+\n";
    qstr_substring(*outfh, qual, 0, l_qseq, reversed);
    (*outfh) << "\n";
}


static int process_bigwig_for_total_auc(const char* fn, double* all_auc, FILE* errfp = stderr) {
    //in part lifted from https://github.com/dpryan79/libBigWig/blob/master/test/testIterator.c
    //this is the buffer
    if(bwInit(BW_READ_BUFFER) != 0) {
        fprintf(errfp, "Error in bwInit, exiting\n");
        return -1;
    }
    bigWigFile_t *fp = bwOpen((char *)fn, NULL, "r");
    if(!fp) {
        fprintf(errfp, "Error in opening %s as BigWig file, exiting\n", fn);
        return -1;
    }
    fprintf(stdout,"opened %s, BW read buffer is %u\n",fn, BW_READ_BUFFER);
    fflush(stdout);
    uint32_t i, tid, blocksPerIteration;
    //better to ask for a few blocks for better memory and time stats
    blocksPerIteration = 10;
    bwOverlapIterator_t *iter = nullptr;
    uint64_t total_num_intervals = 0;
    //loop through all the chromosomes in the BW
    for(tid = 0; tid < fp->cl->nKeys; tid++)
    {
        if(fp->cl->len[tid] < 1)
            continue;
        iter = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[tid], 0, fp->cl->len[tid], blocksPerIteration);

        if(!iter->data)
        {
            fprintf(errfp, "WARNING: no intervals for chromosome %s in %s as BigWig file, skipping\n", fp->cl->chrom[tid], fn);
            goto next;
            continue;
        }
        while(iter->data)
        {
            uint32_t num_intervals = iter->intervals->l;
            total_num_intervals+=num_intervals;
            uint32_t istart = 0;
            uint32_t iend = 0;
            for(int j = 0; j < num_intervals; j++)
            {
                istart = iter->intervals->start[j];
                iend = iter->intervals->end[j];
                double value = (iend-istart) * iter->intervals->value[j];
                (*all_auc) += value;
            }
            iter = bwIteratorNext(iter);
        }
        next: // To ensure that we are destroying for cases where no intervals are available (1115)
              // Could replace with RAII, but this is simpler and fits the style better
        bwIteratorDestroy(iter);
    }

    bwClose(fp);
    bwCleanup();
    return 0;
}


using chr2bool = hashset<std::string>;
template <typename T>
static int process_bigwig(const char* fn, double* annotated_auc, annotation_map_t<T>* amap, chr2bool* annotation_chrs_seen, FILE* afp, int keep_order_idx = -1, Op op = csum, FILE* errfp = stderr, str2dblist* store_local=nullptr) {
    //in part lifted from https://github.com/dpryan79/libBigWig/blob/master/test/testIterator.c
    if(bwInit(BW_READ_BUFFER) != 0) {
        fprintf(errfp, "Error in bwInit, exiting\n");
        return -1;
    }
    bigWigFile_t *fp = bwOpen((char *)fn, NULL, "r");
    if(!fp) {
        fprintf(errfp, "Error in opening %s as BigWig file, exiting\n", fn);
        return -1;
    }
    int (*printPtr) (char* buf, const char*, long, long, T, double*, long) = &print_shared;
    int (*outputFunc)(void* fh, char* buf, uint32_t buf_len) = &my_write;
    if(SUMS_ONLY)
        printPtr = &print_shared_sums_only;
    char* buf = new char[1024];
    uint32_t tid, blocksPerIteration;
    //ask for huge # of blocks per chromosome to ensure we get all in one go
    //(this is for convenience, not performance)
    blocksPerIteration = 4000000;
    //blocksPerIteration = 1;
    bwOverlapIterator_t *iter = nullptr;
    //for certain modes we only want to process BW intervals once
    std::vector<bool> intervals_seen(STARTING_NUM_INTERVALS,false);
    //loop through all the chromosomes in the BW
    for(tid = 0; tid < fp->cl->nKeys; tid++)
    {
        //only process the chromosome if it's in the annotation
        if(amap->find(fp->cl->chrom[tid]) != amap->end()) {
            iter = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[tid], 0, fp->cl->len[tid], blocksPerIteration);
            if(!iter->data)
            {
                fprintf(errfp, "WARNING: no interval data for chromosome %s in %s as BigWig file, skipping\n", fp->cl->chrom[tid], fn);
                continue;
            }
            uint32_t num_intervals = iter->intervals->l;
            if(num_intervals == 0) {
                fprintf(errfp, "WARNING: 0 intervals for chromosome %s in %s as BigWig file, skipping\n", fp->cl->chrom[tid], fn);
                continue;
            }
            if(op == csum) {
                if(num_intervals > intervals_seen.size())
                    intervals_seen.resize(num_intervals, false);
                for(unsigned int i = 0; i < intervals_seen.size(); i++)
                    intervals_seen[i] = false;
            }
            uint32_t istart = iter->intervals->start[0];
            uint32_t iend = iter->intervals->end[num_intervals-1];
            std::vector<T*>& annotations = amap->operator[](fp->cl->chrom[tid]);
            long z, j, k;
            long last_j = 0;
            long asz = annotations.size();
            double* local_vals;
            //if running in multithreaded mode, want to store the values locally
            //but also don't want to reallocate for every new bigwig file, so
            //we allocate once per thread per chromosome
            if(store_local) {
                if(store_local->find(fp->cl->chrom[tid]) == store_local->end())
                    local_vals = new double[asz];
                else
                    local_vals = (*store_local)[fp->cl->chrom[tid]];
                std::fill(local_vals, local_vals + asz, 0.);
            }
            //loop through annotation intervals as outer loop
            for(z = 0; z < asz; z++) {
                const auto &az = annotations[z];
                double sum = 0;
                double min = MAX_INT;
                double max = 0;
                T start = az[0];
                T ostart = start;
                T end = az[1];
                //find the first BW interval starting *before* our annotation interval
                //this is if we have overlapping/out-of-order intervals in the annotation
                while(start < iter->intervals->start[last_j])
                    last_j--;
                for(j = last_j; j < num_intervals; j++)
                {
                    istart = iter->intervals->start[j];
                    iend = iter->intervals->end[j];
                    //is our start overlapping?
                    if(start >= istart && start < iend)
                    {
                        long last_k = end > iend ? iend : end;
                        //stat mode
                        //avoid having if's in the inner loops as much as possible
                        switch(op) {
                            case csum:
                            case cmean:
                                for(k = start; k < last_k; k++)
                                    sum += iter->intervals->value[j];
                                break;
                            case cmin:
                                for(k = start; k < last_k; k++)
                                    min = iter->intervals->value[j] < min ? iter->intervals->value[j]:min;
                                break;
                            case cmax:
                                for(k = start; k < last_k; k++)
                                    max = iter->intervals->value[j] > max ? iter->intervals->value[j]:max;
                                break;
                        }

                        //move start up
                        if(k < end)
                            start = k;
                        //break out if we've hit the end of this annotation interval
                        if(k >= end)
                            break;
                    }
                }
                last_j = j;
                if(op == csum)
                    (*annotated_auc) += sum;
                //0-based start
                double annot_length = end - ostart;
                T value = sum;
                switch(op) {
                    case cmean:
                        value = (double)sum / (double)annot_length;
                        break;
                    case cmin:
                        value = min;
                        break;
                    case cmax:
                        value = max;
                        break;
                    case csum:; // do nothing
                }
                //not trying to keep the order in the BED file, just print them as we find them
                if(keep_order_idx == -1) {
                    int buf_len = (*printPtr)(buf, fp->cl->chrom[tid], (long) ostart, (long) end, value, nullptr, 0);
                    (*outputFunc)(afp, buf, buf_len);
                }
                else if(store_local)
                    local_vals[z] = value;
                else
                    az[keep_order_idx] = value;
            }
            annotation_chrs_seen->insert(fp->cl->chrom[tid]);
            if(store_local)
                (*store_local)[fp->cl->chrom[tid]] = local_vals;
            bwIteratorDestroy(iter);
        }
    }

    bwClose(fp);
    bwCleanup();
    return 0;
}


template <typename T>
static void output_missing_annotations(const annotation_map_t<T>* annotations, const chr2bool* annotations_seen, FILE* ofp, Op op = csum) {
    //check if we're doing means output doubles, otherwise output longs
    T val = 0;
    int (*printPtr) (char* buf, const char*, long, long, T, double*, long) = &print_shared;
    int (*outputFunc)(void* fh, char* buf, uint32_t buf_len) = &my_write;
    if(SUMS_ONLY)
        printPtr = &print_shared_sums_only;
    char* buf = new char[1024];
    for(auto const& kv : *annotations) {
        if(annotations_seen->find(kv.first) == annotations_seen->end()) {
            const auto &ants = kv.second;
            for(unsigned long z = 0; z < kv.second.size(); z++) {
                const auto p = ants[z];
                int buf_len = (*printPtr)(buf, kv.first.c_str(), p[0], p[1], val, nullptr, z);
                (*outputFunc)(ofp, buf, buf_len);
            }
        }
    }
}

template <typename T>
void output_all_coverage_ordered_by_BED(const strlist* chrm_order, annotation_map_t<T>* annotations, FILE* afp, BGZF* afpz, FILE* uafp,BGZF* uafpz, Op op = csum, str2dblist* store_local = nullptr) {
    int (*outputFunc)(void* fh, char* buf, uint32_t buf_len) = &my_write;
    void* out_fh = afp;
    void* uout_fh = uafp;
    if(afpz) {
        outputFunc = &my_gzwrite;
        out_fh = afpz;
    }
    if(uafpz)
        uout_fh = uafpz;
    double* local_vals = nullptr;
    for(auto const c : *chrm_order) {
        if(!c)
            continue;
        std::vector<T*>& annotations_for_chr = (*annotations)[c];
        int (*printPtr) (char*, const char*, long, long, T, double*, long) = &print_shared;
        if(SUMS_ONLY)
            printPtr = &print_shared_sums_only;
        if(store_local) {
            local_vals = (*store_local)[c];
            printPtr = &print_local;
            if(SUMS_ONLY)
                printPtr = &print_local_sums_only;
        }
        //check if we're doing means output doubles, otherwise output longs
        char* buf = new char[OUT_BUFF_SZ];
        char* bufptr = buf;
        int buf_len = 0;
        int buf_written = 0;
        //unique
        char* ubuf = nullptr;
        if(uafp)
            ubuf = new char[OUT_BUFF_SZ];
        char* ubufptr = ubuf;
        int ubuf_len = 0;
        int ubuf_written = 0;
        int num_lines_per_buf = round(OUT_BUFF_SZ / COORD_STR_LEN) - 3;
        for(long z = 0; z < annotations_for_chr.size(); z++) {
            const auto &item = annotations_for_chr[z];
            const T start = item[0], end = item[1];
            T val = item[2];
            if(buf_written >= num_lines_per_buf) {
                bufptr[0]='\0';
                (*outputFunc)(out_fh, buf, buf_len);
                bufptr = buf;
                buf_written = 0;
                buf_len = 0;
            }
            int written = (*printPtr)(bufptr, c, (long) start, (long) end, val, local_vals, z);
            bufptr += written;
            buf_len += written;
            buf_written++;
            //do uniques if asked to
            if(uafp) {
                val = item[3];
                if(ubuf_written >= num_lines_per_buf) {
                    ubufptr[0]='\0';
                    (*outputFunc)(uout_fh, ubuf, ubuf_len);
                    ubufptr = ubuf;
                    ubuf_written = 0;
                    ubuf_len = 0;
                }
                written = (*printPtr)(ubufptr, c, (long) start, (long) end, val, local_vals, z);
                ubufptr += written;
                ubuf_len += written;
                ubuf_written++;
            }
        }
        char last_line[1024];
        if(buf_written > 0) {
            bufptr[0]='\0';
            (*outputFunc)(out_fh, buf, buf_len);
        }
        if(ubuf_written > 0) {
            ubufptr[0]='\0';
            (*outputFunc)(uout_fh, ubuf, ubuf_len);
        }
    }
}

//multiple sources for this kind of tokenization, one which was useful was:
//https://yunmingzhang.wordpress.com/2015/07/14/how-to-read-file-line-by-lien-and-split-a-string-in-c/
void split_string(std::string line, char delim, strvec* tokens) {
    tokens->clear();
    std::stringstream ss(line);
    std::string token;
    while(getline(ss, token, delim))
    {
        tokens->push_back(token);
    }
}

template <typename T>
void process_bigwig_worker(strvec& bwfns, annotation_map_t<T>* annotations, strlist* chrm_order, int keep_order_idx, Op op) {
    //want to just get the filename itself, no path
    str2dblist store_local;
    for(auto bwfn_ : bwfns) {
        strvec tokens;
        const char* bwfn = bwfn_.c_str();
        fprintf(stderr, "about to process %s\n", bwfn);
        std::string str(bwfn_);
        split_string(str, '/', &tokens);
        char afn[1024];
        FILE* afp = nullptr;
        sprintf(afn, "%s.err", tokens.back().c_str());
        FILE* errfp = fopen(afn, "w");
        sprintf(afn, "%s.all.tsv", tokens.back().c_str());
        afp = fopen(afn, "w");
        chr2bool annotation_chrs_seen;
        double annotated_auc = 0.0;

        int ret = process_bigwig(bwfn, &annotated_auc, annotations, &annotation_chrs_seen, afp, keep_order_idx, op = op, errfp = errfp, &store_local);
        if(ret != 0) {
            fprintf(errfp,"FAILED to process bigwig %s\n", bwfn);
            if(afp)
                fclose(afp);
            fclose(errfp);
            return;
        }
        //if we wanted to keep the chromosome order of the annotation output matching the input BED file
        if(keep_order_idx == 2)
            output_all_coverage_ordered_by_BED(chrm_order, annotations, afp, nullptr, nullptr, nullptr, op, &store_local);
        else
            output_missing_annotations(annotations, &annotation_chrs_seen, afp, op = op);
        if(afp)
            fclose(afp);
        //fprintf(aucfp, "AUC\t%" PRIu64 "\n", annotated_auc);
        fprintf(stdout, "AUC_ANNOTATED_BASES\t%.3f\t%s\n", annotated_auc, bwfn);
        //fprintf(errfp, "AUC\t%.3f\n", annotated_auc);
        fprintf(errfp,"SUCCESS processing bigwig %s\n", bwfn);
        fclose(errfp);
    }
    //hold off on final deletion, for performance
    /*for( auto mitr : store_local)
        delete mitr.second;*/
}

Op get_operation(const char* opstr) {
    if(strcmp(opstr, "mean") == 0)
        return cmean;
    if(strcmp(opstr, "min") == 0)
        return cmin;
    if(strcmp(opstr, "max") == 0)
        return cmax;
    return csum;
}


typedef hashmap<std::string, uint8_t*> str2str;
static const uint64_t frag_lens_mask = 0x00000000FFFFFFFF;
static const int FRAG_LEN_BITLEN = 32;
template <typename T>
int go_bw(const char* bw_arg, int argc, const char** argv, Op op, htsFile *bam_fh, int nthreads, bool keep_order, bool has_annotation, FILE* afp, BGZF* afpz, annotation_map_t<T>* annotations, chr2bool* annotation_chrs_seen, const char* prefix, bool sum_annotation, strlist* chrm_order, FILE* auc_file, uint64_t num_annotations) {
    //only calculate AUC across either the BAM or the BigWig, but could be restricting to an annotation as well
    int err = 0;
    bool LOAD_BALANCE = false;
    int slen = strlen(bw_arg);
    bool is_bw_list_file = strcmp(bw_arg+(slen-3), "txt") == 0;
    fprintf(stderr,"Processing %s\n",bw_arg);
    fflush(stderr);
    //just do all/total AUC if no options are passed in
    if(argc == 1
            || (argc == 2 && has_option(argv, argv+argc, "--auc"))
            || (argc == 3 && has_option(argv, argv+argc, "--bwbuffer"))
            || (argc == 4 && has_option(argv, argv+argc, "--bwbuffer") && has_option(argv, argv+argc, "--auc"))) {
        //should be the same as "all_auc" except support possibility of continuous values
        //in the BigWig (but not in the BAM, since we control how we count)
        double total_auc = 0.0;
        int ret = process_bigwig_for_total_auc(bw_arg, &total_auc);
        if(ret == 0)
            fprintf(stdout, "AUC_ALL_BASES\t%.3f\n", total_auc);
        return ret;
    }

    double annotated_total_auc = 0.0;
    //process bigwig for annotation/auc
    int keep_order_idx = keep_order?2:-1;
    //TODO: look into implemention multithreaded mode for single BigWig processing (maybe per chromosome?)
    if(is_bw_list_file) {
        strvec* files_per_thread[nthreads];
        uint64_t bytes_per_thread[nthreads];
        for(int i=0; i < nthreads; i++) {
            files_per_thread[i] = new strvec();
            bytes_per_thread[i] = 0;
        }
        FILE* bw_list_fp = fopen(bw_arg, "r");
        if(unlikely(bw_list_fp == nullptr)) assert(false);
        char *bwfn = (char *)std::malloc(LINE_BUFFER_LENGTH);
        size_t length = LINE_BUFFER_LENGTH;
        ssize_t bytes_read = getline(&bwfn, &length, bw_list_fp);
        int file_idx = 0;
        struct stat fstat;
        mate2len file2size;
        strvec files;
        std::vector<uint64_t> fsizes;
        uint64_t total_fsize = 0;
        uint32_t num_files = 0;
        while(bytes_read != -1) {
            char *bp = bwfn;
            bp[bytes_read-1]='\0';
            int thread_i = file_idx++ % nthreads;
            std::string str(bp);
            files.push_back(str);
            if(LOAD_BALANCE) {
                stat(bp, &fstat);
                fsizes.push_back(fstat.st_size);
                total_fsize += fstat.st_size;
            }
            num_files++;
            bytes_read = getline(&bwfn, &length, bw_list_fp);
        }
        //now load balance between threads based on file size
        uint64_t per_thread_size = total_fsize / nthreads;
        int max_num_files_per_thread = num_files / nthreads;
        int thread_i = 0;
        int num_files_current_thread = 0;
        for(int i=0; i < num_files; i++) {
            if((LOAD_BALANCE && bytes_per_thread[thread_i] + fsizes[i] > per_thread_size) ||
                (num_files_current_thread >= max_num_files_per_thread) && thread_i+1 < nthreads) {
                thread_i++;
                num_files_current_thread = 0;
            }
            if(LOAD_BALANCE)
                bytes_per_thread[thread_i] += fsizes[i];
            files_per_thread[thread_i]->push_back(files[i]);
            num_files_current_thread++;
        }
        std::vector<std::thread> threads;
        for(int i=0; i < nthreads; i++) {
                threads.push_back(std::thread(process_bigwig_worker<T>, std::ref(*(files_per_thread[i])), annotations, chrm_order, keep_order_idx, op=op));
        }
        for(auto &t: threads) t.join();
        fclose(bw_list_fp);
        if(afp && afp != stdout)
            fclose(afp);
        if(afpz)
            bgzf_close(afpz);
        std::free(bwfn);
        return 0;
    }
    //don't have a list of BigWigs, so just process the single one
    int ret = process_bigwig(bw_arg, &annotated_total_auc, annotations, annotation_chrs_seen, afp, keep_order_idx, op=op);
    //if we wanted to keep the chromosome order of the annotation output matching the input BED file
    if(keep_order)
        output_all_coverage_ordered_by_BED(chrm_order, annotations, afp, afpz, nullptr, nullptr, op);
    else
        output_missing_annotations(annotations, annotation_chrs_seen, afp, op = op);
    if(afp && afp != stdout)
        fclose(afp);
    if(ret == 0 && auc_file)
        fprintf(auc_file, "AUC_ANNOTATED_BASES\t%.3f\n", annotated_total_auc);
    if(auc_file && auc_file != stdout)
        fclose(auc_file);
    return ret;
}

int sam_index_iterator_wrapper(bam1_t* b, htsFile* bfh, bam_hdr_t* bhdr, hts_itr_t* sam_itr) {
    return sam_itr_next(bfh, sam_itr, b);
}

int sam_scan_iterator_wrapper(bam1_t* b, htsFile* bfh, bam_hdr_t* bhdr, hts_itr_t* sam_itr) {
    return sam_read1(bfh, bhdr, b);
}

int NUM_CHARS_IN_REGION_STR = 1000;
//based on http://www.cplusplus.com/reference/iterator/iterator/
template <typename T>
class BAMIterator : public std::iterator<std::input_iterator_tag, bam1_t>
{
    bam1_t* b;
    htsFile* bfh;
    bam_hdr_t* bhdr;
    hts_idx_t* bidx;
    hts_itr_t* sam_itr;
    int (*itrPtr)(bam1_t* b, htsFile* bfh, bam_hdr_t* bhdr, hts_itr_t* sam_itr) = &sam_scan_iterator_wrapper;
    char* amap;
    char** amap_ptr;

public:
    BAMIterator(bam1_t* z, htsFile* bam_fh, bam_hdr_t* bam_hdr) :b(z),bfh(bam_fh),bhdr(bam_hdr),bidx(nullptr),sam_itr(nullptr) {}
    BAMIterator(bam1_t* z, htsFile* bam_fh, bam_hdr_t* bam_hdr, const char* bam_fn, annotation_map_t<T>* annotations, uint32_t annotations_count, strlist* chrm_order) :b(z),bfh(bam_fh),bhdr(bam_hdr),bidx(nullptr),sam_itr(nullptr) {
        if(annotations_count == 0)
            return;
        //given a set of regions, check to see if we have an accompaning BAM index file (.bai)
        //check if BAI exists, if not proceed with linear scan through BAM iterator
        if((bidx = sam_index_load(bfh, bam_fn)) == 0) {
            fprintf(stderr,"no index for BAM/CRAM file, doing full scan\n");
            return;
        }
        uint32_t amap_count = annotations_count;
        amap = new char[amap_count*NUM_CHARS_IN_REGION_STR];
        amap_ptr = new char*[amap_count];
        uint64_t k = 0;
        for(auto const c : *chrm_order) {
            std::vector<T*>& annotations_for_chr = (*annotations)[c];
            for(long z = 0; z < annotations_for_chr.size(); z++) {
                const auto &item = annotations_for_chr[z];
                const T start = item[0], end = item[1];
                //keep the auto null char
                amap_ptr[k++] = amap;
                amap += (sprintf(amap, "%s:%lu-%lu", c, (long) start, (long) end)+1);
            }
        }
        assert(k==amap_count);
        sam_itr = sam_itr_regarray(bidx, bhdr, amap_ptr, amap_count);
        if(!sam_itr) {
            fprintf(stderr,"failed to create SAM file iterator, exiting\n");
            exit(-1);
        }
        //delete amap;
        //delete amap_ptr;
        itrPtr = &sam_index_iterator_wrapper;
    }
    BAMIterator(const BAMIterator& bitr) : b(bitr.b),bfh(bitr.bfh),bhdr(bitr.bhdr),bidx(bitr.bidx),sam_itr(bitr.sam_itr),itrPtr(bitr.itrPtr) {}

    BAMIterator& operator++() {
        int r = itrPtr(b, bfh, bhdr, sam_itr);
        if(r < 0)
            b = nullptr;
        return *this;
    }

    BAMIterator operator++(int) {BAMIterator temp(*this); operator++(); return temp;}
    bool operator==(const BAMIterator& rhs) const {return b==rhs.b;}
    bool operator!=(const BAMIterator& rhs) const {return b!=rhs.b;}
    bam1_t* operator*() {return b;}
    //~BAMIterator() { if(sam_itr) { hts_itr_destroy(sam_itr); delete amap; delete amap_ptr;} }
    ~BAMIterator() { if(sam_itr) { hts_itr_destroy(sam_itr);} }
};

int finalize_tabix_index(const char* fname, const char* ifname, BGZF* bfh, hts_idx_t* cidx, int* chrms_in_cidx, const bam_hdr_t *hdr) {
    //this function assumes that the chromosome (chrm) order indexes have been tracked while adding
    //intervals to the BGZip file we're finalizing the index for here
    //but now we need to create the final array of chrm order indexes mapped to chrm names
    //while tracking the total length of all chrm names catted together
    //this will serve as part of the index metadata
    //
    //1) track order by which chromosomes are added to index
    //2) create index2chromosome name
    //3) track total chromosome name length (concatenated overall names including the separating '\0's)
    if(hts_idx_finish(cidx, bgzf_tell(bfh)) != 0) {
        fprintf(stderr,"Error finishing BGZF index for base coverage, skipping\n");
        return -1;
    }
    //largely lifted from: https://github.com/samtools/htslib/blob/4162046b28a7d9d8a104ce28086d9467cc05c212/tbx.c#L216
    tbx_t *tbx;
    tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
    tbx->conf = tbx_conf_bed; 
    tbx->idx = cidx;
    //first slot is the number of chromosomes present in the index
    int num_chrms = chrms_in_cidx[0]; 
    int i, all_cnames_len = 0, l_nm; 
    uint32_t x[7];
    memcpy(x, &tbx->conf, 24);
    char** name = new char*[num_chrms];
    int k = 0;
    for(i=1; i < hdr->n_targets+1; i++) {
        if(chrms_in_cidx[i] > 0) {
            all_cnames_len += strlen(hdr->target_name[i-1]) + 1; //+1 for '\0'
            //now copy chrm name into names
            name[k++] = hdr->target_name[i-1];
        }
    }
    assert(k==num_chrms);
    i = 0;
    
    l_nm = x[6] = all_cnames_len;
    
    uint8_t* meta = new uint8_t[l_nm + 28]; 

    if (ed_is_big())
        for (i = 0; i < 7; ++i)
            x[i] = ed_swap_4(x[i]);
    memcpy(meta, x, 28);
    int l = 0;
    for (l = 28, i = 0; i < num_chrms; ++i) {
        int xi = strlen(name[i]) + 1;
        memcpy(meta + l, name[i], xi);
        l += xi;
    }
    //delete name;
    hts_idx_set_meta(tbx->idx, l, meta, 0);

    if(hts_idx_save_as(cidx, fname, ifname, HTS_FMT_CSI) != 0) {
        fprintf(stderr,"Error saving BGZF index for base coverage, skipping\n");
        return -1;
    }
    return 0;
}


template <typename T>
int go_bam(const char* bam_arg, int argc, const char** argv, Op op, htsFile *bam_fh, int nthreads, bool keep_order, bool has_annotation, FILE* afp, BGZF* afpz, annotation_map_t<T>* annotations, chr2bool* annotation_chrs_seen, const char* prefix, bool sum_annotation, strlist* chrm_order, FILE* auc_file, uint64_t num_annotations, uint32_t window_size = 0) {
    //only calculate AUC across either the BAM or the BigWig, but could be restricting to an annotation as well
    uint64_t all_auc = 0;
    uint64_t unique_auc = 0;
    uint64_t annotated_auc = 0;
    uint64_t unique_annotated_auc = 0;

    std::cerr << "Processing BAM: \"" << bam_arg << "\"" << std::endl;

    bam_hdr_t *hdr = sam_hdr_read(bam_fh);
    if(!hdr) {
        std::cerr << "ERROR: Could not read header for " << bam_arg
                  << ": " << std::strerror(errno) << std::endl;
        return -1;
    }
    if(has_option(argv, argv+argc, "--head")) {
        print_header(hdr);
    }
    hts_set_threads(bam_fh, nthreads);


    //setup list of callbacks for the process_cigar()
    //this is so we only have to walk the cigar for each alignment ~1 time
    callback_list process_cigar_callbacks;
    args_list process_cigar_output_args;

    args_list maplen_outlist;
    uint64_t total_number_bases_processed = 0;
    maplen_outlist.push_back(&total_number_bases_processed);
    bool count_bases = has_option(argv, argv+argc, "--num-bases");
    if(count_bases) {
        process_cigar_callbacks.push_back(maplength);
        process_cigar_output_args.push_back(&maplen_outlist);
    }

    bool print_qual = has_option(argv, argv+argc, "--print-qual");
    bool include_sc = false;
    FILE* softclip_file = nullptr;
    uint64_t total_softclip_count = 0;
    uint64_t total_number_sequence_bases_processed = 0;
    if(has_option(argv, argv+argc, "--include-softclip")) {
        include_sc = true;
        char afn[1024];
        sprintf(afn, "%s.softclip.tsv", prefix);
        softclip_file = fopen(afn, "w");
    }
    const bool only_polya_sc = has_option(argv, argv+argc, "--only-polya");
    const bool include_n_mms = has_option(argv, argv+argc, "--include-n");
    const bool double_count = has_option(argv, argv+argc, "--double-count");
    const bool report_end_coord = has_option(argv, argv+argc, "--ends");
    if(has_option(argv, argv+argc, "--test-polya")) {
        SOFTCLIP_POLYA_TOTAL_COUNT_MIN=1;
        SOFTCLIP_POLYA_RATIO_MIN=0.01;
    }

    size_t recs = 0;
    std::vector<MdzOp> mdzbuf;
    bam1_t *rec = bam_init1();
    if(!rec) {
        std::cerr << "ERROR: Could not initialize BAM object: "
                  << std::strerror(errno) << std::endl;
        return -1;
    }
    kstring_t sambuf{ 0, 0, nullptr };
    bool first = true;
    //largest human chromosome is ~249M bases
    //long chr_size = 250000000;
    long chr_size = -1;
    std::unique_ptr<uint32_t[]> coverages, unique_coverages;
    bool compute_coverage = false;
    int bw_unique_min_qual = 0;
    read2len overlapping_mates;
    bigWigFile_t *bwfp = nullptr;
    bigWigFile_t *ubwfp = nullptr;
    //--coverage -> output perbase coverage to STDOUT (compute_coverage=true)
    //--bigwig -> output perbase coverage to bigwig (compute_coverage=true),
    //  this option overrides --coverage=>coverage will be *only* written to the bigwig
    //  even if --coverage is also passed in
    //--auc -> output AUC of coverage (compute_coverage=true)
    //--annotation output annotated regions of coverage (compute_coverage=true)
    bool auc_opt = has_option(argv, argv+argc, "--auc") || argc == 1;
    bool coverage_opt = has_option(argv, argv+argc, "--coverage");
    bool annotation_opt = has_option(argv, argv+argc, "--annotation");
    bool bigwig_opt = has_option(argv, argv+argc, "--bigwig");
#ifdef WINDOWS_MINGW
    if(bigwig_opt) {
        bigwig_opt = false;
        fprintf(stderr,"WARNING: writing BigWigs (--bigwig) is not supported on Windows at this time, no BigWig file(s) will be written, but any other options will still be processed.\n");
    }
#endif
    bool dont_output_coverage = !(coverage_opt || bigwig_opt);
    FILE* cov_fh = stdout;
    bool gzip = has_option(argv, argv+argc, "--gzip");
    bool no_coverage_stdout = gzip || has_option(argv, argv+argc, "--no-coverage-stdout");
    //gzFile gcov_fh;
    BGZF* gcov_fh = nullptr;
    hts_idx_t* cidx = nullptr;
    
    bool unique = has_option(argv, argv+argc, "--min-unique-qual");
    FILE* uafp = nullptr;
    BGZF* uafpz = nullptr;
    if(coverage_opt || auc_opt || annotation_opt || bigwig_opt) {
        compute_coverage = true;
        chr_size = get_longest_target_size(hdr);
        coverages.reset(new uint32_t[chr_size]);
        if(bigwig_opt)
            bwfp = create_bigwig_file(hdr, prefix,"all.bw");
        if(unique) {
            if(annotation_opt && window_size == 0) {
                uafp = stdout;
                if(gzip || has_option(argv, argv+argc, "--no-annotation-stdout")) {
                    char afn[1024];
                    if(gzip) {
                        sprintf(afn, "%s.unique.tsv.gz", prefix);
                        uafpz = bgzf_open(afn,"w10");
                        uafp = nullptr;
                    }
                    else {
                        sprintf(afn, "%s.unique.tsv", prefix);
                        uafp = fopen(afn, "w");
                    }
                }
            }
            if(bigwig_opt)
                ubwfp = create_bigwig_file(hdr, prefix, "unique.bw");
            bw_unique_min_qual = atoi(*(get_option(argv, argv+argc, "--min-unique-qual")));
            unique_coverages.reset(new uint32_t[chr_size]);
        }
        if(coverage_opt && !bigwig_opt && no_coverage_stdout) {
            char cov_fn[1024];
            if(gzip) {
                sprintf(cov_fn, "%s.coverage.tsv.gz", prefix);
                gcov_fh = bgzf_open(cov_fn,"w10");
                cov_fh = nullptr;
                //from https://github.com/samtools/htslib/blob/c9175183c42382f1030503e88ca7e60cb9c08536/sam.c#L923
                //and https://github.com/brentp/hts-nim/blob/0eaa867e747d3bc844b5ecb575796e4688b966f5/src/hts/csi.nim#L34
                int min_shift = 14;
                int n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3;
                int fmt = HTS_FMT_CSI;
                cidx = hts_idx_init(0, fmt, 0, min_shift, n_lvls);
            }
            else {
                sprintf(cov_fn, "%s.coverage.tsv", prefix);
                cov_fh = fopen(cov_fn,"w");
            }
        }
    }
    fraglen2count* frag_dist = new fraglen2count(1);
    mate2len* frag_mates = new mate2len(1);
    char cov_prefix[50]="";
    int32_t ptid = -1;
    std::unique_ptr<uint32_t[]> starts, ends;
    bool compute_ends = false;
    FILE* rsfp = nullptr;
    FILE* refp = nullptr;
    if(has_option(argv, argv+argc, "--read-ends")) {
        compute_ends = true;
        char refn[1024];
        sprintf(refn, "%s.starts.tsv", prefix);
        rsfp = fopen(refn,"w");
        sprintf(refn, "%s.ends.tsv", prefix);
        refp = fopen(refn,"w");
        if(chr_size == -1)
            chr_size = get_longest_target_size(hdr);
        starts.reset(new uint32_t[chr_size]);
        ends.reset(new uint32_t[chr_size]);
    }
    bool print_frag_dist = false;
    FILE* fragdist_file = nullptr;
    if(has_option(argv, argv+argc, "--frag-dist")) {
        char afn[1024];
        sprintf(afn, "%s.frags.tsv", prefix);
        fragdist_file = fopen(afn, "w");
        print_frag_dist = true;
    }
    const bool echo_sam = has_option(argv, argv+argc, "--echo-sam");
    std::fstream alts_file;
    bool compute_alts = false;
    if(has_option(argv, argv+argc, "--alts")) {
        char afn[1024];
        sprintf(afn, "%s.alts.tsv", prefix);
        alts_file.open(afn, std::fstream::out);
        compute_alts = true;
    }
    FILE* jxs_file = nullptr;
    bool extract_junctions = false;
    uint32_t len = 0;
    args_list junctions;
    coords jx_coords;
    str2cstr jx_pairs;
    str2int jx_counts;
    if(has_option(argv, argv+argc, "--junctions")) {
        junctions.push_back(&len);
        junctions.push_back(&jx_coords);
        char afn[1024];
        sprintf(afn, "%s.jxs.tsv", prefix);
        jxs_file = fopen(afn, "w");
        extract_junctions = true;
        process_cigar_callbacks.push_back(extract_junction);
        process_cigar_output_args.push_back(&junctions);
    }
    const bool require_mdz = has_option(argv, argv+argc, "--require-mdz");
    //the number of reads we actually looked at (didn't filter)
    uint64_t reads_processed = 0;

    char* cigar_str = new char[10000];

    bool long_reads = false;
    if(has_option(argv, argv+argc, "--long-reads")) {
        long_reads = true;
    }
    int jx_str_sz = 2048;
    if(long_reads)
        //enough for the cigar string and ~100 junctions
        jx_str_sz = 12048;

    //no filter out by default
    int filter_in_mask = 0xFFFFFFFF;
    if(has_option(argv, argv+argc, "--filter-in")) {
        filter_in_mask = atoi(*(get_option(argv, argv+argc, "--filter-in")));
    }
    //filter out alignments with either BAM_FUNMAP and/or BAM_FSECONDARY flags set by default (260)
    int filter_out_mask = 260;
    if(has_option(argv, argv+argc, "--filter-out")) {
        filter_out_mask = atoi(*(get_option(argv, argv+argc, "--filter-out")));
    }
    bam1_t* rec_ = bam_init1();
    uint64_t num_annotations_ = 0;
    if(dont_output_coverage && !auc_opt)
        num_annotations_ = num_annotations;
    int num_cigar_ops = process_cigar_callbacks.size();

    //init to 0's
    int* chrms_in_cidx = new int[hdr->n_targets+1]{};

    //TODO: also implement automatic detection of >=80% region coverage of genome
    //AND automatically turn this on if we're doing windowed regions
    bool skip_index = has_option(argv, argv+argc, "--no-index");
    int num_annotations_for_index = num_annotations;
    if(skip_index)
       num_annotations_for_index = 0; 
    bool no_region = true;
    if(num_annotations > 0)
        no_region = false;
    BAMIterator<T> bitr(rec_, bam_fh, hdr, bam_arg, annotations, num_annotations_for_index, chrm_order);
    BAMIterator<T> end(nullptr, nullptr, nullptr);
    for(++bitr; bitr != end; ++bitr) {
        recs++;
        rec = *bitr;
        bam1_core_t *c = &rec->core;
        //read name
        char* qname = bam_get_qname(rec);
        //fprintf(stderr, "recs %lu, qname %s\n",recs,qname);
        //*******Main Quantification Conditional (for ref & alt coverage, frag dist)
        //filter OUT unmapped and secondary alignments
        //if((c->flag & BAM_FUNMAP) == 0 && (c->flag & BAM_FSECONDARY) == 0) {
        //catch case where c-flag is 0 and we've specified an all inclusive filter-in option (default)
        if(((c->flag & filter_in_mask) != 0 && (c->flag & filter_out_mask) == 0)
                                        || (c->flag == 0 && filter_in_mask == 0xFFFFFFFF)) {
            reads_processed++;
            //base-0 start coordinate
            int32_t refpos = rec->core.pos;
            //size of aligned portion of the read (start to end on the reference)
            uint32_t maplen = -1;
            //base-1 end coordinate
            int32_t end_refpos = -1;
            //base-0 mate start coordinate
            int32_t mrefpos = rec->core.mpos;
            //used for adjusting the fragment lengths
            int32_t total_intron_len = 0;
            //ref chrm/contig ID
            int32_t tid = rec->core.tid;
            int32_t tlen = rec->core.isize;

            if(tid != ptid && ptid != -1)
                chr_size = hdr->target_len[ptid];
            
            if(softclip_file)
                total_number_sequence_bases_processed += c->l_qseq;

            //*******Reference coverage tracking
            if(compute_coverage) {
                if(tid != ptid) {
                    if(ptid != -1) {
                        overlapping_mates.clear();
                        sprintf(cov_prefix, "cov\t%d", ptid);
                        if(coverage_opt || bigwig_opt || auc_opt || window_size > 0) {
                            if(no_region) {
                                all_auc += print_array<int32_t>(cov_prefix, hdr->target_name[ptid], ptid, (int32_t*) coverages.get(), chr_size, false, bwfp, cov_fh, dont_output_coverage, no_region, gcov_fh, cidx, chrms_in_cidx, afp, afpz, window_size, op);
                                if(unique) {
                                    sprintf(cov_prefix, "ucov\t%d", ptid);
                                    unique_auc += print_array<int32_t>(cov_prefix, hdr->target_name[ptid], ptid, (int32_t*) unique_coverages.get(), chr_size, false, ubwfp, cov_fh, dont_output_coverage, no_region);
                                }
                            }
                            else {
                                all_auc += print_array<uint32_t>(cov_prefix, hdr->target_name[ptid], ptid, coverages.get(), chr_size, false, bwfp, cov_fh, dont_output_coverage, no_region, gcov_fh, cidx, chrms_in_cidx, afp, afpz, window_size, op);
                                if(unique) {
                                    sprintf(cov_prefix, "ucov\t%d", ptid);
                                    unique_auc += print_array<uint32_t>(cov_prefix, hdr->target_name[ptid], ptid, unique_coverages.get(), chr_size, false, ubwfp, cov_fh, dont_output_coverage, no_region);
                                }
                            }
                        }
                        //if we also want to sum coverage across a user supplied file of annotated regions
                        int keep_order_idx = keep_order?2:-1;
                        if(sum_annotation && annotations->find(hdr->target_name[ptid]) != annotations->end()) {
                            sum_annotations(coverages.get(), (*annotations)[hdr->target_name[ptid]], chr_size, hdr->target_name[ptid], afp, &annotated_auc, op, !annotation_opt, keep_order_idx);
                            if(unique) {
                                keep_order_idx = keep_order?3:-1;
                                sum_annotations(unique_coverages.get(), (*annotations)[hdr->target_name[ptid]], chr_size, hdr->target_name[ptid], uafp, &unique_annotated_auc, op, !annotation_opt, keep_order_idx);
                            }
                            if(!keep_order)
                                annotation_chrs_seen->insert(hdr->target_name[ptid]);
                        }
                    }
                    //need to reset the array for the *current* chromosome's size, not the past one
                    reset_array(coverages.get(), hdr->target_len[tid]);
                    if(unique)
                        reset_array(unique_coverages.get(), hdr->target_len[tid]);
                }
                end_refpos = calculate_coverage(rec, coverages.get(), unique_coverages.get(), double_count, bw_unique_min_qual, &overlapping_mates, &total_intron_len, no_region);
            }
            //additional counting options which make use of knowing the end coordinate/maplen
            //however, if we're already running calculate_coverage, we don't need to redo this
            if(end_refpos == -1 && (report_end_coord || print_frag_dist))
                end_refpos = calculate_coverage(rec, nullptr, nullptr, double_count, bw_unique_min_qual, nullptr, &total_intron_len, no_region);

            if(report_end_coord)
                fprintf(stdout, "%s\t%d\n", qname, end_refpos);

            //*******Fragment length distribution (per chromosome)
            if(print_frag_dist) {
                //csaw's getPESizes criteria
                //first, don't count read that's got problems
                if((c->flag & BAM_FSECONDARY) == 0 && (c->flag & BAM_FSUPPLEMENTARY) == 0 &&
                        (c->flag & BAM_FPAIRED) != 0 && (c->flag & BAM_FMUNMAP) == 0 &&
                        ((c->flag & BAM_FREAD1) != 0) != ((c->flag & BAM_FREAD2) != 0) && rec->core.tid == rec->core.mtid) {
                    //are we the later mate? if so we calculate the frag length
                    if(frag_mates->find(qname) != frag_mates->end()) {
                        uint64_t both_lens = (*frag_mates)[qname];
                        int32_t both_intron_lengths = total_intron_len + (both_lens & frag_lens_mask);
                        both_lens = both_lens >> FRAG_LEN_BITLEN;
                        int32_t mreflen = (both_lens & frag_lens_mask);
                        frag_mates->erase(qname);
                        if(((c->flag & BAM_FREVERSE) != 0) != ((c->flag & BAM_FMREVERSE) != 0) &&
                                (((c->flag & BAM_FREVERSE) == 0 && refpos < mrefpos + mreflen) || ((c->flag & BAM_FMREVERSE) == 0 && mrefpos < end_refpos))) {
                            if(both_intron_lengths > abs(rec->core.isize))
                                both_intron_lengths = 0;
                            (*frag_dist)[abs(rec->core.isize)-both_intron_lengths]++;
                        }
                    }
                    else {
                        uint64_t both_lens = end_refpos - refpos;
                        both_lens = both_lens << FRAG_LEN_BITLEN;
                        both_lens |= total_intron_len;
                        (*frag_mates)[qname] = both_lens;
                    }
                }
            }

            //*******Start/end positions (for TSS,TES)
            //track read starts/ends
            //if minimum quality is set, then we only track starts/ends for alignments that pass
            if(compute_ends) {
                int32_t refpos = rec->core.pos;
                if(tid != ptid) {
                    if(ptid != -1) {
                        for(uint32_t j = 0; j < chr_size; j++) {
                            if(starts[j] > 0)
                                fprintf(rsfp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, starts[j]);
                            if(ends[j] > 0)
                                fprintf(refp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, ends[j]);
                        }
                    }
                    reset_array(starts.get(), hdr->target_len[tid]);
                    reset_array(ends.get(), hdr->target_len[tid]);
                }
                if(bw_unique_min_qual == 0 || rec->core.qual >= bw_unique_min_qual) {
                    starts[refpos]++;
                    if(end_refpos == -1)
                        end_refpos = refpos + align_length(rec);
                    //offset by 1
                    ends[end_refpos-1]++;
                }
            }
            ptid = tid;

            //echo back the sam record
            if(echo_sam) {
                int ret = sam_format1(hdr, rec, &sambuf);
                if(ret < 0) {
                    std::cerr << "Could not format SAM record: " << std::strerror(errno) << std::endl;
                    return -1;
                }
                kstring_out(std::cout, &sambuf);
                std::cout << '\n';
            }
            //*******Alternate base coverages, soft clipping output
            //track alt. base coverages
            if(compute_alts) {
                if(first) {
                    if(print_qual) {
                        uint8_t *qual = bam_get_qual(rec);
                        if(qual[0] == 255) {
                            std::cerr << "WARNING: --print-qual specified but quality strings don't seem to be present" << std::endl;
                            print_qual = false;
                        }
                    }
                    first = false;
                }
                const uint8_t *mdz = bam_aux_get(rec, "MD");
                if(!mdz) {
                    if(require_mdz) {
                        std::stringstream ss;
                        ss << "No MD:Z extra field for aligned read \"" << hdr->target_name[c->tid] << "\"";
                        throw std::runtime_error(ss.str());
                    }
                    output_from_cigar(rec, alts_file, &total_softclip_count, include_sc, only_polya_sc); // just use CIGAR
                } else {
                    mdzbuf.clear();
                    parse_mdz(mdz + 1, mdzbuf); // skip type character at beginning
                    output_from_cigar_mdz(
                            rec, mdzbuf, alts_file, &total_softclip_count,
                            print_qual, include_sc, only_polya_sc, include_n_mms); // use CIGAR and MD:Z
                }
            }
            //*******Run various cigar-related functions for 1 pass through the cigar string
            if(num_cigar_ops > 0)
                process_cigar(rec->core.n_cigar, bam_get_cigar(rec), &cigar_str, &process_cigar_callbacks, &process_cigar_output_args);

            //*******Extract jx co-occurrences (not all junctions though)
            if(extract_junctions) {
                bool paired = (c->flag & BAM_FPAIRED) != 0;
                int32_t tlen_orig = tlen;
                int32_t mtid = c->mtid;
                if(tid != mtid)
                    tlen = mtid > tid ? 1000 : -1000;
                //output
                coords* cl = (coords*) junctions[1];
                int sz = cl->size();
                char* jx_str = nullptr;
                //first create jx string for any of the normal conditions
                if(sz >= 4 || (paired && sz >= 2)) {
                    jx_str = new char[jx_str_sz];
                    //coordinates are 1-based chromosome
                    int ix = sprintf(jx_str, "%s\t%d\t%d\t%d\t%s\t", hdr->target_name[tid], refpos+1, (c->flag & 16) != 0, tlen_orig, cigar_str);
                    //int ix = sprintf(jx_str, "%s\t%d\t%d\t%d\t", hdr->target_name[tid], refpos+1, (c->flag & 16) != 0, tlen_orig);
                    for(int jx = 0; jx < sz; jx++) {
                        uint32_t coord = refpos + (*cl)[jx];
                        if(jx % 2 == 0) {
                            if(jx >=2 )
                                ix += sprintf(jx_str+ix, ",");
                            ix += sprintf(jx_str+ix, "%d-", coord+1);
                        }
                        else
                            ix += sprintf(jx_str+ix, "%d", coord);
                    }
                }
                //now determine if we're 1st/2nd/single mate
                if(paired) {
                    //first mate
                    if(tlen > 0 && sz >= 2) {
                        jx_pairs[qname] = jx_str;
                        jx_counts[qname] = sz;
                    }
                    //2nd mate
                    else if(tlen < 0) {
                        bool prev_mate_printed = false;
                        //1st mate with > 0 introns
                        int mate_sz = 0;
                        if(jx_pairs.find(qname) != jx_pairs.end()) {
                            char* pre_jx_str = jx_pairs[qname];
                            mate_sz = jx_counts[qname];
                            //there must be at least 2 introns between the mates
                            if(mate_sz >= 4 || (mate_sz >= 2 && sz >= 2)) {
                                fprintf(jxs_file, "%s", pre_jx_str);
                                prev_mate_printed = true;
                            }
                            delete pre_jx_str;
                            jx_pairs.erase(qname);
                            jx_counts.erase(qname);
                        }
                        //2nd mate with > 0 introns
                        if(sz >= 4 || (mate_sz >= 2 && sz >= 2)) {
                            if(prev_mate_printed)
                                fprintf(jxs_file, "\t");
                            fprintf(jxs_file, "%s", jx_str);
                            prev_mate_printed = true;
                        }
                        if(prev_mate_printed)
                            fprintf(jxs_file,"\n");
                        delete jx_str;
                    }
                }
                //not paired, only care if we have 2 or more introns
                else if(sz >= 4) {
                    fprintf(jxs_file, "%s\n", jx_str);
                    delete jx_str;
                }
                //reset for next alignment
                *((uint32_t*) junctions[0]) = 0;
                cl->clear();
            }
        }
    }
    if(ptid != -1)
        chr_size = hdr->target_len[ptid];
    delete(cigar_str);
    if(jxs_file) {
        fclose(jxs_file);
    }
    if(print_frag_dist) {
        if(ptid != -1)
            print_frag_distribution(frag_dist, fragdist_file);
        fclose(fragdist_file);
    }
    if(compute_coverage) {
        if(ptid != -1) {
            sprintf(cov_prefix, "cov\t%d", ptid);
            if(coverage_opt || bigwig_opt || auc_opt || window_size > 0) {
                if(no_region)
                    all_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, (int32_t*) coverages.get(), chr_size, false, bwfp, cov_fh, dont_output_coverage, no_region, gcov_fh, cidx, chrms_in_cidx, afp, afpz, window_size, op);
                else
                    all_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, coverages.get(), chr_size, false, bwfp, cov_fh, dont_output_coverage, no_region, gcov_fh, cidx, chrms_in_cidx, afp, afpz, window_size, op);
                if(coverage_opt || window_size > 0) {
                    //now print out all contigs/chrms in header which had 0 coverage, only do this for the "all reads" coverage
                    char* last_interval_line = new char[1024];
                    int line_len = 0;
                    int ret = 0;
                    int (*printPtr) (void* fh, char* buf, uint32_t buf_len) = &my_write;
                    void* wcfh = afp; 
                    if(!afp) {
                        printPtr = &my_gzwrite;
                        wcfh = afpz; 
                    }
                    uint32_t wi = 0;
                    char* val = new char[10];
                    sprintf(val,"%d",0);
                    if(op == cmean)
                        sprintf(val,"%.2f",0.00);
                    uint32_t wend = 0;
                    for(int ci=0; ci < hdr->n_targets; ci++) {
                        uint32_t chr_len = hdr->target_len[ci];
                        char* chr_name = hdr->target_name[ci];
                        if(chrms_in_cidx[ci+1] == 0) {
                            chrms_in_cidx[ci+1] = ++chrms_in_cidx[0];
                            if(window_size > 0) {
                                for(wi=0; wi < chr_len; wi+=window_size) {
                                    wend = wi+window_size; 
                                    if(wend > chr_len)
                                        wend = chr_len;
                                    line_len = sprintf(last_interval_line, "%s\t%u\t%u\t%s\n", chr_name, wi, wend, val); 
                                    (*printPtr)(wcfh, last_interval_line, line_len);
                                }
                            }
                            if(coverage_opt) {
                                line_len = sprintf(last_interval_line, "%s\t0\t%u\t0\n", chr_name, chr_len); 
                                if(gcov_fh) {
                                    ret = bgzf_write(gcov_fh, last_interval_line, line_len);
                                    if(cidx) {
                                        if(hts_idx_push(cidx, chrms_in_cidx[ci+1]-1, 0, hdr->target_len[ci], bgzf_tell(gcov_fh), 1) < 0) {
                                            fprintf(stderr,"error writing line in index at coordinates: %s:%u-%u, tid: %d idx tid: %d exiting\n", hdr->target_name[ci], 0, hdr->target_len[ci], ci, chrms_in_cidx[ci+1]-1);
                                            exit(-1);
                                        }
                                    }
                                }
                                else
                                    ret = fwrite(last_interval_line, sizeof(char), line_len, cov_fh);
                            }
                        }
                    }
                }
                if(unique) {
                    sprintf(cov_prefix, "ucov\t%d", ptid);
                    if(no_region)
                        unique_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, (int32_t*) unique_coverages.get(), chr_size, false, ubwfp, cov_fh, dont_output_coverage, no_region);
                    else
                        unique_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, unique_coverages.get(), chr_size, false, ubwfp, cov_fh, dont_output_coverage, no_region);
                }
            }
            if(sum_annotation && annotations->find(hdr->target_name[ptid]) != annotations->end()) {
                int keep_order_idx = keep_order?2:-1;
                sum_annotations(coverages.get(), (*annotations)[hdr->target_name[ptid]], chr_size, hdr->target_name[ptid], afp, &annotated_auc, op, false, keep_order_idx);
                if(unique) {
                    keep_order_idx = keep_order?3:-1;
                    sum_annotations(unique_coverages.get(), (*annotations)[hdr->target_name[ptid]], chr_size, hdr->target_name[ptid], uafp, &unique_annotated_auc, op, false, keep_order_idx);
                }
                if(!keep_order)
                    annotation_chrs_seen->insert(hdr->target_name[ptid]);
            }
            //if we wanted to keep the chromosome order of the annotation output matching the input BED file
            //assert(afpz == uafpz || (afpz != nullptr && uafpz != nullptr));
            if(keep_order)
                output_all_coverage_ordered_by_BED(chrm_order, annotations, afp, afpz, uafp, uafpz);
        }
        if(sum_annotation && auc_file) {
            fprintf(auc_file, "ALL_READS_ANNOTATED_BASES\t%" PRIu64 "\n", annotated_auc);
            if(unique)
                fprintf(auc_file, "UNIQUE_READS_ANNOTATED_BASES\t%" PRIu64 "\n", unique_annotated_auc);
        }
        if(sum_annotation && !keep_order) {
            output_missing_annotations(annotations, annotation_chrs_seen, afp);
            if(unique)
                output_missing_annotations(annotations, annotation_chrs_seen, uafp);
        }
        if(auc_file) {
            fprintf(auc_file, "ALL_READS_ALL_BASES\t%" PRIu64 "\n", all_auc);
            if(unique)
                fprintf(auc_file, "UNIQUE_READS_ALL_BASES\t%" PRIu64 "\n", unique_auc);
        }
    }
    if(compute_ends) {
        if(ptid != -1) {
            for(uint32_t j = 0; j < chr_size; j++) {
                if(starts[j] > 0)
                    fprintf(rsfp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, starts[j]);
                if(ends[j] > 0)
                    fprintf(refp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, ends[j]);
            }
        }
    }
    if(bwfp) {
        bwClose(bwfp);
        if(!ubwfp)
            bwCleanup();
    }
    if(ubwfp) {
        bwClose(ubwfp);
        bwCleanup();
    }
    //for writing out an index for BGZipped coverage BED files
    char temp_afn[1024];
    int min_shift = 14;
    tbx_conf_t tconf = tbx_conf_bed;
    if(cov_fh && cov_fh != stdout)
        fclose(cov_fh);
    if(gzip && gcov_fh) {
        sprintf(temp_afn, "%s.coverage.tsv.gz", prefix);
        char temp_afni[1024];
        sprintf(temp_afni, "%s.coverage.tsv.gz.csi", prefix);
        int check = finalize_tabix_index(temp_afn, temp_afni, gcov_fh, cidx, chrms_in_cidx, hdr);
        bgzf_close(gcov_fh);
    }
    if(gzip && afpz) {
        sprintf(temp_afn, "%s.annotation.tsv.gz", prefix);
        if(window_size > 0)
            sprintf(temp_afn, "%s.window.tsv.gz", prefix);
        bgzf_close(afpz);
        if(tbx_index_build(temp_afn, min_shift, &tconf) != 0) {
            fprintf(stderr,"Error dumping BGZF index for annotation coverage (all alignments), skipping\n");
        }
    }
    if(gzip && uafpz) {
        sprintf(temp_afn, "%s.unique.tsv.gz", prefix);
        bgzf_close(uafpz);
        if(tbx_index_build(temp_afn, min_shift, &tconf) != 0) {
            fprintf(stderr,"Error dumping BGZF index for annotation coverage (unique alignments), skipping\n");
        }
    }
    if(rsfp)
        fclose(rsfp);
    if(refp)
        fclose(refp);
    if(compute_alts && alts_file)
        alts_file.close();
    if(auc_file && auc_file != stdout)
        fclose(auc_file);
    if(afp && afp != stdout)
        fclose(afp);
    if(uafp)
        fclose(uafp);
    fprintf(stderr,"Read %" PRIu64 " records\n",recs);
    if(count_bases) {
        fprintf(stdout,"%" PRIu64 " records passed filters\n",reads_processed);
        fprintf(stdout,"%" PRIu64 " bases in alignments which passed filters\n",*((uint64_t*) maplen_outlist[0]));
        //fprintf(stdout,"%lu bases in alignments which passed filters\n",total_number_bases_processed);
    }
    if(softclip_file) {
        fprintf(softclip_file,"%" PRIu64 " bases softclipped\n",total_softclip_count);
        fprintf(softclip_file,"%" PRIu64 " total number of processed sequence bases\n",total_number_sequence_bases_processed);
        fclose(softclip_file);
    }
    fprintf(stderr,"# of overlapping pairs: %" PRIu64 "\n", num_overlapping_pairs);
    return 0;
}

template <typename T>
int go(const char* fname_arg, int argc, const char** argv, Op op, htsFile *bam_fh, bool is_bam) {
    //number of bam decompression threads
    //0 == 1 thread for the whole program,fname_arg//decompression shares a single core with processing
    //This can also indicate the number of parallel threads to process a list of BigWigs for
    //the purpose of summing over a passed in annotation
    int nthreads = 0;
    if(has_option(argv, argv+argc, "--threads")) {
        const char** nthreads_ = get_option(argv, argv+argc, "--threads");
        nthreads = atoi(*nthreads_);
    }
    bool keep_order = !has_option(argv, argv+argc, "--keep-order");
    strlist chrm_order;
    FILE* afp = nullptr;
    annotation_map_t<T> annotations;
    bool sum_annotation = false;
    chr2bool annotation_chrs_seen;
    //setup hashmap to store BED file of *non-overlapping* annotated intervals to sum coverage across
    //maps chromosome to vector of uint arrays storing start/end of annotated intervals
    int err = 0;
    bool has_annotation = has_option(argv, argv+argc, "--annotation");
    bool gzip = has_option(argv, argv+argc, "--gzip");
    bool no_annotation_stdout = has_option(argv, argv+argc, "--no-annotation-stdout");
    const char* prefix = fname_arg;
    uint64_t num_annotations = 0;
    if(has_option(argv, argv+argc, "--prefix"))
            prefix = *(get_option(argv, argv+argc, "--prefix"));
    BGZF* afpz = nullptr;
    uint32_t window_size = 0;
    if(has_annotation) {
        const char* afile = *(get_option(argv, argv+argc, "--annotation"));
        if(!afile) {
            std::cerr << "No argument to --annotation" << std::endl;
            return -1;
        }
        //TODO: parse afile for a window size (e.g. 200) if doing windowed regions
        char* output_prefix = new char[100];
        sprintf(output_prefix, "window");
        window_size = strtol(afile, nullptr, 10);
        if(window_size == 0) {
            afp = fopen(afile, "r");
            err = read_annotation(afp, &annotations, &chrm_order, keep_order, &num_annotations);
            fclose(afp);
            assert(!annotations.empty());
            std::cerr << annotations.size() << " chromosomes for annotated regions read\n";
            sprintf(output_prefix, "annotation");
            sum_annotation = true;
        }
        else
            fprintf(stderr, "computing coverage windows of length %u\n", window_size);

        afp = stdout;
        if(gzip || no_annotation_stdout) {
            char afn[1024];
            if(gzip) {
                sprintf(afn, "%s.%s.tsv.gz", prefix, output_prefix);
                afpz = bgzf_open(afn,"w10");
                afp = nullptr;
            }
            else {
                sprintf(afn, "%s.%s.tsv", prefix, output_prefix);
                afp = fopen(afn, "w");
            }
        }
    }
    //if no args are passed in other than a file (BAM or BW)
    //then just compute the auc
    FILE* auc_file = nullptr;
    //if we 1) have no params OR 2) we have no params but --bwbuffer OR 3) --auc with/wo any other options
    if(argc == 1
            || has_option(argv, argv+argc, "--auc")
            || (argc == 3 && has_option(argv, argv+argc, "--bwbuffer"))) {
        auc_file = stdout;
        if(has_option(argv, argv+argc, "--no-auc-stdout")) {
            char afn[1024];
            sprintf(afn, "%s.auc.tsv", prefix);
            auc_file = fopen(afn, "w");
        }
    }

    assert(err == 0);
    if(is_bam)
        return go_bam(fname_arg, argc, argv, op, bam_fh, nthreads, keep_order, has_annotation, afp, afpz, &annotations, &annotation_chrs_seen, prefix, sum_annotation, &chrm_order, auc_file, num_annotations, window_size = window_size);
    else
        return go_bw(fname_arg, argc, argv, op, bam_fh, nthreads, keep_order, has_annotation, afp, afpz, &annotations, &annotation_chrs_seen, prefix, sum_annotation, &chrm_order, auc_file, num_annotations);
}

int get_file_format_extension(const char* fname) {
    int slen = strlen(fname);
    if(strcmp("bam", &(fname[slen-3])) == 0 || strcmp("sam", &(fname[slen-3])) == 0)
        return BAM_FORMAT;
    if(strcmp("cram", &(fname[slen-4])) == 0)
        return CRAM_FORMAT;
    if(strcmp("bw", &(fname[slen-2])) == 0
            || strcmp("BW", &(fname[slen-2])) == 0
            || strcmp("bigwig", &(fname[slen-6])) == 0
            || strcmp("bigWig", &(fname[slen-6])) == 0
            || strcmp("BigWig", &(fname[slen-6])) == 0)
        return BW_FORMAT;
    return UNKNOWN_FORMAT;
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
    if(has_option(argv, argv+argc, "--bwbuffer")) {
        const char* opstr = *(get_option(argv, argv+argc, "--bwbuffer"));
        BW_READ_BUFFER = atol(opstr);
    }
    if(has_option(argv, argv+argc, "--sums-only")) {
        SUMS_ONLY = true;
    }
    const char *fname_arg = get_positional_n(argv, argv+argc, 0);
    if(!fname_arg) {
        std::cerr << "ERROR: Could not find <bam|bw> positional arg" << std::endl;
        return -1;
    }

    int format_code = get_file_format_extension(fname_arg);
    if(format_code == UNKNOWN_FORMAT) {
        std::cerr << "ERROR: Could determine format of " << fname_arg << " exiting" << std::endl;
        return -1;
    }

    bool is_bam = (format_code == BAM_FORMAT || format_code == CRAM_FORMAT);
    htsFile* bam_fh = nullptr;
    if(is_bam) {
        bam_fh = sam_open(fname_arg, "r");
        if(!bam_fh) {
            std::cerr << "ERROR: Could not open " << fname_arg << ": "
                      << std::strerror(errno) << std::endl;
            return -1;
        }
        const htsFormat* format = hts_get_format(bam_fh);
        const char* hts_format_ex = hts_format_file_extension(format);
        if(CRAM_FORMAT) {
            //from https://github.com/samtools/samtools/pull/299/files
            //and https://github.com/brentp/mosdepth/blob/389ca702c5709654a5d4c1608073d26315ce3e35/mosdepth.nim#L867
            //turn off decoding of unused base qualities and other unused fields for just base coverage
            //but only if --alts isn't passed in
            hts_set_opt(bam_fh, CRAM_OPT_DECODE_MD, 0);
            hts_set_opt(bam_fh, CRAM_OPT_REQUIRED_FIELDS, SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_RNEXT | SAM_PNEXT);
            if(has_option(argv, argv+argc, "--alts")) {
                //we want everything decoded
                hts_set_opt(bam_fh, CRAM_OPT_DECODE_MD, 1);
                hts_set_opt(bam_fh, CRAM_OPT_REQUIRED_FIELDS, SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_MAPQ | SAM_RNEXT | SAM_PNEXT | SAM_TLEN | SAM_QUAL | SAM_AUX | SAM_RGAUX | SAM_SEQ);
            }
            if(has_option(argv, argv+argc, "--fasta")) {
                const char* fasta_file = *(get_option(argv, argv+argc, "--fasta"));
                int ret = hts_set_fai_filename(bam_fh, fasta_file);
                if(ret != 0) {
                    std::cerr << "ERROR: Could not use the passed in FASTA index " << fasta_file << " exiting" << std::endl;
                    return -1;
                }
            }
        }
    }
    Op op = csum;
    if(has_option(argv, argv+argc, "--op")) {
        const char* opstr = *(get_option(argv, argv+argc, "--op"));
        op = get_operation(opstr);
    }
    std::ios::sync_with_stdio(false);
    if(!is_bam || op == cmean)
        return go<double>(fname_arg, argc, argv, op, bam_fh, is_bam);
    else
        return go<long>(fname_arg, argc, argv, op, bam_fh, is_bam);
}
