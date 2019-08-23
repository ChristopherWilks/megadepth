//
// Created by Ben Langmead on 2018-12-12.
// modified by cwilks 2019-01-09.
//

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
#include "bigWig.h"


typedef std::vector<std::string> strvec;
typedef std::unordered_map<std::string, uint64_t> mate2len;
typedef std::unordered_map<std::string, double*> str2dblist;

uint64_t MAX_INT = (2^63);
uint64_t BLOCKS_PER_BW_ITERATION = 1000000000;
//how many intervals to start with for a chromosome in a BigWig file
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

static const void print_version() {
    //fprintf(stderr, "bamcount %s\n", string(BAMCOUNT_VERSION).c_str());
    std::cout << "bamcount " << std::string(BAMCOUNT_VERSION) << std::endl;
}

static const char USAGE[] = "BAM and BigWig utility.\n"
    "\n"
    "Usage:\n"
    "  bamcount <bam|bw> [options]\n"
    "\n"
    "Options:\n"
    "  -h --help            Show this screen.\n"
    "  --version            Show version.\n"
    "  --threads            # of threads to do: BAM decompression OR compute sums over multiple BigWigs in parallel\n"
    "                       if the 2nd is intended then a TXT file listing the paths to the BigWigs to process in parallel\n"
    "                       should be passed in as the main input file instead of a single BigWig file.\n"
    "\n"
    "Extract regions and their counts from a BigWig outputting BED format if a BigWig file is detected as input (exclusive of the other BAM modes):\n"
    "                       Extracts all reads from the passed in BigWig and output as BED format and output to STDOUT,\n"
    "                       this will also report total AUC to STDERR,\n"
    "                       if --just-auc is also passed in, it will *only* report AUC\n"
    "  --annotation <bed> <prefix>     Only output the regions in this BED\n"
    "                                  Uses prefix to name the BED file to output to (similar to BAM processing)\n"
    "  --op <raw, sum[default], mean, min, max>     Statistic to run on the intervals provided by --annotation\n"
    "                                               \"raw\" will simply output the original intervals and their values\n"
    "                                               from the BigWig which overlap with the regions in --annotation\n"
    "\n"
    "Extract basic junction information from the BAM, including co-occurrence\n"
    "  --junctions <prefix> Extract jx coordinates, strand, and anchor length, per read\n"
    "                       Writes to a TSV file <prefix>.jxs.tsv\n"
    "\n"
    "Non-reference summaries:\n"
    "  --alts <prefix>              Print differing from ref per-base coverages\n"
    "                               Writes to a CSV file <prefix>.alts.tsv\n"
    "  --include-softclip <prefix>  Print a record to the alts CSV for soft-clipped bases\n"
    "                               Writes total counts to a separate TSV file <prefix>.softclip.tsv\n"
    "  --only-polya                 If --include-softclip, only print softclips which are mostly A's or T's\n"
    "  --include-n                  Print mismatch records when mismatched read base is N\n"
    "  --print-qual                 Print quality values for mismatched bases\n"
    "  --delta                      Print POS field as +/- delta from previous\n"
    "  --require-mdz                Quit with error unless MD:Z field exists everywhere it's\n"
    "                               expected\n"
    "  --no-head                    Don't print sequence names and lengths in header\n"
    "\n"
    "Coverage and quantification:\n"
    "  --coverage           Print per-base coverage (slow but totally worth it)\n"
    "  --auc <prefix>       Print per-base area-under-coverage, will generate it for the genome\n"
    "                       and for the annotation if --annotation is also passed in\n"
    "                       Writes to a TSV file <prefix>.auc.tsv\n"
    "  --bigwig <prefix>    Output coverage as BigWig file(s).  Writes to <prefix>.all.bw\n"
    "                       (also <prefix>.unique.bw when --min-unique-qual is specified).\n"
    "                       Requires libBigWig.\n"
    "  --annotation <bed> <prefix>\n"
    "                       Path to BED file containing list of regions to sum coverage over\n"
    "                       (tab-delimited: chrm,start,end)\n"
    "  --keep-bam-order     Output annotation coverage in the order chromosomes appear in the BAM file.\n"
    "                       The default is to output annotation coverage in the order chromosomes appear in the BED file.\n"
    "  --min-unique-qual <int>\n"
    "                       Output second bigWig consisting built only from alignments\n"
    "                       with at least this mapping quality.  --bigwig must be specified.\n"
    "                       Also produces second set of annotation sums based on this coverage\n"
    "                       if --annotation is enabled\n"
    "  --double-count       Allow overlapping ends of PE read to count twice toward\n"
    "                       coverage\n"
    "  --num-bases          Report total sum of bases in alignments processed (that pass filters)\n"
    "\n"
    "Other outputs:\n"
    "  --read-ends          Print counts of read starts/ends, if --min-unique-qual is set\n"
    "                       then only the alignments that pass that filter will be counted here\n"
    "                       Writes to 2 TSV files: <prefix>.starts.tsv, <prefix>.ends.tsv\n"
    "  --frag-dist <prefix> Print fragment length distribution across the genome\n"
    "                       Writes to a TSV file <prefix>.frags.tsv\n"
    "  --echo-sam           Print a SAM record for each aligned read\n"
    "  --ends               Report end coordinate for each read (useful for debugging)\n"
    "  --test-polya         Lower Poly-A filter minimums for testing (only useful for debugging/testing)\n"
    "\n";

static const char* get_positional_n(const char ** begin, const char ** end, size_t n) {
    size_t i = 0;
    for(const char **itr = begin; itr != end; itr++) {
        if((*itr)[0] != '-') {
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
//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
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

//const char seq_nt16_str[] = "=ACMGRSVTWYHKDBN";
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
    for(int k = 0; k < rec->core.n_cigar; k++) {
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
    for(long i = 0; i < arr_sz; i++)
        arr[i] = 0;
}

//used for buffering up text/gz output
int OUT_BUFF_SZ=4000000;
int COORD_STR_LEN=34;
static uint64_t print_array(const char* prefix, 
                        char* chrm,
                        const uint32_t* arr, 
                        const long arr_sz,
                        const bool skip_zeros,
                        bigWigFile_t* bwfp,
                        const bool just_auc = false) {
    bool first = true;
    bool first_print = true;
    float running_value = 0;
    uint32_t last_pos = 0;
    uint64_t auc = 0;
    //from https://stackoverflow.com/questions/27401388/efficient-gzip-writing-with-gzprintf
    int chrnamelen = strlen(chrm);
    int total_line_len = chrnamelen + COORD_STR_LEN;
    int num_lines_per_buf = round(OUT_BUFF_SZ / total_line_len) - 3;
    //fprintf(stdout, "num_lines_per_buf %d\n", num_lines_per_buf);
    char buf[OUT_BUFF_SZ];
    int buf_written = 0;
    char* bufptr = buf;
    //TODO: speed this up, maybe keep a separate vector
    //which tracks where the runs of the same value stop
    //then only loop through that one w/o the if's
    //
    //OR
    //
    //create modules (functions/classes) which determine
    //the type of output ahead of time to get rid of the if's
    //
    //this will print the coordinates in base-0
    for(uint32_t i = 0; i < arr_sz; i++) {
        if(first || running_value != arr[i]) {
            if(!first) {
                if(running_value > 0 || !skip_zeros) {
                    //based on wiggletools' AUC calculation
                    auc += (i - last_pos) * ((long) running_value);
                    if(not just_auc) {
                        if(bwfp && first_print)
                            bwAddIntervals(bwfp, &chrm, &last_pos, &i, &running_value, 1);
                        else if(bwfp)
                            bwAppendIntervals(bwfp, &last_pos, &i, &running_value, 1);
                        else {
                            if(buf_written >= num_lines_per_buf) {
                                bufptr[0]='\0';
                                fprintf(stdout, "%s", buf); 
                                bufptr = buf;
                                buf_written = 0;
                            }
                            bufptr += sprintf(bufptr, "%s\t%u\t%u\t%.0f\n", chrm, last_pos, i, running_value);
                            buf_written++;
                        } 
                        first_print = false;
                    }
                }
            }
            first = false;
            running_value = arr[i];
            last_pos = i;
        }
    }
    if(!first) {
        if(running_value > 0 || !skip_zeros) {
            auc += (arr_sz - last_pos) * ((long) running_value);
            if(not just_auc) {
                if(bwfp && first_print)
                    bwAddIntervals(bwfp, &chrm, &last_pos, (uint32_t*) &arr_sz, &running_value, 1);
                else if(bwfp)
                    bwAppendIntervals(bwfp, &last_pos, (uint32_t*) &arr_sz, &running_value, 1);
                else {
                    if(buf_written > 0) {
                        bufptr[0]='\0';
                        fprintf(stdout, "%s", buf); 
                    }
                    fprintf(stdout, "%s\t%u\t%lu\t%.0f\n", chrm, last_pos, arr_sz, running_value);
                }
            }
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
    if ((type & 1) && (type & 2)) *((uint32_t*) (*out)[0]) += len;
}

static void end_genomic_coord(const int op, const int len, args_list* out) {
    int type = bam_cigar_type(op);
    if (type & 2) *((uint32_t*) (*out)[0]) += len;
}

static const int32_t align_length(const bam1_t *rec) {
    //bam_endpos processes the whole cigar string
    return bam_endpos(rec) - rec->core.pos;
}

typedef std::unordered_map<std::string, char*> str2cstr;
typedef std::unordered_map<std::string, int> str2int;
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


typedef std::unordered_map<std::string, uint32_t*> read2len;
static const int32_t calculate_coverage(const bam1_t *rec, uint32_t* coverages, 
                                        uint32_t* unique_coverages, const bool double_count, 
                                        const int min_qual, read2len* overlapping_mates,
                                        int32_t* total_intron_length) {
    int32_t refpos = rec->core.pos;
    int32_t mrefpos = rec->core.mpos;
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
    int32_t** mspans = nullptr;
    int mspans_idx = 0;
    std::string tn(qname);
    int32_t end_pos = bam_endpos(rec);
    uint32_t mate_passes_quality = 0;
    //-----First Mate Check
    //if we're the first mate and
    //we're avoiding double counting and we're a proper pair
    //and we overlap with our mate, then store our cigar + length
    //for the later mate to adjust its coverage appropriately
    if(coverages && !double_count && (rec->core.flag & BAM_FPROPER_PAIR) == 2) {
        if(rec->core.tid == rec->core.mtid &&
                end_pos > mrefpos && 
                refpos < mrefpos) {
            const uint32_t* mcigar = bam_get_cigar(rec);
            uint32_t n_cigar = rec->core.n_cigar;
            uint32_t* mate_info = new uint32_t[n_cigar+3];
            mate_info[0] = n_cigar;
            mate_info[1] = refpos;
            mate_info[2] = unique && passing_qual;
            std::memcpy(mate_info+3, mcigar, 4*n_cigar);
            (*overlapping_mates)[tn] = mate_info;
        }
        //-------Second Mate Check
        else if(overlapping_mates->find(tn) != overlapping_mates->end()) {
            uint32_t* mate_info = (*overlapping_mates)[tn];
            uint32_t mn_cigar = mate_info[0];
            int32_t real_mate_pos = mate_info[1];
            mate_passes_quality = mate_info[2];
            uint32_t* mcigar = &(mate_info[3]);
            //bash cigar to get spans of overlap
            int32_t malgn_end_pos = real_mate_pos;
            mspans = new int32_t*[mn_cigar];
            for (k = 0; k < mn_cigar; ++k) {
                const int cigar_op = bam_cigar_op(mcigar[k]);
                if(bam_cigar_type(cigar_op)&2) {
                    const int32_t len = bam_cigar_oplen(mcigar[k]);
                    if(bam_cigar_type(cigar_op)&1) {
                        mspans[mspans_idx] = new int32_t(2);
                        mspans[mspans_idx][0] = malgn_end_pos;
                        mspans[mspans_idx][1] = malgn_end_pos + len;
                        mspans_idx++;
                    }
                    malgn_end_pos += len;
                }
            }
            delete[] mate_info;
            int nerased = overlapping_mates->erase(tn);
            assert(nerased == 1);
            n_mspans = mspans_idx;
            mendpos = malgn_end_pos;
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
                    for(z = algn_end_pos; z < algn_end_pos + len; z++) {
                        coverages[z]++;
                        unique_coverages[z]++;
                    }
                    //now fixup overlapping segment but only if mate passed quality
                    if(n_mspans > 0 && algn_end_pos < mendpos) {
                        //loop until we find the next overlapping span
                        //if are current segment is too early we just keep the span index where it is
                        while(mspans_idx < n_mspans && algn_end_pos >= mspans[mspans_idx][1])
                            mspans_idx++;
                        int32_t cur_end = algn_end_pos + len;
                        int32_t left_end = algn_end_pos;
                        if(left_end < mspans[mspans_idx][0])
                            left_end = mspans[mspans_idx][0];
                        //check 1) we've still got mate spans 2) current segment overlaps the current mate span
                        while(mspans_idx < n_mspans && left_end < mspans[mspans_idx][1] 
                                                    && cur_end > mspans[mspans_idx][0]) {
                            //set right end of segment to decrement
                            int32_t right_end = cur_end;
                            int32_t next_left_end = left_end;
                            if(right_end >= mspans[mspans_idx][1]) {
                                right_end = mspans[mspans_idx][1];
                                //if our segment is greater than the previous mate's
                                //also increment the mate spans index
                                mspans_idx++;
                                if(mspans_idx < n_mspans)
                                    next_left_end = mspans[mspans_idx][0];
                            }
                            else {
                                next_left_end = mspans[mspans_idx][1];
                            }
                            for(z = left_end; z < right_end; z++) {
                                coverages[z]--;
                                if(mate_passes_quality)
                                    unique_coverages[z]--;
                            }
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
                    for(z = algn_end_pos; z < algn_end_pos + len; z++) {
                        coverages[z]++;
                    }
                    //now fixup overlapping segment
                    if(n_mspans > 0 && algn_end_pos < mendpos) {
                        //loop until we find the next overlapping span
                        //if are current segment is too early we just keep the span index where it is
                        while(mspans_idx < n_mspans && algn_end_pos >= mspans[mspans_idx][1])
                            mspans_idx++;
                        int32_t cur_end = algn_end_pos + len;
                        int32_t left_end = algn_end_pos;
                        if(left_end < mspans[mspans_idx][0])
                            left_end = mspans[mspans_idx][0];
                        //check 1) we've still got mate spans 2) current segment overlaps the current mate span
                        while(mspans_idx < n_mspans && left_end < mspans[mspans_idx][1] 
                                                    && cur_end > mspans[mspans_idx][0]) {
                            //set right end of segment to decrement
                            int32_t right_end = cur_end;
                            int32_t next_left_end = left_end;
                            if(right_end >= mspans[mspans_idx][1]) {
                                right_end = mspans[mspans_idx][1];
                                //if our segment is greater than the previous mate's
                                //also increment the mate spans index
                                //delete[] mspans[mspans_idx];
                                mspans_idx++;
                                if(mspans_idx < n_mspans)
                                    next_left_end = mspans[mspans_idx][0];
                            }
                            else {
                                next_left_end = mspans[mspans_idx][1];
                            }
                            for(z = left_end; z < right_end; z++) {
                                coverages[z]--;
                            }
                            left_end = next_left_end;
                        }
                    }    
                }
                algn_end_pos += len;
            }
        }
    }
    if(mspans) {
        for(k = 0; k < n_mspans; k++)
            delete[] mspans[k];
        delete[] mspans;
    }
    return algn_end_pos;
}

typedef std::unordered_map<std::string, std::vector<double*>*> annotation_map_t;
typedef std::vector<char*> strlist;
//about 3x faster than the sstring/string::getline version
static const int process_region_line(char* line, const char* delim, annotation_map_t* amap, strlist* chrm_order, bool keep_order) {
	char* line_copy = strdup(line);
	char* tok = strtok(line_copy, delim);
	int i = 0;
	char* chrm = nullptr;
    long start = -1;
    long end = -1;
    int ret = 0;
	int last_col = END_COL;
	while(tok != nullptr) {
		if(i > last_col)
			break;
		if(i == CHRM_COL) {
			chrm = strdup(tok);
		}
		if(i == START_COL)
			start = atol(tok);
		if(i == END_COL)
			end = atol(tok);
		i++;
		tok = strtok(nullptr, delim);
	}
    //if we need to keep the order, then we'll store values here
    int alen = keep_order?4:2;
    double* coords = new double[alen];
    coords[0] = start;
    coords[1] = end;
    if(alen >= 4) {
        coords[2] = 0.0;
        coords[3] = 0.0;
    }
    if(amap->find(chrm) == amap->end()) {
        auto* v = new std::vector<double*>();
        (*amap)[chrm] = v;
        chrm_order->push_back(chrm);
    }
    (*amap)[chrm]->push_back(coords);
	if(line_copy)
		free(line_copy);
	if(line)
		free(line);
    return ret;
}
    
static const int read_annotation(FILE* fin, annotation_map_t* amap, strlist* chrm_order, bool keep_order) {
	char* line = new char[LINE_BUFFER_LENGTH];
	size_t length = LINE_BUFFER_LENGTH;
	assert(fin);
	ssize_t bytes_read = getline(&line, &length, fin);
	int err = 0;
	while(bytes_read != -1) {
	    err = process_region_line(strdup(line), "\t", amap, chrm_order, keep_order);
        assert(err==0);
		bytes_read = getline(&line, &length, fin);
    }
	std::cerr << "building whole annotation region map done\n";
    return err;
}

static void sum_annotations(const uint32_t* coverages, const std::vector<double*>* annotations, const long chr_size, const char* chrm, FILE* ofp, uint64_t* annotated_auc, bool just_auc = false, int keep_order_idx = -1) {
    long z, j;
    for(z = 0; z < annotations->size(); z++) {
        long sum = 0;
        long start = (*annotations)[z][0];
        long end = (*annotations)[z][1];
        for(j = start; j < end; j++) {
            assert(j < chr_size);
            sum += (uint32_t) coverages[j];
        }
        (*annotated_auc) = (*annotated_auc) + sum;
        if(!just_auc) {
            if(keep_order_idx == -1)
                fprintf(ofp, "%s\t%lu\t%lu\t%lu\n", chrm, start, end, sum);
            else
                (*annotations)[z][keep_order_idx] = sum;
        }
    }
}


static bigWigFile_t* create_bigwig_file(const bam_hdr_t *hdr, const char* out_fn, const char *suffix) {
    if(bwInit(1<<BIGWIG_INIT_VAL) != 0) {
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
typedef std::unordered_map<int32_t, uint32_t> fraglen2count;
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
    fprintf(outfn, "STAT\tCOUNT\t%lu\n", count);
    fprintf(outfn, "STAT\tMEAN_LENGTH\t%.3f\n", mean);
    fprintf(outfn, "STAT\tMODE_LENGTH\t%lu\n", mode);
    fprintf(outfn, "STAT\tMODE_LENGTH_COUNT\t%lu\n", mode_count);
    fprintf(outfn, "STAT\tKALLISTO_COUNT\t%lu\n", kcount);
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


typedef std::unordered_map<std::string, bool> chr2bool;
enum Op { craw, csum, cmean, cmin, cmax };
typedef std::unordered_map<std::string, int> str2op;
static int process_bigwig(const char* fn, uint64_t* all_auc, uint64_t* annotated_auc, annotation_map_t* amap, chr2bool* annotation_chrs_seen, FILE* afp, bool just_auc, int keep_order_idx = -1, Op op = csum, FILE* errfp = stderr, str2dblist* store_local=nullptr) {
    //in part lifted from https://github.com/dpryan79/libBigWig/blob/master/test/testIterator.c
    if(bwInit(1<<17) != 0) {
        fprintf(errfp, "Error in bwInit, exiting\n");
        return 1;
    }
    char* fn_copy = strdup(fn);
    bigWigFile_t *fp = bwOpen(fn_copy, NULL, "r");
    if(!fp) {
        fprintf(errfp, "Error in opening %s as BigWig file, exiting\n", fn);
        return 1;
    }
    uint32_t i, tid, blocksPerIteration;
    //ask for huge # of blocks to ensure we get all in one go
    blocksPerIteration = BLOCKS_PER_BW_ITERATION;
    //blocksPerIteration = 1;
    bwOverlapIterator_t *iter = nullptr;
    //raw mode we only want to process BW intervals once
    std::vector<bool> intervals_seen(STARTING_NUM_INTERVALS,false);
    //loop through all the chromosomes in the BW
    for(tid = 0; tid < fp->cl->nKeys; tid++)
    {
        //only process the chromosome if it's in the annotation 
        if(amap->find(fp->cl->chrom[tid]) != amap->end()) {
            iter = bwOverlappingIntervalsIterator(fp, fp->cl->chrom[tid], 0, fp->cl->len[tid], blocksPerIteration);
            if(!iter->data)
            {
                fprintf(errfp, "WARNING: no intervals for chromosome %s in %s as BigWig file, skipping\n", fp->cl->chrom[tid], fn);
                continue;
            }
            uint32_t num_intervals = iter->intervals->l;
            if(op == csum) {
                if(num_intervals > intervals_seen.size())
                    intervals_seen.resize(num_intervals, false);
                for(int i = 0; i < intervals_seen.size(); i++)
                    intervals_seen[i] = false;
            }
            uint32_t istart = iter->intervals->start[0];
            uint32_t iend = iter->intervals->end[num_intervals-1];
            std::vector<double*>* annotations = (*amap)[fp->cl->chrom[tid]];
            long z, j, k;
            long last_j = 0;
            long asz = annotations->size();
            double* local_vals;
            //if running in multithreaded mode, want to store the values locally
            //but also don't want to reallocate for every new bigwig file, so
            //we allocate once per thread per chromosome
            if(store_local) {
                if(store_local->find(fp->cl->chrom[tid]) == store_local->end())
                    local_vals = new double[asz];
                else
                    local_vals = (*store_local)[fp->cl->chrom[tid]];
                for(int i = 0; i < asz; i++)
                    local_vals[i] = 0.0;
            }
            //loop through annotation intervals as outer loop
            for(z = 0; z < asz; z++) {
                double sum = 0;
                double min = MAX_INT;
                double max = 0;
                long start = (*annotations)[z][0];
                long ostart = start;
                long end = (*annotations)[z][1];
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
                        //raw mode
                        if(op == craw) {
                            k = iend;
                            //only want to output intervals once
                            if(intervals_seen[j])
                                continue;
                            intervals_seen[j] = true;
                            //don't try to keep order, just output the intervals as they are in the bigwig
                            fprintf(afp, "%s\t%u\t%u\t%.3f\n", fp->cl->chrom[tid], istart, iend, iter->intervals->value[j]);
                            continue;
                        }
                        //stat mode
                        else {
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
                    (*annotated_auc) = (*annotated_auc) + ((uint32_t) sum);
                if(!just_auc || op != csum) {
                    //0-based start
                    double annot_length = end - ostart;
                    double value = sum;
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
                    }
                    //not trying to keep the order in the BED file, just print them as we find them
                    if(keep_order_idx == -1)
                        fprintf(afp, "%s\t%lu\t%lu\t%.0f\n", fp->cl->chrom[tid], ostart, end, value);
                    else
                        if(store_local)
                            local_vals[z] = value;
                        else
                            (*annotations)[z][keep_order_idx] = value;
                }
            }
            (*annotation_chrs_seen)[fp->cl->chrom[tid]] = true;
            if(store_local)
                (*store_local)[fp->cl->chrom[tid]] = local_vals;
            bwIteratorDestroy(iter);
        }
    }

    bwClose(fp);
    bwCleanup();
    if(fn_copy);
        delete(fn_copy);
    return 0;
}

void print_long_local(FILE* afp,const char* c, long start, long end, double val, double* local_vals, long z) {
        fprintf(afp, "%s\t%lu\t%lu\t%lu\n", c, start, end, (long) local_vals[z]);
}

void print_long_shared(FILE* afp,const char* c, long start, long end, double val, double* local_vals, long z) {
        fprintf(afp, "%s\t%lu\t%lu\t%lu\n", c, start, end, (long) val);
}

void print_double_local(FILE* afp, const char* c, long start, long end, double val, double* local_vals, long z) {
        fprintf(afp, "%s\t%lu\t%lu\t%.3f\n", c, start, end, local_vals[z]);
}

void print_double_shared(FILE* afp, const char* c, long start, long end, double val, double* local_vals, long z) {
        fprintf(afp, "%s\t%lu\t%lu\t%.3f\n", c, start, end, val);
}

static void output_missing_annotations(const annotation_map_t* annotations, const std::unordered_map<std::string, bool>* annotations_seen, FILE* ofp, Op op = csum) {
    //check if we're doing means output doubles, otherwise output longs
    void (*printPtr) (FILE*, const char*, long, long, double, double*, long) = op == cmean?&print_double_shared:&print_long_shared;
    for(auto const& kv : *annotations) {
        if(annotations_seen->find(kv.first) == annotations_seen->end()) {
            std::vector<double*>* annotations_for_chr = kv.second;
            for(long z = 0; z < annotations_for_chr->size(); z++) {
                long start = (*annotations_for_chr)[z][0];
                long end = (*annotations_for_chr)[z][1];
                (*printPtr)(ofp, kv.first.c_str(), start, end, 0.000, nullptr, z);
            }
        }
    }
}

void output_all_coverage_ordered_by_BED(const strlist* chrm_order, annotation_map_t* annotations, FILE* afp, FILE* uafp, Op op = csum, str2dblist* store_local = nullptr) {
    double* local_vals = nullptr;
    for(auto const c : *chrm_order) {
        if(!c)
            continue;
        std::vector<double*>* annotations_for_chr = (*annotations)[c];
        void (*printPtr_long) (FILE*, const char*, long, long, double, double*, long) = &print_long_shared;
        void (*printPtr_double) (FILE*, const char*, long, long, double, double*, long) = &print_double_shared;
        if(store_local) {
            local_vals = (*store_local)[c];
            printPtr_long = &print_long_local;
            printPtr_double = &print_double_local;
        }
        //check if we're doing means output doubles, otherwise output longs
        void (*printPtr) (FILE*, const char*, long, long, double, double*, long) = op == cmean?printPtr_double:printPtr_long;
        for(long z = 0; z < annotations_for_chr->size(); z++) {
            long start = (*annotations_for_chr)[z][0];
            long end = (*annotations_for_chr)[z][1];
            double val = (*annotations_for_chr)[z][2];
            (*printPtr)(afp, c, start, end, val, local_vals, z);
            //do uniques if asked to
            if(uafp) {
                val = (*annotations_for_chr)[z][3];
                (*printPtr)(uafp, c, start, end, val, local_vals, z);
            }
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
void process_bigwig_worker(strvec& bwfns, annotation_map_t* annotations, strlist* chrm_order, bool just_auc, int keep_order_idx, Op op) {
    //want to just get the filename itself, no path
    str2dblist store_local;
    char* bwfn_copy;
    for(auto bwfn_ : bwfns) {
        strvec tokens;
        const char* bwfn = bwfn_.c_str();
        fprintf(stderr, "about to process %s\n", bwfn);
        bwfn_copy = strdup(bwfn);
        std::string str(bwfn_copy);
        split_string(str, '/', &tokens);
        char afn[1024];
        /*sprintf(afn, "%s.auc.tsv", tokens.back().c_str());
        FILE* aucfp = fopen(afn, "w");*/
        FILE* afp = nullptr;
        sprintf(afn, "%s.err", tokens.back().c_str());
        FILE* errfp = fopen(afn, "w");
        if(!just_auc) {
            sprintf(afn, "%s.all.tsv", tokens.back().c_str());
            afp = fopen(afn, "w");
        }
        chr2bool annotation_chrs_seen;
        uint64_t annotated_auc = 0;

        int ret = process_bigwig(bwfn, nullptr, &annotated_auc, annotations, &annotation_chrs_seen, afp, just_auc, keep_order_idx, op = op, errfp = errfp, &store_local);
        if(ret != 0) {
            fprintf(errfp,"FAILED to process bigwig %s\n", bwfn);
            if(afp)
                fclose(afp);
            fclose(errfp);
            //fclose(aucfp);
            return;
        }
        //if we wanted to keep the chromosome order of the annotation output matching the input BED file
        if(keep_order_idx == 2)
            output_all_coverage_ordered_by_BED(chrm_order, annotations, afp, nullptr, op, &store_local);
        else
            output_missing_annotations(annotations, &annotation_chrs_seen, afp, op = op);
        if(afp)
            fclose(afp);
        //fprintf(aucfp, "AUC\t%" PRIu64 "\n", annotated_auc);
        fprintf(errfp, "AUC\t%" PRIu64 "\n", annotated_auc);
        fprintf(errfp,"SUCCESS processing bigwig %s\n", bwfn);
        fclose(errfp);
        //fclose(aucfp);
        if(bwfn_copy)
            delete(bwfn_copy);
    }
    //hold off on final deletion, for performance
    /*for( auto mitr : store_local)
        delete mitr.second;*/
}

Op get_operation(const char* opstr) {
    if(strcmp(opstr, "raw") == 0)
        return craw;
    if(strcmp(opstr, "mean") == 0)
        return cmean;
    if(strcmp(opstr, "min") == 0)
        return cmin;
    if(strcmp(opstr, "max") == 0)
        return cmax;
    return csum;
}

typedef std::unordered_map<std::string, uint8_t*> str2str;
static const uint64_t frag_lens_mask = 0x00000000FFFFFFFF;
static const int FRAG_LEN_BITLEN = 32;
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
    std::ios::sync_with_stdio(false);


    FILE* auc_file = nullptr;
    bool just_auc = false;
    char afn[1024];
    if(has_option(argv, argv+argc, "--auc")) {
        just_auc = !has_option(argv, argv+argc, "--coverage") && !has_option(argv, argv+argc, "--annotation");
        const char *prefix = *get_option(argv, argv+argc, "--auc");
        sprintf(afn, "%s.auc.tsv", prefix);
        auc_file = fopen(afn, "w");
    }
    //setup hashmap to store BED file of *non-overlapping* annotated intervals to sum coverage across
    //maps chromosome to vector of uint arrays storing start/end of annotated intervals
    FILE* afp = nullptr;
    annotation_map_t annotations; 
    chr2bool annotation_chrs_seen;
    bool sum_annotation = false;
    uint64_t all_auc = 0;
    uint64_t annotated_auc = 0;
    bool keep_order = true;
    strlist chrm_order;
    int err = 0;
    const char* annot_prefix = nullptr;
    if(has_option(argv, argv+argc, "--annotation")) {
        sum_annotation = true;
        const char* afile = *(get_option(argv, argv+argc, "--annotation"));
        if(!afile) {
            std::cerr << "No argument to --annotation" << std::endl;
            return 1;
        }
        annot_prefix = *(get_option(argv, argv+argc, "--annotation", 1));
        if(!annot_prefix) {
            std::cerr << "No argument to --annotation" << std::endl;
            return 1;
        }
        keep_order = !has_option(argv, argv+argc, "--keep-bam-order");
        afp = fopen(afile, "r");
        err = read_annotation(afp, &annotations, &chrm_order, keep_order);
        fclose(afp);
        afp = nullptr;
        if(!just_auc) {
            char afn[1024];
            sprintf(afn, "%s.all.tsv", annot_prefix);
            afp = fopen(afn, "w");
        }
        assert(!annotations.empty());
        std::cerr << annotations.size() << " chromosomes for annotated regions read\n";
    }
    assert(err == 0);
    const char *bam_arg = get_positional_n(argv, argv+argc, 0);
    if(!bam_arg) {
        std::cerr << "ERROR: Could not find <bam|bw> positional arg" << std::endl;
        return 1;
    }

    htsFile *bam_fh = sam_open(bam_arg, "r");
    if(!bam_fh) {
        std::cerr << "ERROR: Could not open " << bam_arg << ": "
                  << std::strerror(errno) << std::endl;
        return 1;
    }
    
    //number of bam decompression threads
    //0 == 1 thread for the whole program, i.e.
    //decompression shares a single core with processing
    //This can also indicate the number of parallel threads to process a list of BigWigs for
    //the purpose of summing over a passed in annotation
    int nthreads = 0;
    if(has_option(argv, argv+argc, "--threads")) {
        const char** nthreads_ = get_option(argv, argv+argc, "--threads");
        nthreads = atoi(*nthreads_);
    }
   
    bool LOAD_BALANCE = false;
    const htsFormat* format = hts_get_format(bam_fh);
    const char* hts_format_ex = hts_format_file_extension(format);
    if(strcmp(hts_format_ex, "bam") != 0)
    {
        int slen = strlen(bam_arg);
        bool is_bw_list_file = strcmp(bam_arg+(slen-3), "txt") == 0;
        if(is_bw_list_file || bwIsBigWig(strdup(bam_arg), nullptr))
        {
            Op op = csum;
            if(has_option(argv, argv+argc, "--op")) {
                const char* opstr = *(get_option(argv, argv+argc, "--op"));
                op = get_operation(opstr);
            }
            std::cerr << "Processing BigWig(s): \"" << bam_arg << "\"" << std::endl;
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
                FILE* bw_list_fp = fopen(bam_arg, "r");
                char* bwfn = new char[LINE_BUFFER_LENGTH];
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
                    bwfn[bytes_read-1]='\0';
                    int thread_i = file_idx++ % nthreads;
                    std::string str(strdup(bwfn));
                    files.push_back(str);
                    if(LOAD_BALANCE) {
                        stat(bwfn, &fstat);
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
                for(int i=0; i < nthreads; i++)
                    threads.push_back(std::thread(process_bigwig_worker, std::ref(*(files_per_thread[i])), &annotations, &chrm_order, just_auc, keep_order_idx, op=op));
                delete(bwfn);
                for(int i=0; i < threads.size(); i++)
                    threads[i].join();
                fclose(bw_list_fp);
                if(afp)
                    fclose(afp);
                fclose(auc_file);
                return 0; 
            }
            //don't have a list of BigWigs, so just process the single one
            int ret = process_bigwig(bam_arg, &all_auc, &annotated_auc, &annotations, &annotation_chrs_seen, afp, just_auc, keep_order_idx, op=op);
            //if we wanted to keep the chromosome order of the annotation output matching the input BED file
            if(keep_order)
                output_all_coverage_ordered_by_BED(&chrm_order, &annotations, afp, nullptr, op);
            else
                output_missing_annotations(&annotations, &annotation_chrs_seen, afp, op = op);
            if(afp)
                fclose(afp);
            if(ret == 0) {
                //fprintf(auc_file, "ALL_READS_ANNOTATED_BASES\t%" PRIu64 "\n", annotated_auc);
                fprintf(stderr, "AUC\t%" PRIu64 "\n", annotated_auc);
                //fclose(auc_file);
            }
            return ret;
        }
        else
        {
            std::cerr << "ERROR: Input file is not a BAM (.bam)/BigWig (.bw)/BigWig list file (.txt) " << bam_arg << std::endl;
            return 1;
        }
    }
    std::cerr << "Processing BAM: \"" << bam_arg << "\"" << std::endl;

    bam_hdr_t *hdr = sam_hdr_read(bam_fh);
    if(!hdr) {
        std::cerr << "ERROR: Could not read header for " << bam_arg
                  << ": " << std::strerror(errno) << std::endl;
        return 1;
    }
    if(!has_option(argv, argv+argc, "--no-head")) {
        print_header(hdr);
    }
    hts_set_threads(bam_fh, nthreads);
    
    //setup list of callbacks for the process_cigar()
    //this is so we only have to walk the cigar for each alignment ~1 time
    callback_list process_cigar_callbacks;
    args_list process_cigar_output_args;

    int filter_out = -1;
    int filter_in = -1;
    std::fstream fq0_file;
    std::fstream fq1_file;
    std::fstream fq2_file;
    bool bam2fastq = false;
    bool one_file = false;
    bool re_reversed = false;
    str2str read_mate;
    if(has_option(argv, argv+argc, "--bam2fastq")) {
        bam2fastq = true;
        re_reversed = has_option(argv, argv+argc, "--re-reverse");
        one_file = has_option(argv, argv+argc, "--one-file");
        const char *prefix = *get_option(argv, argv+argc, "--bam2fastq");
        char afn[1024];
        if(!one_file) {
            sprintf(afn, "%s_1.fastq", prefix);
            fq1_file.open(afn, std::fstream::out);
            sprintf(afn, "%s_2.fastq", prefix);
            fq2_file.open(afn, std::fstream::out);
        }
        sprintf(afn, "%s.fastq", prefix);
        fq0_file.open(afn, std::fstream::out);
        if(has_option(argv, argv+argc, "--filter-out")) {
            const char** filter_out_ = get_option(argv, argv+argc, "--filter-out");
            filter_out = atoi(*filter_out_);
        }
        if(has_option(argv, argv+argc, "--filter-in")) {
            const char** filter_in_ = get_option(argv, argv+argc, "--filter-in");
            filter_in = atoi(*filter_in_);
        }
    }

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
        const char *prefix = *get_option(argv, argv+argc, "--include-softclip");
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
        return 1;
    }
    kstring_t sambuf{ 0, 0, nullptr };
    bool first = true;
    //largest human chromosome is ~249M bases
    //long chr_size = 250000000;
    long chr_size = -1;
    uint32_t* coverages = nullptr;
    uint32_t* unique_coverages = nullptr;
    bool compute_coverage = false;
    bool unique = false;
    int bw_unique_min_qual = 0;
    FILE* uafp = nullptr;
    read2len overlapping_mates;
    bigWigFile_t *bwfp = nullptr;
    bigWigFile_t *ubwfp = nullptr;
    uint64_t unique_auc = 0;
    uint64_t unique_annotated_auc = 0;
    if(has_option(argv, argv+argc, "--coverage") || has_option(argv, argv+argc, "--auc")) {
        compute_coverage = true;
        chr_size = get_longest_target_size(hdr);
        coverages = new uint32_t[chr_size];
        const char *bw_fn = nullptr;
        if(has_option(argv, argv+argc, "--bigwig")) {
            bw_fn = *get_option(argv, argv+argc, "--bigwig");
            bwfp = create_bigwig_file(hdr, bw_fn, "all.bw");
        }
        if(has_option(argv, argv+argc, "--min-unique-qual")) {
            unique = true;
            if(bw_fn)
                ubwfp = create_bigwig_file(hdr, bw_fn, "unique.bw");
            bw_unique_min_qual = atoi(*(get_option(argv, argv+argc, "--min-unique-qual")));
            unique_coverages = new uint32_t[chr_size];
            if(sum_annotation && !just_auc) {
                char afn[1024];
                sprintf(afn, "%s.unique.tsv", annot_prefix);
                uafp = fopen(afn, "w");
            }
        }
    }
    fraglen2count* frag_dist = new fraglen2count(1);
    mate2len* frag_mates = new mate2len(1);
    char prefix[50]="";
    int32_t ptid = -1;
    uint32_t* starts = nullptr;
    uint32_t* ends = nullptr;
    bool compute_ends = false;
    FILE* rsfp = nullptr;
    FILE* refp = nullptr;
    if(has_option(argv, argv+argc, "--read-ends")) {
        compute_ends = true;
        const char *prefix = *get_option(argv, argv+argc, "--read-ends");
        char refn[1024];
        sprintf(refn, "%s.starts.tsv", prefix);
        rsfp = fopen(refn,"w");
        sprintf(refn, "%s.ends.tsv", prefix);
        refp = fopen(refn,"w");
        if(chr_size == -1) 
            chr_size = get_longest_target_size(hdr);
        starts = new uint32_t[chr_size];
        ends = new uint32_t[chr_size];
    }
    bool print_frag_dist = false;
    FILE* fragdist_file = nullptr;
    if(has_option(argv, argv+argc, "--frag-dist")) {
        const char *prefix = *get_option(argv, argv+argc, "--frag-dist");
        char afn[1024];
        sprintf(afn, "%s.frags.tsv", prefix);
        fragdist_file = fopen(afn, "w");
        print_frag_dist = true;
    }
    const bool echo_sam = has_option(argv, argv+argc, "--echo-sam");
    std::fstream alts_file;
    bool compute_alts = false;
    if(has_option(argv, argv+argc, "--alts")) {
        const char *prefix = *get_option(argv, argv+argc, "--alts");
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
        const char *prefix = *get_option(argv, argv+argc, "--junctions");
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

    while(sam_read1(bam_fh, hdr, rec) >= 0) {
        recs++;
        bam1_core_t *c = &rec->core;
        //read name
        char* qname = bam_get_qname(rec);
        //*******Main Quantification Conditional (for ref & alt coverage, frag dist)
        //filter OUT unmapped and secondary alignments
        if((c->flag & BAM_FUNMAP) == 0 && (c->flag & BAM_FSECONDARY) == 0) {
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

            
            if(softclip_file)
                total_number_sequence_bases_processed += c->l_qseq;

            //*******Reference coverage tracking
            if(compute_coverage) {
                if(tid != ptid) {
                    if(ptid != -1) {
                        overlapping_mates.clear();
                        sprintf(prefix, "cov\t%d", ptid);
                        all_auc += print_array(prefix, hdr->target_name[ptid], coverages, hdr->target_len[ptid], false, bwfp, just_auc);
                        if(unique) {
                            sprintf(prefix, "ucov\t%d", ptid);
                            unique_auc += print_array(prefix, hdr->target_name[ptid], unique_coverages, hdr->target_len[ptid], false, ubwfp, just_auc);
                        }
                        //if we also want to sum coverage across a user supplied file of annotated regions
                        int keep_order_idx = keep_order?2:-1;
                        //TODO: streamline hashmap accesses to annotations
                        if(sum_annotation && annotations.find(hdr->target_name[ptid]) != annotations.end()) {
                            sum_annotations(coverages, annotations[hdr->target_name[ptid]], hdr->target_len[ptid], hdr->target_name[ptid], afp, &annotated_auc, just_auc, keep_order_idx);
                            if(unique) {
                                keep_order_idx = keep_order?3:-1;
                                sum_annotations(unique_coverages, annotations[hdr->target_name[ptid]], hdr->target_len[ptid], hdr->target_name[ptid], uafp, &unique_annotated_auc, just_auc, keep_order_idx);
                            }
                            if(!keep_order)
                                annotation_chrs_seen[strdup(hdr->target_name[ptid])] = true;
                        }
                    }
                    reset_array(coverages, chr_size);
                    if(unique)
                        reset_array(unique_coverages, chr_size);
                }
                end_refpos = calculate_coverage(rec, coverages, unique_coverages, double_count, bw_unique_min_qual, &overlapping_mates, &total_intron_len);
            }
            //additional counting options which make use of knowing the end coordinate/maplen
            //however, if we're already running calculate_coverage, we don't need to redo this
            if(end_refpos == -1 && (report_end_coord || print_frag_dist))
                end_refpos = calculate_coverage(rec, nullptr, nullptr, double_count, bw_unique_min_qual, nullptr, &total_intron_len);

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
                        for(uint32_t j = 0; j < hdr->target_len[ptid]; j++) {
                            if(starts[j] > 0)
                                fprintf(rsfp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, starts[j]);
                            if(ends[j] > 0)
                                fprintf(refp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, ends[j]);
                        }
                    }
                    reset_array(starts, chr_size);
                    reset_array(ends, chr_size);
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
                    return 1;
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
            if(process_cigar_callbacks.size() > 0)
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
                    jx_str = new char[2048];
                    //coordinates are 1-based chromosome
                    int ix = sprintf(jx_str, "%s\t%d\t%d\t%d\t%s\t", hdr->target_name[tid], refpos+1, (c->flag & 16) != 0, tlen_orig, cigar_str);
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
            sprintf(prefix, "cov\t%d", ptid);
            all_auc += print_array(prefix, hdr->target_name[ptid], coverages, hdr->target_len[ptid], false, bwfp, just_auc);
            if(unique) {
                sprintf(prefix, "ucov\t%d", ptid);
                unique_auc += print_array(prefix, hdr->target_name[ptid], unique_coverages, hdr->target_len[ptid], false, ubwfp, just_auc);
            }
            if(sum_annotation && annotations.find(hdr->target_name[ptid]) != annotations.end()) {
                int keep_order_idx = keep_order?2:-1;
                sum_annotations(coverages, annotations[hdr->target_name[ptid]], hdr->target_len[ptid], hdr->target_name[ptid], afp, &annotated_auc, just_auc, keep_order_idx);
                if(unique) {
                    keep_order_idx = keep_order?3:-1;
                    sum_annotations(unique_coverages, annotations[hdr->target_name[ptid]], hdr->target_len[ptid], hdr->target_name[ptid], uafp, &unique_annotated_auc, just_auc, keep_order_idx);
                }
                if(!keep_order)
                    annotation_chrs_seen[strdup(hdr->target_name[ptid])] = true;
            }
            //if we wanted to keep the chromosome order of the annotation output matching the input BED file
            if(keep_order)
                output_all_coverage_ordered_by_BED(&chrm_order, &annotations, afp, uafp);
        }
        if(sum_annotation && auc_file) {
            fprintf(auc_file, "ALL_READS_ANNOTATED_BASES\t%" PRIu64 "\n", annotated_auc);
            if(unique)
                fprintf(auc_file, "UNIQUE_READS_ANNOTATED_BASES\t%" PRIu64 "\n", unique_annotated_auc);
        }
        delete[] coverages;
        if(unique)
            delete[] unique_coverages;
        if(sum_annotation && !keep_order) {
            output_missing_annotations(&annotations, &annotation_chrs_seen, afp);
            if(unique)
                output_missing_annotations(&annotations, &annotation_chrs_seen, uafp);
        }
        if(auc_file) {
            fprintf(auc_file, "ALL_READS_ALL_BASES\t%" PRIu64 "\n", all_auc);
            if(unique)
                fprintf(auc_file, "UNIQUE_READS_ALL_BASES\t%" PRIu64 "\n", unique_auc);
        }
    }
    if(compute_ends) {
        if(ptid != -1) {
            for(uint32_t j = 0; j < hdr->target_len[ptid]; j++) {
                if(starts[j] > 0)
                    fprintf(rsfp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, starts[j]);
                if(ends[j] > 0)
                    fprintf(refp,"%s\t%d\t%d\n", hdr->target_name[ptid], j+1, ends[j]);
            }
        }
        delete[] starts;
        delete[] ends;
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
    if(rsfp)
        fclose(rsfp);
    if(refp)
        fclose(refp);
    if(compute_alts && alts_file)
        alts_file.close();
    if(auc_file)
        fclose(auc_file);
    if(afp)
        fclose(afp);
    if(uafp)
        fclose(uafp);
    if(bam2fastq) {
        fq0_file.close();
        if(!one_file) {
            fq1_file.close();
            fq2_file.close();
        }
    }
    fprintf(stdout,"Read %lu records\n",recs);
    if(count_bases) {
        fprintf(stdout,"%lu records passed filters\n",reads_processed);
        fprintf(stdout,"%lu bases in alignments which passed filters\n",*((uint64_t*) maplen_outlist[0]));
        //fprintf(stdout,"%lu bases in alignments which passed filters\n",total_number_bases_processed);
    }
    if(softclip_file) {
        fprintf(softclip_file,"%lu bases softclipped\n",total_softclip_count);
        fprintf(softclip_file,"%lu total number of processed sequence bases\n",total_number_sequence_bases_processed);
        fclose(softclip_file);
    }
    return 0;
}
