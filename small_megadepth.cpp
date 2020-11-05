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
#include <zlib.h>
#include <iterator>

#include <sys/stat.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include "bigWig.h"
//#include "tabix_set_meta.h"
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

#if defined(__GNUC__) || defined(__clang__)
#  ifndef unlikely
#    define unlikely(x) __builtin_expect(!!(x), 0)
#  endif
#  ifndef likely
#    define likely(x) __builtin_expect(!!(x), 1)
#  endif
#endif

enum Op { csum, cmean, cmin, cmax };
typedef hashmap<std::string, int> str2op;

template <typename T>
using annotation_map_t = hashmap<std::string, std::vector<T*>>;
typedef std::vector<char*> strlist;

using chr2bool = hashset<std::string>;

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
//typedef hashmap<std::string, uint64_t> mate2len;
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
//static const int OUT_BUFF_SZ=4000000;
static const int OUT_BUFF_SZ=50;
static const int COORD_STR_LEN=34;

static const void print_version() {
    //fprintf(stderr, "megadepth %s\n", string(MEGADEPTH_VERSION).c_str());
    std::cout << "megadepth " << std::string(MEGADEPTH_VERSION) << std::endl;
}

int my_write(void* fh, char* buf, uint32_t buf_len) {
    return fprintf((FILE*) fh, "%s", buf); 
}

int my_gzwrite(void* fh, char* buf, uint32_t buf_len) {
    return bgzf_write((BGZF*)fh, buf, buf_len);
    //return gzwrite(*((gzFile*) fh), buf, buf_len); 
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

static const long get_longest_target_size(const bam_hdr_t * hdr) {
    long max = 0;
    for(int32_t i = 0; i < hdr->n_targets; i++) {
        if(hdr->target_len[i] > max)
            max = hdr->target_len[i];
    }
    return max;
}

static void reset_array(int32_t* arr, const long arr_sz) {
/*#if __AVX2__
    __m256i zero = _mm256_setzero_si256();
    const size_t nsimd = arr_sz / 4;
    const size_t nsimd4 = (nsimd / 4) * 4;
    size_t i = 0;
    for(; i < nsimd4; i += 4) {
        _mm256_storeu_si256((__m256i *)(arr + 4 * i), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 1)), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 2)), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 3)), zero);
    }
    for(;i < nsimd; ++i) {
        _mm256_storeu_si256((__m256i *)(arr + 4 * i), zero);
    }
    for(i *= 4; i < arr_sz; ++i) {
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
#else*/
    std::memset(arr, 0, sizeof(int32_t) * arr_sz);
//#endif
}

static void reset_array(uint32_t* arr, const long arr_sz) {
/*#if __AVX2__
    __m256i zero = _mm256_setzero_si256();
    const size_t nsimd = arr_sz / 4;
    const size_t nsimd4 = (nsimd / 4) * 4;
    size_t i = 0;
    for(; i < nsimd4; i += 4) {
        _mm256_storeu_si256((__m256i *)(arr + 4 * i), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 1)), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 2)), zero);
        _mm256_storeu_si256((__m256i *)(arr + 4 * (i + 3)), zero);
    }
    for(;i < nsimd; ++i) {
        _mm256_storeu_si256((__m256i *)(arr + 4 * i), zero);
    }
    for(i *= 4; i < arr_sz; ++i) {
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
#else*/
    std::memset(arr, 0, sizeof(uint32_t) * arr_sz);
//#endif
}
                                
typedef hashmap<uint32_t,uint32_t> int2int;
static uint64_t print_array(const char* prefix, 
                        char* chrm,
                        int32_t tid,
                        const int32_t* arr, 
                        const long arr_sz,
                        const bool skip_zeros,
                        bigWigFile_t* bwfp,
                        FILE* cov_fh,
                        //gzFile& gcov_fh,
                        BGZF* gcov_fh,
                        hts_idx_t* cidx,
                        const bool dont_output_coverage = false) {
    bool first = true;
    bool first_print = true;
    uint32_t running_value = 0;
    uint32_t last_pos = 0;
    uint64_t auc = 0;
    //from https://stackoverflow.com/questions/27401388/efficient-gzip-writing-with-gzprintf
    int chrnamelen = strlen(chrm);
    int total_line_len = chrnamelen + COORD_STR_LEN;
    int num_lines_per_buf = round(OUT_BUFF_SZ / total_line_len) - 3;
    //fprintf(stdout, "num_lines_per_buf %d\n", num_lines_per_buf);
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
    uint32_t buf_len = 0;
    int bytes_written = 0;
    char* startp = new char[32];
    char* endp = new char[32];
    char* valuep = new char[32];
    float running_value_ = 0.0;
    int check = 0;
    for(uint32_t i = 0; i < arr_sz; i++) {
        if(first || arr[i] != 0) {
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
                            /*if(buf_written >= num_lines_per_buf) {
                                bufptr[0]='\0';
                                (*printPtr)(cfh, buf, buf_len);
                                bufptr = buf;
                                buf_written = 0;
                                buf_len = 0;
                            }*/
                            std::memcpy(bufptr, chrm, chrnamelen);
                            bufptr += chrnamelen;
                            //start bytes_written from here
                            bytes_written = chrnamelen;
                            
                            *bufptr='\t';
                            bufptr+=1;
                            bytes_written++;
                            //idea from https://github.com/brentp/mosdepth/releases/tag/v0.2.9
                            uint32_t digits = u32toa_countlut(last_pos, bufptr, '\t');
                            bufptr+=digits+1;
                            bytes_written+=digits+1;
                            
                            digits = u32toa_countlut(i, bufptr, '\t');
                            bufptr+=digits+1;
                            bytes_written+=digits+1;
                            
                            digits = u32toa_countlut(running_value, bufptr, '\n');
                            bufptr+=digits+1;
                            bytes_written+=digits+1;
                            bufptr[0]='\0';
                            (*printPtr)(cfh, buf, bytes_written);
                            //TODO add this for last line
                            if(cidx) {
                                if(hts_idx_push(cidx, tid, last_pos, i, bgzf_tell((BGZF*) cfh), 1) < 0)
                                    fprintf(stderr,"error writing line in index at coordinates: %s:%u-%u, exiting\n",chrm,last_pos,i);
                            }
                            
                            buf_len += bytes_written;
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
            running_value += arr[i];
            last_pos = i;
        }
    }
    char last_line[1024];
    if(!first) {
        if(running_value > 0 || !skip_zeros) {
            auc += (arr_sz - last_pos) * ((long) running_value);
            if(not dont_output_coverage) {
                if(bwfp && first_print) {
                    running_value_ = static_cast<float>(running_value);
                    bwAddIntervals(bwfp, &chrm, &last_pos, (uint32_t*) &arr_sz, &running_value_, 1);
                }
                else if(bwfp) {
                    running_value_ = static_cast<float>(running_value);
                    bwAppendIntervals(bwfp, &last_pos, (uint32_t*) &arr_sz, &running_value_, 1);
                }
                else {
                    /*if(buf_written > 0) {
                        bufptr[0]='\0';
                        (*printPtr)(cfh, buf, buf_len);
                    }*/
                    buf_len = sprintf(last_line, "%s\t%u\t%lu\t%u\n", chrm, last_pos, arr_sz, running_value);
                    //buf_len = sprintf(last_line, "%s\t%u\t%lu\t%.0f\n", chrm, last_pos, arr_sz, running_value);
                    (*printPtr)(cfh, last_line, buf_len);
                }
            }
        }
    }
    return auc;
}


typedef hashmap<std::string, uint32_t*> read2len;
static const int32_t calculate_coverage(const bam1_t *rec, int32_t* coverages, 
                                        int32_t* unique_coverages, const bool double_count, 
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
    const std::string tn(qname);
    int32_t end_pos = bam_endpos(rec);
    uint32_t mate_passes_quality = 0;
    //-----First Mate Check
    //if we're the first mate and
    //we're avoiding double counting and we're a proper pair
    //and we overlap with our mate, then store our cigar + length
    //for the later mate to adjust its coverage appropriately
    if(coverages && !double_count && (rec->core.flag & BAM_FPROPER_PAIR) == 2) {
        auto mit = overlapping_mates->find(tn);
    
        if(rec->core.tid == rec->core.mtid &&
                end_pos > mrefpos && 
                refpos <= mrefpos &&
                mit == overlapping_mates->end()) {
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
        else if(mit != overlapping_mates->end()) {
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
            overlapping_mates->erase(mit);
            n_mspans = mspans_idx;
            mendpos = malgn_end_pos;
        }
    }
    mspans_idx = 0;
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
                coverages[algn_end_pos]++;
                coverages[algn_end_pos+len]--;
                if(unique && passing_qual) {
                    unique_coverages[algn_end_pos]++;
                    unique_coverages[algn_end_pos+len]--;
                }
                //TODO: skip updating the coverage at every value
                /*for(z = algn_end_pos; z < algn_end_pos + len; z++) {
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
                }*/    
            }
            algn_end_pos += len;
        }
    }
    if(mspans) {
        for(k = 0; k < n_mspans; k++)
            delete[] mspans[k];
        delete[] mspans;
    }
    return algn_end_pos;
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



template <typename T>
int go_bam(const char* bam_arg, int argc, const char** argv, Op op, htsFile *bam_fh, int nthreads, bool keep_order, bool has_annotation, FILE* afp, BGZF* afpz, annotation_map_t<T>* annotations, chr2bool* annotation_chrs_seen, const char* prefix, bool sum_annotation, strlist* chrm_order, FILE* auc_file, uint64_t num_annotations) {
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
    hts_set_threads(bam_fh, nthreads);
    /*hFile* hfs = hts_hfile(bam_fh);
    hfile_set_blksize(hfs, 1024*1024*3);*/
    //runs slower
    //hts_set_opt(bam_fh, HTS_OPT_BLOCK_SIZE, 1024*1024*100);
    
    //setup list of callbacks for the process_cigar()
    //this is so we only have to walk the cigar for each alignment ~1 time
    size_t recs = 0;
    const bool double_count = has_option(argv, argv+argc, "--double-count");
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
    int32_t* coverages = nullptr;
    bool compute_coverage = false;
    read2len overlapping_mates;
    //--coverage -> output perbase coverage to STDOUT (compute_coverage=true)
    bool coverage_opt = has_option(argv, argv+argc, "--coverage");
    bool bigwig_opt = has_option(argv, argv+argc, "--bigwig");
    bool auc_opt = has_option(argv, argv+argc, "--auc") || argc == 1;
    bool dont_output_coverage = !(coverage_opt || bigwig_opt);
    FILE* cov_fh = stdout;
    bool gzip = has_option(argv, argv+argc, "--gzip");
    bool no_coverage_stdout = gzip || has_option(argv, argv+argc, "--no-coverage-stdout");
    //gzFile gcov_fh;
    BGZF* gcov_fh = nullptr;
    hts_idx_t* cidx = nullptr;
    
    if(coverage_opt) {
        compute_coverage = true;
        chr_size = get_longest_target_size(hdr);
        coverages = new int32_t[chr_size];
        if(coverage_opt && no_coverage_stdout) {
            char cov_fn[1024];
            if(gzip) {
                sprintf(cov_fn, "%s.coverage.tsv.gz", prefix);
                //gzFile gcov_fh_ = gzopen(cov_fn,"w");
                //gcov_fh = gzopen(cov_fn,"w1");
                gcov_fh = bgzf_open(cov_fn,"w10");
                cov_fh = nullptr;
                //from https://github.com/samtools/htslib/blob/c9175183c42382f1030503e88ca7e60cb9c08536/sam.c#L923
                //and https://github.com/brentp/hts-nim/blob/0eaa867e747d3bc844b5ecb575796e4688b966f5/src/hts/csi.nim#L34
                int min_shift = 14;
                int n_lvls = (TBX_MAX_SHIFT - min_shift + 2) / 3;
                int fmt = HTS_FMT_CSI;
                //n_lvls = adjust_n_lvls(min_shift, n_lvls, chr_size);
                cidx = hts_idx_init(0, fmt, 0, min_shift, n_lvls);
            }
            else {
                sprintf(cov_fn, "%s.coverage.tsv", prefix);
                cov_fh = fopen(cov_fn,"w");
            }
        }
    }
    mate2len* frag_mates = new mate2len(1);
    char cov_prefix[50]="";
    int32_t ptid = -1;
    uint32_t len = 0;
    int filter_in_mask = 0xFFFFFFFF;
    uint64_t reads_processed = 0;
    //filter out alignments with either BAM_FUNMAP and/or BAM_FSECONDARY flags set by default (260)
    int filter_out_mask = 260;
    bam1_t* rec_ = bam_init1();
    uint64_t num_annotations_ = 0;
    if(dont_output_coverage && !auc_opt)
        num_annotations_ = num_annotations;
    BAMIterator<T> bitr(rec_, bam_fh, hdr, bam_arg, annotations, num_annotations, chrm_order);
    BAMIterator<T> end(nullptr, nullptr, nullptr);
    //while(sam_read1(bam_fh, hdr, rec) >= 0) {
    for(++bitr; bitr != end; ++bitr) {
        recs++;
        rec = *bitr;
        bam1_core_t *c = &rec->core;
        //read name
        char* qname = bam_get_qname(rec);
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

            //*******Reference coverage tracking
            if(compute_coverage) {
                if(tid != ptid) {
                    if(ptid != -1) {
                        overlapping_mates.clear();
                        chr_size = hdr->target_len[ptid];
                        //sprintf(cov_prefix, "cov\t%d", ptid);
                        if(coverage_opt)
                            all_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, coverages, chr_size, false, nullptr, cov_fh, gcov_fh, cidx, dont_output_coverage);
                    }
                    reset_array(coverages, chr_size);
                }
                end_refpos = calculate_coverage(rec, coverages, nullptr, double_count, 0, &overlapping_mates, &total_intron_len);
            }
            ptid = tid;
        }
    }
    if(compute_coverage) {
        if(ptid != -1) {
            chr_size = hdr->target_len[ptid];
            //sprintf(cov_prefix, "cov\t%d", ptid);
            if(coverage_opt || bigwig_opt || auc_opt)
                all_auc += print_array(cov_prefix, hdr->target_name[ptid], ptid, coverages, chr_size, false, nullptr, cov_fh, gcov_fh, cidx, dont_output_coverage);
        }
        delete[] coverages;
    }
    //for writing out an index for BGZipped coverage BED files
    //mostly from Tabix source: https://github.com/samtools/htslib/blob/develop/tabix.c
    char temp_afn[1024];
    tbx_conf_t tconf = tbx_conf_bed;
    int min_shift = 14; 
    if(cov_fh && cov_fh != stdout)
        fclose(cov_fh);
    if(gzip && gcov_fh) {
        sprintf(temp_afn, "%s.coverage.tsv.gz", prefix);
        char temp_afni[1024];
        sprintf(temp_afni, "%s.coverage.tsv.gz.csi", prefix);
        if(hts_idx_finish(cidx, bgzf_tell(gcov_fh)) != 0)
            fprintf(stderr,"Error finishing BGZF index for base coverage, skipping\n");
        else {
            /*tbx_t *tbx;
            tbx = (tbx_t*)calloc(1, sizeof(tbx_t));
            tbx->conf = tconf; 
            tbx->idx = cidx;
            tbx_set_meta(tbx);*/
            if(hts_idx_save_as(cidx, temp_afn, temp_afni, HTS_FMT_CSI) != 0)
                fprintf(stderr,"Error saving BGZF index for base coverage, skipping\n");
        }
        bgzf_close(gcov_fh);
    }

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
    int err = 0;
    bool gzip = has_option(argv, argv+argc, "--gzip");
    const char* prefix = fname_arg;
    uint64_t num_annotations = 0;
    if(has_option(argv, argv+argc, "--prefix"))
            prefix = *(get_option(argv, argv+argc, "--prefix"));
    FILE* afp = nullptr;
    BGZF* afpz = nullptr;
    //if no args are passed in other than a file (BAM or BW)
    //then just compute the auc 
    FILE* auc_file = nullptr;
    bool keep_order = true;
    bool has_annotation = false;
    bool sum_annotation = false;
    annotation_map_t<T>* annotations;
    chr2bool* annotation_chrs_seen;
    return go_bam(fname_arg, argc, argv, op, bam_fh, nthreads, keep_order, has_annotation, afp, afpz, annotations, annotation_chrs_seen, prefix, sum_annotation, nullptr, nullptr, num_annotations);
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
    if(has_option(argv, argv + argc, "--version")) {
        print_version();
        return 0;
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
