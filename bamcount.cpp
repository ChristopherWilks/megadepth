//
// Created by Ben Langmead on 2018-12-12.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>
#include <string>
#include <cerrno>
#include <htslib/sam.h>
#include <htslib/bgzf.h>

#include "docopt/docopt.h"

static const char USAGE[] =
R"(BAM and BigWig utility.

Usage:
  bamcount nonref <bam> [options]
  bamcount bigwig <bam> [options]
  bamcount both <bam> [options]

Options:
  -h --help           Show this screen.
  --version           Show version.
  --include-softclip  Print a record for soft-clipped bases.
  --include-n         Print mismatch records when mismatched read base is N
  --no-head           Don't print sequence names and lengths in header
  --delta             Print POS field as +/- delta from previous
  --echo-sam          Print a SAM record for each aligned read
  --double-count      Allow overlapping ends of PE read to count twice toward
                      coverage
  --require-mdz       Quit with error unless MD:Z field exists everywhere it's
                      expected
  --print-qual        Print quality values for mismatched bases
)";

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::tuple;
using std::get;
using std::make_tuple;
using std::stringstream;

/**
 * Holds an MDZ "operation"
 * op can be 
 */
struct MdzOp {
    char op;
    int run;
    char str[1024];
};

static inline std::ostream& seq_substring(std::ostream& os, const uint8_t *str, size_t off, size_t run) {
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
            memcpy(ops.back().str, mdz + st, (size_t)(i - st));
            ops.back().str[i - st] = '\0';
        } else if(mdz[i] == '^') {
            i++;
            int st = i;
            while (i < mdz_len && isalpha(mdz[i])) i++;
            assert(i > st);
            ops.emplace_back(MdzOp{'^', i - st, ""});
            memcpy(ops.back().str, mdz + st, (size_t)(i - st));
            ops.back().str[i - st] = '\0';
        } else {
            stringstream ss;
            ss << "Unknown MD:Z operation: \"" << mdz[i] << "\"";
            throw std::runtime_error(ss.str());
        }
    }
}

#if 0
/**
 * Prints a stacked version of an alignment
 */
static void cigar_mdz_to_stacked(
        const BamAlignmentRecord& rec,
        vector<tuple<char, int, CharString>>& mdz,
        IupacString& rds,
        IupacString& rfs)
{
    const IupacString& seq{rec.seq};
    size_t mdz_off = 0, seq_off = 0;
    for (CigarElement<> e : rec.cigar) {
        const char &op = e.operation;
        bool ref_consuming = (strchr("DNMX=", op) != nullptr);
        if(ref_consuming && mdz_off >= mdz.size()) {
            stringstream ss;
            ss << "Found read-consuming CIGAR op after MD:Z had been exhausted" << endl;
            throw std::runtime_error(ss.str());
        }
        if (op == 'M' || op == 'X' || op == '=') {
            // Look for block matches and mismatches in MD:Z string
            size_t runleft = e.count;
            while (runleft > 0 && mdz_off < mdz.size()) {
                char mdz_op;
                size_t mdz_run;
                CharString mdz_str;
                std::tie(mdz_op, mdz_run, mdz_str) = mdz[mdz_off];
                size_t run_comb = std::min(runleft, mdz_run);
                runleft -= run_comb;
                assert(mdz_op == 'X' or mdz_op == '=');
                append(rds, infix(seq, seq_off, seq_off + run_comb));
                if (mdz_op == '=') {
                    append(rfs, infix(seq, seq_off, seq_off + run_comb));
                } else {
                    assert(length(mdz_str) == run_comb);
                    append(rfs, mdz_str);
                }
                seq_off += run_comb;
                if (run_comb < mdz_run) {
                    assert(mdz_op == '=');
                    get<1>(mdz[mdz_off]) -= run_comb;
                } else {
                    mdz_off++;
                }
            }
        } else if (op == 'I') {
            append(rds, infix(seq, seq_off, seq_off + e.count));
            for (size_t i = 0; i < e.count; i++) {
                append(rfs, '-');
            }
        } else if (op == 'D') {
            char mdz_op;
            size_t mdz_run;
            CharString mdz_str;
            std::tie(mdz_op, mdz_run, mdz_str) = mdz[mdz_off];
            assert(mdz_op == '^');
            assert(e.count == mdz_run);
            assert(length(mdz_str) == e.count);
            mdz_off++;
            for (size_t i = 0; i < e.count; i++) {
                append(rds, '-');
            }
            append(rfs, mdz_run);
        } else if (op == 'N') {
            for (size_t i = 0; i < e.count; i++) {
                append(rds, '-');
                append(rfs, '-');
            }
        } else if (op == 'S') {
            append(rds, infix(seq, seq_off, seq_off + e.count));
            for (size_t i = 0; i < e.count; i++) {
                append(rfs, '-');
            }
            seq_off += e.count;
        } else if (op == 'H') {
        } else if (op == 'P') {
        } else {
            stringstream ss;
            ss << "No such CIGAR operation as \"" << op << "\"";
            throw std::runtime_error(ss.str());
        }
    }
    assert(mdz_off == mdz.size());
}
#endif

static void output_from_cigar_mdz(
        const bam1_t *rec,
        std::vector<MdzOp>& mdz,
        bool print_qual = false,
        bool include_ss = false,
        bool include_n_mms = false,
        bool delta = false)
{
    uint8_t *seq = bam_get_seq(rec);
    uint8_t *qual = bam_get_qual(rec);
    uint32_t *cigar = bam_get_cigar(rec);
    size_t mdzi = 0, seq_off = 0;
    int32_t ref_off = rec->core.pos;
    for(int k = 0; k < rec->core.n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        int run = bam_cigar_oplen(cigar[k]);
        if((strchr("DNMX=", BAM_CIGAR_STR[op]) != nullptr) && mdzi >= mdz.size()) {
            stringstream ss;
            ss << "Found read-consuming CIGAR op after MD:Z had been exhausted" << endl;
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
                    uint8_t cread = bam_seqi(seq, seq_off);
                    if(!include_n_mms && run_comb == 1 && cread == 'N') {
                        // skip
                    } else {
                        cout << rec->core.tid << ',' << ref_off << ",X,";
                        seq_substring(cout, seq, seq_off, run_comb);
                        if(print_qual) {
                            cout << ',';
                            cstr_substring(cout, qual, seq_off, run_comb);
                        }
                        cout << endl;
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
            cout << rec->core.tid << ',' << ref_off << ",I,";
            seq_substring(cout, seq, seq_off, run) << endl;
            seq_off += run;
        } else if(op == BAM_CSOFT_CLIP) {
            if(include_ss) {
                cout << rec->core.tid << ',' << ref_off << ",S,";
                seq_substring(cout, seq, seq_off, run) << endl;
                seq_off += run;
            }
        } else if (op == BAM_CDEL) {
            assert(mdz[mdzi].op == '^');
            assert(run == mdz[mdzi].run);
            assert(strlen(mdz[mdzi].str) == run);
            mdzi++;
            cout << rec->core.tid << ',' << ref_off << ",D," << run << endl;
            ref_off += run;
        } else if (op == BAM_CREF_SKIP) {
            ref_off += run;
        } else if (op == BAM_CHARD_CLIP) {
        } else if (op == BAM_CPAD) {
        } else {
            stringstream ss;
            ss << "No such CIGAR operation as \"" << op << "\"";
            throw std::runtime_error(ss.str());
        }
    }
    assert(mdzi == mdz.size());
}

static void output_from_cigar(const bam1_t *rec) {
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
                cout << rec->core.tid << ',' << refpos << ",D," << run << endl;
                refpos += run;
                break;
            }
            case BAM_CSOFT_CLIP:
            case BAM_CINS: {
                cout << rec->core.tid << ',' << refpos << ',' << BAM_CIGAR_STR[op] << ',';
                seq_substring(cout, seq, seqpos, run) << endl;
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
                stringstream ss;
                ss << "No such CIGAR operation as \"" << op << "\"";
                throw std::runtime_error(ss.str());
            }
        }
    }
}

static int get_header(const std::string& bam_fn, bam_hdr_t *& hdr) {
    BGZF *bam_fh = bgzf_open(bam_fn.c_str(), "r");
    if(!bam_fh) {
        std::cerr << "Couldn't open " << bam_fn << ": " << strerror(errno) << endl;
        return 1;
    }
    hdr = bam_hdr_read(bam_fh);
    if(!hdr) {
        std::cerr << "Couldn't read header from " << bam_fn << ": " << strerror(errno) << endl;
        return 1;
    }
    return 0;
}

static void print_header(const bam_hdr_t * hdr) {
    for(int32_t i = 0; i < hdr->n_targets; i++) {
        cout << '@' << i << ','
             << hdr->target_name[i] << ','
             << hdr->target_len[i] << endl;
    }
}

int main(int argc, const char** argv) {
    std::map<std::string, docopt::value> args
            = docopt::docopt(USAGE,
                             { argv + 1, argv + argc },
                             true,                      // show help if requested
                             "Bamcount 0.1");               // version string
    
    if(args["nonref"].asBool()) {
        assert(args.find("<bam>") != args.end());
        std::string bam_arg{args["<bam>"].asString()};
    
        htsFile *bam_fh = sam_open(bam_arg.c_str(), "r");
        if(!bam_fh) {
            std::cerr << "Couldn't open " << bam_arg << ": " << strerror(errno) << endl;
            return 1;
        }

        const bool print_qual = args["--print-qual"].asBool();
        const bool include_ss = args["--include-softclip"].asBool();
        const bool include_n_mms = args["--include-n"].asBool();

        size_t recs = 0;
        bam_hdr_t *hdr = sam_hdr_read(bam_fh);
        if(!hdr) {
            std::cerr << "Couldn't read header for " << bam_arg << std::endl;
            return 1;
        }
        if(!args["--no-head"].asBool()) {
            print_header(hdr);
        }
        std::vector<MdzOp> mdzbuf;
        bam1_t *rec = bam_init1();
        if (!rec) {
            perror(nullptr);
            return 1;
        }
        kstring_t sambuf{ 0, 0, nullptr };

        while(sam_read1(bam_fh, hdr, rec) >= 0) {
            recs++;
            bam1_core_t *c = &rec->core;
            if((c->flag & BAM_FUNMAP) == 0) {
                if(args["--echo-sam"].asBool()) {
                    int ret = sam_format1(hdr, rec, &sambuf);
                    if(ret < 0) {
                        std::cerr << "Could not format SAM record: " << strerror(errno) << std::endl;
                        return 1;
                    }
                    kstring_out(std::cout, &sambuf);
                    std::cout << std::endl;
                }
                const uint8_t *mdz = bam_aux_get(rec, "MD");
                if(!mdz) {
                    if(args["--require-mdz"].asBool()) {
                        stringstream ss;
                        ss << "No MD:Z extra field for aligned read \"" << hdr->target_name[c->tid] << "\"";
                        throw std::runtime_error(ss.str());
                    }
                    output_from_cigar(rec); // just use CIGAR
                } else {
                    mdzbuf.clear();
                    parse_mdz(mdz + 1, mdzbuf); // skip type character at beginning
                    output_from_cigar_mdz(
                            rec, mdzbuf, print_qual,
                            include_ss, include_n_mms); // use CIGAR and MD:Z
                }
            }
        }
    
        cout << "Read " << recs << " records" << endl;
    } else if(args["bigwig"].asBool()) {
        cout << "BigWig mode not implemented yet!" << endl;
    }
    
    return 0;
}
