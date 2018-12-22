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

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/bam_io.h>

#include "docopt/docopt.h"

static const char USAGE[] =
R"(BAM and BigWig utility.

Usage:
  bwam nonref <bam> [options]
  bwam bigwig <bam> [options]
  bwam both <bam> [options]

Options:
  -h --help           Show this screen.
  --version           Show version.
  --include-softclip  Print a record for soft-clipped bases.
  --include-n         Print mismatch records when mismatched read base is N
  --echo-sam          Print a SAM record for each aligned read
  --double-count      Allow overlapping ends of PE read to count twice toward
                      coverage
  --require-mdz       Quit with error unless MD:Z field exists everywhere it's
                      expected
  --print-qual        Print quality values for mismatched bases
)";

using namespace seqan;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::tuple;
using std::get;
using std::make_tuple;
using std::stringstream;


/**
 * Parse given MD:Z extra field into a vector of MD:Z operations.
 */
static void parse_mdz(
        const CharString& mdz,
        vector<tuple<char, int, CharString>>& ops)
{
    int i = 0;
    const size_t ln = length(mdz);
    bool saw_d = false;
    while(i < ln) {
        if(isdigit(mdz[i])) {
            int run = 0;
            while(i < ln && isdigit(mdz[i])) {
                run *= 10;
                run += (int)(mdz[i] - '0');
                i++;
            }
            if(run > 0) {
                ops.emplace_back(make_tuple('=', run, CharString("")));
            }
        } else if(isalpha(mdz[i])) {
            size_t st = i;
            while(i < ln && isalpha(mdz[i])) i++;
            assert(i > st);
            ops.emplace_back(make_tuple('X', i - st, infix(mdz, st, i)));
        } else if(mdz[i] == '^') {
            saw_d = true;
            i++;
            size_t st = i;
            while (i < ln && isalpha(mdz[i])) i++;
            assert(i > st);
            ops.emplace_back(make_tuple('^', i - st, infix(mdz, st, i)));
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
        const BamAlignmentRecord& rec,
        vector<tuple<char, int, CharString>>& mdz,
        bool print_qual = false,
        bool include_ss = false,
        bool include_n_mms = false)
{
    const IupacString& seq{rec.seq};
    assert(rec.beginPos >= 0);
    size_t mdzi = 0, seq_off = 0, ref_off = rec.beginPos;
    for (CigarElement<> e : rec.cigar) {
        const char &op = e.operation;
        if((strchr("DNMX=", op) != nullptr) && mdzi >= mdz.size()) {
            stringstream ss;
            ss << "Found read-consuming CIGAR op after MD:Z had been exhausted" << endl;
            throw std::runtime_error(ss.str());
        }
        if (op == 'M' || op == 'X' || op == '=') {
            // Look for block matches and mismatches in MD:Z string
            size_t runleft = e.count;
            while (runleft > 0 && mdzi < mdz.size()) {
                char mdz_op;
                size_t mdz_run;
                CharString mdz_str;
                std::tie(mdz_op, mdz_run, mdz_str) = mdz[mdzi];
                size_t run_comb = std::min(runleft, mdz_run);
                runleft -= run_comb;
                assert(mdz_op == 'X' or mdz_op == '=');
                if (mdz_op == '=') {
                    // nop
                } else {
                    assert(mdz_op == 'X');
                    assert(length(mdz_str) == run_comb);
                    if(!include_n_mms && run_comb == 1 && rec.seq[seq_off] == 'N') {
                        // skip
                    } else {
                        cout << rec.rID << ',' << ref_off << ",X,"
                             << infix(rec.seq, seq_off, seq_off + run_comb);
                        if(print_qual) {
                            cout << ',' << infix(rec.qual, seq_off, seq_off + run_comb);
                        }
                        cout << endl;
                    }
                }
                seq_off += run_comb;
                ref_off += run_comb;
                if (run_comb < mdz_run) {
                    assert(mdz_op == '=');
                    get<1>(mdz[mdzi]) -= run_comb;
                } else {
                    mdzi++;
                }
            }
        } else if(op == 'I') {
            cout << rec.rID << ',' << ref_off << ",I,"
                 << infix(rec.seq, seq_off, seq_off + e.count) << endl;
            seq_off += e.count;
        } else if(op == 'S') {
            if(include_ss) {
                cout << rec.rID << ',' << ref_off << ",S,"
                     << infix(rec.seq, seq_off, seq_off + e.count) << endl;
                seq_off += e.count;
            }
        } else if (op == 'D') {
            char mdz_op;
            size_t mdz_run;
            CharString mdz_str;
            std::tie(mdz_op, mdz_run, mdz_str) = mdz[mdzi];
            assert(mdz_op == '^');
            assert(e.count == mdz_run);
            assert(length(mdz_str) == e.count);
            mdzi++;
            cout << rec.rID << ',' << ref_off << ",D," << e.count << endl;
            ref_off += e.count;
        } else if (op == 'N') {
            ref_off += e.count;
        } else if (op == 'H') {
        } else if (op == 'P') {
        } else {
            stringstream ss;
            ss << "No such CIGAR operation as \"" << op << "\"";
            throw std::runtime_error(ss.str());
        }
    }
    assert(mdzi == mdz.size());
}

static void output_from_cigar(const BamAlignmentRecord& rec) {
    bool simple_cigar = length(rec.cigar) == 1;
    assert(!simple_cigar || rec.cigar[0].operation == 'M');
    if(simple_cigar) {
        return;
    }
    int32_t refpos = rec.beginPos;
    int32_t seqpos = 0;
    for(CigarElement<> e : rec.cigar) {
        switch(e.operation) {
            case 'D': {
                cout << rec.rID << ',' << refpos << ",D," << e.count << endl;
                refpos += e.count;
                break;
            }
            case 'S':
            case 'I': {
                cout << rec.rID << ',' << refpos << ',' << e.operation
                     << ',' << infix(rec.seq, seqpos, seqpos+e.count) << endl;
                seqpos += e.count;
                break;
            }
            case 'N': {
                refpos += e.count;
                break;
            }
            case 'M':
            case 'X':
            case '=': {
                seqpos += e.count;
                refpos += e.count;
                break;
            }
            case 'H':
            case 'P': { break; }
            default: {
                stringstream ss;
                ss << "No such CIGAR operation as \"" << e.operation << "\"";
                throw std::runtime_error(ss.str());
            }
        }
    }
}

int main(int argc, const char** argv) {
    std::map<std::string, docopt::value> args
            = docopt::docopt(USAGE,
                             { argv + 1, argv + argc },
                             true,                      // show help if requested
                             "BWAM 0.1");               // version string
    
    if(args["nonref"].asBool()) {
        assert(args.find("<bam>") != args.end());
        std::string bam_arg{args["<bam>"].asString()};
        // Open input file, BamFileIn can read SAM and BAM files.
        BamFileIn bamFileIn;
        if(!open(bamFileIn, bam_arg.c_str())) {
            std::cerr << "ERROR: Could not open \"" << bam_arg << "\"" << std::endl;
            return 1;
        }
    
        BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());
    
        const bool print_qual = args["--print-qual"].asBool();
        const bool include_ss = args["--include-softclip"].asBool();
        const bool include_n_mms = args["--include-n"].asBool();
        size_t recs = 0;
        try {
            BamHeader header;
            BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());
            readHeader(header, bamFileIn);
            BamAlignmentRecord record;
            vector<tuple<char, int, CharString>> mdzbuf;
            CharString mdz;
            while(!atEnd(bamFileIn)) {
                readRecord(record, bamFileIn);
                recs++;
                bool had_d = false;
                if(!hasFlagUnmapped(record)) {
                    if(args["--echo-sam"].asBool()) {
                        writeRecord(bamFileOut, record);
                    }
                    unsigned tagIdx = 0;
                    BamTagsDict tagsDict(record.tags);
                    if(!findTagKey(tagIdx, tagsDict, "MD")) {
                        if(args["--require-mdz"].asBool()) {
                            stringstream ss;
                            ss << "No MD:Z extra field for aligned read \"" << record.qName << "\"";
                            throw std::runtime_error(ss.str());
                        }
                        output_from_cigar(record); // just use CIGAR
                    } else {
                        clear(mdz);
                        extractTagValue(mdz, tagsDict, tagIdx);
                        mdzbuf.clear();
                        parse_mdz(mdz, mdzbuf);
                        output_from_cigar_mdz(
                                record, mdzbuf, print_qual,
                                include_ss, include_n_mms); // use CIGAR and MD:Z
                    }
                }
            }
        }
        catch (Exception const & e) {
            cout << "ERROR: " << e.what() << endl;
            return 1;
        }
    
        cout << "Read " << recs << " records" << endl;
    } else if(args["bigwig"].asBool()) {
        cout << "BigWig mode not implemented yet!" << endl;
    }
    
    return 0;
}
