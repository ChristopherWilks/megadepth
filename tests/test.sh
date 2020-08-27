#!/usr/bin/env bash
set -xe

static=$1

rm -f test.bam.auc.tsv test.bam.annotation.tsv test.bam.unique.tsv test.bam.frags.tsv test.bam.alts.tsv test.bam.softclip.tsv

if [[ -z $static ]]; then
    time ./md_runner http://stingray.cs.jhu.edu/data/temp/test.bam --prefix test.bam --threads 4 --no-head --bigwig --auc --min-unique-qual 10 --annotation tests/test_exons.bed --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout > test_run_out 2>&1
else
    time ./md_runner tests/test.bam --prefix test.bam --threads 4 --no-head --bigwig --auc --min-unique-qual 10 --annotation tests/test_exons.bed --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout > test_run_out 2>&1
fi

diff <(sort tests/test.bam.orig.frags.tsv) <(sort test.bam.frags.tsv)
diff tests/test.bam.orig.alts.tsv test.bam.alts.tsv
diff tests/test.bam.orig.softclip.tsv test.bam.softclip.tsv
for f in annotation unique; do
    diff tests/test.bam.mosdepth.${f}.per-base.exon_sums.tsv test.bam.${f}.tsv
done
diff tests/test.bam.mosdepth.bwtool.all_aucs test.bam.auc.tsv

cat test.bam.starts.tsv test.bam.ends.tsv | sort -k1,1 -k2,2n -k3,3n > test_starts_ends.tsv
diff test_starts_ends.tsv <(sort -k1,1 -k2,2n -k3,3n tests/test.bam.read_ends.both.unique.tsv)

time ./md_runner tests/test2.bam --threads 4 --no-head --junctions --prefix test2.bam > test2_run_out 2>&1

diff tests/test2.bam.jxs.tsv test2.bam.jxs.tsv

#test just total auc
time ./md_runner test.bam.all.bw | grep "AUC" > test.bw1.total_auc
diff test.bw1.total_auc tests/testbw1.total_auc

#test bigwig2sums/auc
time ./md_runner test.bam.all.bw --annotation tests/testbw1.bed --prefix test.bam.bw1 --no-annotation-stdout | fgrep AUC > test.bw1.annot_auc
diff test.bam.bw1.annotation.tsv tests/testbw1.bed.out.tsv
diff test.bw1.annot_auc tests/testbw1.annot_auc

##use different order in BED file from what's in BW to test keep_order == true
time ./md_runner test.bam.all.bw --annotation tests/testbw2.bed --prefix test.bam.bw2 --no-annotation-stdout | fgrep AUC > test.bw2.annot_auc
diff test.bam.bw2.annotation.tsv tests/testbw2.bed.out.tsv
diff test.bw2.annot_auc tests/testbw2.annot_auc

#test bigwig2mean
time ./md_runner test.bam.all.bw --op mean --annotation tests/testbw2.bed --prefix bw2.mean --no-annotation-stdout >> test_run_out 2>&1
diff bw2.mean.annotation.tsv tests/testbw2.bed.mean

#test bigwig2min
time ./md_runner test.bam.all.bw --op min --annotation tests/testbw2.bed --prefix bw2.min --no-annotation-stdout >> test_run_out 2>&1
diff bw2.min.annotation.tsv tests/testbw2.bed.min

#test bigwig2max
time ./md_runner test.bam.all.bw --op max --annotation tests/testbw2.bed --prefix bw2.max --no-annotation-stdout >> test_run_out 2>&1
diff bw2.max.annotation.tsv tests/testbw2.bed.max

#now test same-start alignments for overlapping pairs
./md_runner tests/test3.bam --coverage --no-head > t3.tsv
diff tests/test3.out.tsv t3.tsv

#with uniques
./md_runner tests/test3.bam --coverage --no-head --min-unique-qual 10 --bigwig --auc --prefix test3
diff tests/test3.auc.out.tsv test3.auc.tsv

#long reads support for junctions
./md_runner tests/long_reads.bam --junctions --prefix long_reads.bam --long-reads
diff tests/long_reads.bam.jxs.tsv long_reads.bam.jxs.tsv

#test bigwig2sum on remote BW
if [[ -z $static ]]; then
    time ./md_runner http://stingray.cs.jhu.edu/data/temp/megadepth.test.bam.all.bw --op mean --annotation tests/testbw2.bed --prefix bw2.remote.mean --no-annotation-stdout >> test_run_out 2>&1
    diff bw2.remote.mean.annotation.tsv tests/testbw2.bed.mean
fi

##only print sums use different order in BED file from what's in BW to test keep_order == true
time ./md_runner test.bam.all.bw --sums-only --annotation tests/testbw2.bed --auc --prefix test.bam.bw2 > test.bam.bw2.annotation.tsv
diff test.bam.bw2.annotation.tsv <(cut -f 4 tests/testbw2.bed.out.tsv)

#clean up
rm -f test*tsv test*auc bw2* test3* test2* t3.tsv long_reads.bam.jxs.tsv test_run_out *null*.unique.tsv test.*.bw
