#!/usr/bin/env bash

rm -f test.bam.auc.tsv test.bam.all.tsv test.bam.unique.tsv test.bam.frags.tsv test.bam.alts.tsv test.bam.softclip.tsv

time ./bcd_runner tests/test.bam --threads 4 --coverage --no-head --bigwig test.bam --auc test.bam --min-unique-qual 10 --annotation tests/test_exons.bed test.bam --frag-dist test.bam --alts test.bam --include-softclip test.bam --only-polya --read-ends test --test-polya > test_run_out 2>&1

diff <(sort tests/test.bam.orig.frags.tsv) <(sort test.bam.frags.tsv)
diff tests/test.bam.orig.alts.tsv test.bam.alts.tsv
diff tests/test.bam.orig.softclip.tsv test.bam.softclip.tsv
for f in all unique; do
    diff tests/test.bam.mosdepth.${f}.per-base.exon_sums.tsv test.bam.${f}.tsv
done
diff tests/test.bam.mosdepth.bwtool.all_aucs test.bam.auc.tsv

cat test.starts.tsv test.ends.tsv | sort -k1,1 -k2,2n -k3,3n > test_starts_ends.tsv
diff test_starts_ends.tsv <(sort -k1,1 -k2,2n -k3,3n tests/test.bam.read_ends.both.unique.tsv)

time ./bcd_runner tests/test2.bam --threads 4 --no-head --junctions test2.bam > test2_run_out 2>&1

diff tests/test2.bam.jxs.tsv test2.bam.jxs.tsv

#test just total auc
time ./bcd_runner test.bam.all.bw --auc  2>&1 > test_run_out | grep "AUC" > test.bw1.total_auc
diff test.bw1.total_auc tests/testbw1.total_auc

#test bigwig2sums/auc
time ./bcd_runner test.bam.all.bw --annotation tests/testbw1.bed test.bam.bw1 | fgrep AUC > test.bw1.annot_auc
diff test.bam.bw1.all.tsv tests/testbw1.bed.out.tsv
diff test.bw1.annot_auc tests/testbw1.annot_auc

##use different order in BED file from what's in BW to test keep_order == true
time ./bcd_runner test.bam.all.bw --annotation tests/testbw2.bed test.bam.bw2 | fgrep AUC > test.bw2.annot_auc
diff test.bam.bw2.all.tsv tests/testbw2.bed.out.tsv
diff test.bw2.annot_auc tests/testbw2.annot_auc

#test bigwig2mean
time ./bcd_runner test.bam.all.bw --op mean --annotation tests/testbw2.bed bw2.mean >> test_run_out 2>&1
diff bw2.mean.all.tsv tests/testbw2.bed.mean

#test bigwig2min
time ./bcd_runner test.bam.all.bw --op min --annotation tests/testbw2.bed bw2.min >> test_run_out 2>&1
diff bw2.min.all.tsv tests/testbw2.bed.min

#test bigwig2max
time ./bcd_runner test.bam.all.bw --op max --annotation tests/testbw2.bed bw2.max >> test_run_out 2>&1
diff bw2.max.all.tsv tests/testbw2.bed.max

#now test same-start alignments for overlapping pairs
./bcd_runner tests/test3.bam --coverage --no-head > t3.tsv
diff tests/test3.out.tsv t3.tsv

#with uniques
./bcd_runner tests/test3.bam --coverage --no-head --min-unique-qual 10 --bigwig test3 --auc test3
diff tests/test3.auc.out.tsv test3.auc.tsv

#long reads support for junctions
./bcd_runner tests/long_reads.bam --junctions long_reads.bam --long-reads
diff tests/long_reads.bam.jxs.tsv long_reads.bam.jxs.tsv
