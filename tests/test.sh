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
diff test_starts_ends.tsv tests/test.bam.read_ends.both.unique.tsv

time ./bcd_runner tests/test2.bam --threads 4 --no-head --junctions test2.bam > test2_run_out 2>&1

diff tests/test2.bam.jxs.tsv test2.bam.jxs.tsv

#test bigwig2sums/auc
time ./bcd_runner test.bam.all.bw --annotation tests/testbw1.bed test.bam.bw1 >> test_run_out 2>&1
diff test.bam.bw1.all.tsv tests/testbw1.bed.out.tsv
tail -n1 test_run_out > test.bam.bw1.auc.tsv 
diff test.bam.bw1.auc.tsv tests/testbw1.bed.auc

##use different order in BED file from what's in BW to test keep_order == true
time ./bcd_runner test.bam.all.bw --annotation tests/testbw2.bed test.bam.bw2 >> test_run_out 2>&1
diff test.bam.bw2.all.tsv tests/testbw2.bed.out.tsv
tail -n1 test_run_out > test.bam.bw2.auc.tsv 
diff test.bam.bw2.auc.tsv tests/testbw2.bed.auc

#test bigwig2mean
time ./bcd_runner test.bam.all.bw --op mean --annotation tests/testbw2.bed bw2.mean >> test_run_out 2>&1
diff bw2.mean.all.tsv tests/testbw2.bed.mean

#test bigwig2min
time ./bcd_runner test.bam.all.bw --op min --annotation tests/testbw2.bed bw2.min >> test_run_out 2>&1
diff bw2.min.all.tsv tests/testbw2.bed.min

#test bigwig2max
time ./bcd_runner test.bam.all.bw --op max --annotation tests/testbw2.bed bw2.max >> test_run_out 2>&1
diff bw2.max.all.tsv tests/testbw2.bed.max
