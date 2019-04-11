#!/usr/bin/env bash

rm -f test.bam.auc.tsv test.bam.all.tsv test.bam.unique.tsv test.bam.frags.tsv test.bam.alts.tsv test.bam.softclip.tsv

time ./bcd_runner tests/test.bam --threads 4 --coverage --no-head --bigwig test.bam --auc test.bam --min-unique-qual 10 --annotation tests/test_exons.bed test.bam --frag-dist test.bam --alts test.bam --include-softclip test.bam --only-polya --read-ends test --test-polya --keep-annot-order > test_run_out 2>&1

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
