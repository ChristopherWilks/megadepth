#!/usr/bin/env bash

time ./bamcount tests/test.bam --threads 4 --coverage --no-head --bigwig test.bam --auc test.bam --min-unique-qual 10 --annotation tests/test_exons.bed test.bam --frag-dist test.bam --alts test.bam > test_run_out 2>&1

diff tests/test_run_out.txt test_run_out

diff tests/test.bam.orig.frags.tsv test.bam.frags.tsv
diff tests/test.bam.orig.alts.tsv test.bam.alts.tsv
for f in all unique; do
    #need to write an option to output bigwig instead of relying on kent tools for this test
    #for t in tsv bw.bg; do
    #    diff <(zcat tests/test.bam.mosdepth.${f}.per-base.${t}) test.bam.${f}.${t}
    #done
    #could do wc on the bw for now, but maybe not
    #diff <(wc -c tests/test.bam.mosdepth.${f}.per-base.bw | cut -d' ' -f 1) <(wc -c test.bam.${f}.bw| cut -d' ' -f 1)
    diff tests/test.bam.mosdepth.${f}.per-base.exon_sums.tsv <(cat test.bam.${f}.tsv | sort -k1,1 -k2,2n -k3,3n)
done
diff tests/test.bam.mosdepth.bwtool.all_aucs test.bam.auc.tsv


