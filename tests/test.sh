#!/usr/bin/env bash
set -xe

static=$1

if [[ -z $static ]]; then
    time ./megadepth http://stingray.cs.jhu.edu/data/temp/test.bam --prefix test.bam --threads 4 --bigwig --auc --min-unique-qual 10 --annotation tests/test_exons.bed --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout --no-auc-stdout --filter-out 260 --add-chr-prefix human > test_run_out 2>&1

    time ./megadepth http://stingray.cs.jhu.edu/data/temp/test.cram --prefix test.cram --threads 4 --coverage --no-coverage-stdout --auc --min-unique-qual 10 --annotation 400 --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout --no-auc-stdout --filter-out 260 > test_cram_run_out 2>&1
else
    time ./megadepth tests/test_noprefix.bam --prefix test.bam --threads 4 --bigwig --auc --min-unique-qual 10 --annotation tests/test_exons.bed --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout --no-auc-stdout --filter-out 260 --add-chr-prefix human > test_run_out 2>&1
    time ./megadepth tests/test.cram --prefix test.cram --threads 4 --coverage --no-coverage-stdout --auc --min-unique-qual 10 --annotation 400 --frag-dist --alts --include-softclip --only-polya --read-ends --test-polya --no-annotation-stdout --no-auc-stdout --filter-out 260 > test_cram_run_out 2>&1
fi


diff <(sort tests/test.bam.orig.frags.tsv) <(sort test.bam.frags.tsv)
diff tests/test.bam.orig.alts.tsv test.bam.alts.tsv
diff tests/test.bam.orig.softclip.tsv test.bam.softclip.tsv
for f in annotation unique; do
    diff tests/test.bam.mosdepth.${f}.per-base.exon_sums.tsv test.bam.${f}.tsv
done
diff tests/test.bam.mosdepth.bwtool.all_aucs test.bam.auc.tsv

#test base coverage other than BigWigs
diff tests/test.cram.coverage.tsv test.cram.coverage.tsv

#test --annotation <window_bp_length>
cut -f 1,4 test.cram.window.tsv | perl -ne 'chomp; ($c,$v)=split(/\t/,$_); $h{$c}+=$v; END { for $c (sort keys %h) { print "$c\t".$h{$c}."\n"; }}' > test.cram.window.summed.tsv
diff tests/test.cram.window.summed.tsv test.cram.window.summed.tsv

#check --op mean with BAMs
./megadepth tests/test.bam --annotation tests/test_exons.bed --op mean --add-chr-prefix human > test.bam.mean
paste <(cut -f 4 test.bam.annotation.tsv) <(cut -f 2- test.bam.mean) | perl -ne 'chomp; $f=$_; ($sum,$s,$e,$m)=split(/\t/,$_); $d=($e-$s); $m2=$sum/$d; $m2=sprintf("%.2f",$m2); if($m != $m2) { print "$f\n"; $ret=1;} END { exit($ret); }'

./megadepth tests/test.bam | fgrep "ALL_READS_ALL_BASES" > auc.single
diff auc.single <(fgrep "ALL_READS_ALL_BASES" tests/test.bam.mosdepth.bwtool.all_aucs)

cat test.bam.starts.tsv test.bam.ends.tsv | sort -k1,1 -k2,2n -k3,3n > test_starts_ends.tsv
diff test_starts_ends.tsv <(sort -k1,1 -k2,2n -k3,3n tests/test.bam.read_ends.both.unique.tsv)

time ./megadepth tests/test2.bam --threads 4 --junctions --all-junctions --prefix test2.bam > test2_run_out 2>&1

diff tests/test2.bam.jxs.tsv test2.bam.jxs.tsv
diff tests/test2.bam.all_jxs.tsv test2.bam.all_jxs.tsv

#test just total auc
time ./megadepth test.bam.all.bw | grep "AUC" > test.bw1.total_auc
diff test.bw1.total_auc tests/testbw1.total_auc

#test bigwig2sums/auc
time ./megadepth test.bam.all.bw --annotation tests/testbw1.bed --auc --prefix test.bam.bw1 --no-annotation-stdout --no-auc-stdout
diff test.bam.bw1.annotation.tsv tests/testbw1.bed.out.tsv
diff test.bam.bw1.auc.tsv tests/testbw1.annot_auc

##use different order in BED file from what's in BW to test keep_order == true
time ./megadepth test.bam.all.bw --annotation tests/testbw2.bed --auc --prefix test.bam.bw2 --no-annotation-stdout --no-auc-stdout
diff test.bam.bw2.annotation.tsv tests/testbw2.bed.out.tsv
diff test.bam.bw2.auc.tsv tests/testbw2.annot_auc

#test bigwig2mean
time ./megadepth test.bam.all.bw --op mean --annotation tests/testbw2.bed --prefix bw2.mean --no-annotation-stdout >> test_run_out 2>&1
diff bw2.mean.annotation.tsv tests/testbw2.bed.mean

#test bigwig2min
time ./megadepth test.bam.all.bw --op min --annotation tests/testbw2.bed --prefix bw2.min --no-annotation-stdout >> test_run_out 2>&1
diff bw2.min.annotation.tsv tests/testbw2.bed.min

#test bigwig2max
time ./megadepth test.bam.all.bw --op max --annotation tests/testbw2.bed --prefix bw2.max --no-annotation-stdout >> test_run_out 2>&1
diff bw2.max.annotation.tsv tests/testbw2.bed.max

#now test same-start alignments for overlapping pairs
./megadepth tests/test3.bam --auc --coverage --prefix t3 --no-auc-stdout > t3.tsv
diff <(head -197 tests/test3.out.tsv) t3.tsv
diff <(tail -n1 tests/test3.out.tsv) t3.auc.tsv

#with uniques
./megadepth tests/test3.bam --coverage  --min-unique-qual 10 --bigwig --auc --prefix test3 --no-auc-stdout
diff tests/test3.auc.out.tsv test3.auc.tsv

#long reads support for junctions
./megadepth tests/long_reads.bam --junctions --prefix long_reads.bam --long-reads
diff tests/long_reads.bam.jxs.tsv long_reads.bam.jxs.tsv

#test bigwig2sum on remote BW
if [[ -z $static ]]; then
    time ./megadepth http://stingray.cs.jhu.edu/data/temp/megadepth.test.bam.all.bw --op mean --annotation tests/testbw2.bed --prefix bw2.remote.mean --no-annotation-stdout >> test_run_out 2>&1
    diff bw2.remote.mean.annotation.tsv tests/testbw2.bed.mean
fi

##only print sums use different order in BED file from what's in BW to test keep_order == true
time ./megadepth test.bam.all.bw --sums-only --annotation tests/testbw2.bed --prefix test.bam.bw2 > test.bam.bw2.annotation.tsv
diff test.bam.bw2.annotation.tsv <(cut -f 4 tests/testbw2.bed.out.tsv)

./megadepth http://stingray.cs.jhu.edu/data/temp/test.bam --prefix test.bam.names --threads 4 --alts --write-names --include-softclip --only-polya --test-polya --no-annotation-stdout --no-auc-stdout --filter-out 260 --add-chr-prefix human > test_run_out2 2>&1
diff test.bam.names.alts.tsv tests/test.bam.names.alts.tsv

#test multiple overlap types for BigWig annotation processing
./megadepth tests/bw.all_overlap_types.test_input.bw --annotation tests/gh_bug_9.bed --auc > test.bw.all_overlap_types.test_output.bed
diff test.bw.all_overlap_types.test_output.bed tests/bw.all_overlap_types.test_output.bed

#test faster mode with collapsed intervals in BigWig annotation processing
./megadepth tests/TCGA_BLCA_A13J.vcf.gz_cg_cov5.bw.bg.gz.chr1.60379.62229.bw --annotation tests/chr1.61863.62160.bed --no-annotation-stdout --prefix TCGA_BLCA_A13J_vs_chr1.61863.62160
diff TCGA_BLCA_A13J_vs_chr1.61863.62160.annotation.tsv tests/TCGA_BLCA_A13J_vs_chr1.61863.62160.annotation.tsv

#test that we're catching out of order chromosomes
output=$(./megadepth tests/TCGA_BLCA_A13J.vcf.gz_cg_cov5.bw.bg.gz.chr1.60379.62229.bw --annotation tests/chr1.61863.62160.bad_chrm_order.bed 2>&1)
fgrep "falling back" <(echo "$output")

#clean up any previous test files
rm -f test*tsv test*auc bw2* test3* test2* t3.* long_reads.bam.jxs.tsv test_run_out *null*.unique.tsv test.*.bw auc.single test.bam.mean test.cram.coverage.tsv test_cram_run_out test.cram.coverage.tsv.summed test.bw.all_overlap_types.test_output.bed TCGA_BLCA_A13J_vs_chr1.61863.62160.annotation.tsv
