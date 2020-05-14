#!/usr/bin/env bash
export LD_LIBRARY_PATH=/data7/megadepth/libBigWig:/data7/megadepth/htslib:$LD_LIBRARY_PATH
MD=../mosdepth
BT=~/bwtool_mine2/bwtool
B2B='python3 ../bam2bigwig.py'
KB2B=/data/kent_tools/bigWigToBedGraph
WT=wiggletools
EXONS=/data7/megadepth/exons.bed

EXONS=$1
BAM=$2

echo -n "" > ${BAM}.per-base.bed.gz.bw.aucs
echo -n "" > ${BAM}.annotated_aucs
for t in default unique; do
    qarg=`perl -e '$t='${t}'; print "0" if($t eq "default"); print "10" if($t eq "unique");'`
    time $MD -F260 -t 4 -Q${qarg} ${BAM}.${t} $BAM
    time cat <(samtools view -H $BAM) <(zcat ${BAM}.${t}.per-base.bed.gz) | $B2B ${BAM}.${t}.per-base.bed.gz.bw
    time $BT summary $EXONS ${BAM}.${t}.per-base.bed.gz.bw /dev/stdout -fill=0 -with-sum -keep-bed -decimals=0 | cut -f1-3,10 | sort -k1,1 -k2,2n -k3,3n > ${BAM}.${t}.per-base.bed.gz.bw.sums.bed
    $KB2B ${BAM}.${t}.per-base.bed.gz.bw ${BAM}.${t}.per-base.bed.gz.bw.sums.bed.bg
    $WT AUC ${BAM}.${t}.per-base.bed.gz.bw | perl -ne 'BEGIN { $t='${t}'; } chomp; $f=$_; $f=~s/\.0+$//; print "ALL_READS_ALL_BASES\t$f\n" if($t eq "default"); print "UNIQUE_READS_ALL_BASES\t$f\n" if($t eq "unique");' >> ${BAM}.per-base.bed.gz.bw.aucs
    cat ${BAM}.${t}.per-base.bed.gz.bw.sums.bed | perl -ne 'BEGIN { $t='${t}'; } chomp; ($c,$s,$e,$v)=split(/\t/,$_); $f+=$v; END { print "ALL_READS_ANNOTATED_BASES\t$f\n" if($t eq "default"); print "UNIQUE_READS_ANNOTATED_BASES\t$f\n" if($t eq "unique"); }' >> ${BAM}.annotated_aucs
done
cat ${BAM}.per-base.bed.gz.bw.aucs >> ${BAM}.annotated_aucs
