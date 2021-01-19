#!/usr/bin/env bash
set -x -o pipefail -o errexit -o nounset
#this will produce a sorted, STAR-compatible formatted file of splice junction calls from BAM processed with Megadepth --all-junctions
#unlike STAR:
#1) splice junction strand is *not* determined (always 0 in column 4)
#2) splice junction motifs are *not* determineda (always 0 in column 5)
#both of these can be determined by extracting the dinucleotide motifs for the given reference coordinates for canonical motifs

#for faster sorting
export LC_ALL=C
#sort & format for compatibility unifier using STAR's SJ.out file format:
#chromosome  1based_intron_start 1based_intron_end   strand_0:undefined,1:+,2:-  intron_motif:0:non-canonical;1:GT/AG,2:CT/AC,3:GC/AG,4:CT/GC,5:AT/AC,6:GT/AT    0:unannotated,1:annotated_in_spliceDB   #uniquely_mapping_reads #multi_mapping_reads    maximum_spliced_alignment_overhang
#this determines maximum_spliced_alignment_overhang via STAR's method as described here (last post by Dobin):
#https://groups.google.com/g/rna-star/c/XN0cWBxVFcM/m/ywcUg_s3CQAJ
sort -k2,2 -k3,3n -k4,4n -k1,1 -u $jx_file | cut -f 1,2-4,6,7 | perl -ne 'chomp; ($qname,$c,$s,$e,$cigar,$is_unique)=split(/\t/,$_); @ns=split(/\d+N/,$cigar); if(scalar(@ns) > 2) { $n=$h{$qname}++; $cigar=$ns[$n]."1N".$ns[$n+1]; } $o=-1; if($cigar=~/(\d+)M\d+N(\d+)M/) { $overhangL=$1; $overhangR=$2; $o=$overhangL<$overhangR?$overhangL:$overhangR; } if($pc) { if($s == $ps && $e == $pe) { if($is_unique == 1) { $ucnt++; } else { $cnt++; } $po=$po>$o?$po:$o;  next; } else { print "$pc\t$ps\t$pe\t0\t0\t0\t$ucnt\t$cnt\t$po\n"; }} $ucnt=0; $cnt=0; if($is_unique == 1) { $ucnt=1; } else { $cnt=1; } $pc=$c; $ps=$s; $pe=$e; $po=$o; END { if($pc) { print "$pc\t$ps\t$pe\t0\t0\t0\t$ucnt\t$cnt\t$po\n"; }}' > ${jx_file}.sjout
