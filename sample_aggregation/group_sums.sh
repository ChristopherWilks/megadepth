#!/bin/bash -x
set -o pipefail -o nounset -o errexit 

manifest=$1
output_sentinel=$2
convert_to_int=$3
#run one group job at a time, but run it for all of the files listed in the manifest
if [[ -n $convert_to_int ]]; then
    cut -f 3 $manifest | xargs -n 1 -P 1 -I{} sh -c 'cut -f 4 $1 | sed "s/\.0*$//" > ${1}.unc' -- {}
else
    cut -f 3 $manifest | xargs -n 1 -P 1 -I{} sh -c 'cut -f 4 $1 > ${1}.unc' -- {}
fi
cut -f 3 $manifest | perl -ne 'chomp; $f=$_.".unc"; @s=stat($f); if($s[7]==0) { `cp blank_exon_sums $f`;}'
#when done, write a "done" sentinel
touch $output_sentinel
