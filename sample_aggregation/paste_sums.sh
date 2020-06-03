#!/bin/bash
set -o pipefail -o errexit 

manifest=$1
output=$2
threads=$3
dont_get_ids=$4
existing_sums=$5
#[optional] number of threads to use when pigz'ing output, default 4
if [[ -z $threads ]]; then
	threads=1
fi

if [ -n "$dont_get_ids" ]; then
#if set, we're doing a final paste of all previously pasted (or copied) sample groups so no need to handle sample IDs
    #for debugging, dont delete pasted intermediates
    cat $manifest | perl -ne 'chomp; $f=$_; $s.=" $f"; $c++; END { if($c > 1) { print "paste $s\n"; `paste $s > '${output}'`; `rm $s`; } else { `cat $f > '${output}'`; }}'
    #cat $manifest | perl -ne 'chomp; $f=$_; $s.=" $f"; $c++; END { if($c > 1) { print "paste $s\n"; `paste $s > '${output}'`; `rm $s`; } else { `cat $f > '${output}'`; `rm $f`; }}'
else
#here we need to make sure we output the correct order of samples IDs as the column header
#acc_header is an env var
#optional flag to switch to accessions instad of rail_ids as the header
    if [[ -z $acc_header ]]; then
    	cut -f 2,3 $manifest | perl -ne 'chomp; ($rid,$f)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f.unc"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; } else { `cat $f.unc >> '${output}'`;  }}'
    	#cut -f 2,3 $manifest | perl -ne 'chomp; ($rid,$f)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f.unc"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; `rm $s`; } else { `cat $f.unc >> '${output}'`; `rm $f.unc`; }}'
    else
    	cut -f 3,4 $manifest | perl -ne 'chomp; ($f,$rid)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f.unc"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; } else { `cat $f.unc >> '${output}'`; }}'
    	#cut -f 3,4 $manifest | perl -ne 'chomp; ($f,$rid)=split(/\t/,$_); $h.="$rid\t"; $s.=" $f.unc"; $c++; END { $h=~s/\t$//; open(OUT,">'${output}'"); print OUT "$h\n"; close(OUT); if($c > 1) { print "paste $s\n"; `paste $s >> '${output}'`; `rm $s`; } else { `cat $f.unc >> '${output}'`; `rm $f.unc`; }}'
    fi
fi


if [ -n "$existing_sums" ]; then
    mv ${output} ${output}.pre_existing
    if [[ $threads -eq 0 ]]; then
        #assumes if we dont want to compress the output then the intermediate input is no compressed either
        paste ${existing_sums} ${output}.pre_existing > ${output}
    else
        paste ${existing_sums} ${output}.pre_existing | pigz --fast -p $threads > ${output}
    fi
fi
