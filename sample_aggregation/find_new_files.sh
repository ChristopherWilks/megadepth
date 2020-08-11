#!/bin/bash
set -o pipefail -o errexit 

#top level of incoming dir
#(where fully done analyses were copied/symlinked from original set of attempts)
search_dir=$1
#map between external sample accessions/UUIDs and study format: study<TAB>sample<TAB>secondary_sample_id (or just sample)
sample_ids_file=$2
destination_dir=$3
#e.g. "all"
analysis_type=$4
#.e.g ".tsv"
additional_suffix=$5

#now uses both study loworder digits *and* run loworder digits in manifest file names
find -L $search_dir -name "*.${analysis_type}${additional_suffix}" -size +0 | perl -ne 'BEGIN { open(IN,"<'${sample_ids_file}'"); %rids; while($line=<IN>) { chomp($line); @f=split(/\t/,$line); $rid=$f[2]; $run_acc=$f[1]; $rids{$run_acc}=[$rid, $f[0]]; } close(IN); } chomp; $f=$_; @f=split(/\//,$f); $fname=pop(@f); `pushd '${destination_dir}' && ln -fs ../$f && popd`; $fname=~/^([^\.]+)\./; $run_acc=$1; $run_acc=~/(..)$/; $acc_loworder=$1; ($rid, $study)=@{$rids{$run_acc}}; $study=~/(..)$/; $study_loworder=$1; push(@{$h{$study_loworder.".".$acc_loworder}},["'${destination_dir}'/$fname",$rid,"'${destination_dir}'/$fname",$run_acc]); END { open(ALL_OUT,">'${destination_dir}'/'${analysis_type}'.groups.manifest"); for $k (keys %h) { open(OUT,">'${destination_dir}'/'${analysis_type}'.$k.manifest"); for $a (@{$h{$k}}) { print ALL_OUT "".$a->[3]."\t'${analysis_type}'.$k.manifest\n"; print OUT "".join("\t",@$a)."\n"; } close(OUT); } close(ALL_OUT); }'

#touch $destination_dir/${analysis_type}.groups.manifest
