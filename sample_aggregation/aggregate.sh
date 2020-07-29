mkdir -p paste

/bin/bash -x find_new_files.sh ./sums srav3h.ids.all.final.tsv paste all .tsv > find.run 2>&1 &

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "/
bin/bash -x sample_aggregation/group_sums.sh paste/all.$s.$r.manifest paste/all.$s.$r.grouped convert_to_int > runs/all.$s.$r.grouped.run 2>&1\n";
' | sort -u > group_sums.u.jobs

/usr/bin/time -v parallel -j 33 < group_sums.u.jobs > group_sums.u.jobs.run 2>&1

ls paste/all.??.??.manifest | perl -ne 'chomp; print "./remove_extra_samples.sh $_\n";' > remove_extra_samples.jobs

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "/
bin/bash -x sample_aggregation/paste_sums.sh paste/all.$s.$r.manifest paste/all.$s.$r.pasted 0 > runs/all.$s.$r.pasted.run 2>&1\n";' | sort -u > p
aste_sums_per_group.u.jobs

/usr/bin/time -v parallel -j33 < paste_sums_per_group.u.jobs > paste_sums_per_group.u.jobs.run 2>&1

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "$s\t$r\n";' | sort -u | perl -ne 'chomp;  ($s,$r)=split(/\t/,$_); `ls paste/all.$s.??.pasted > paste/all.$s.pasted.files.list`;'
