#!/usr/bin/env bash
set -xe

sdir=$(dirname $0)

mkdir -p paste
mkdir -p runs

/bin/bash -x ${sdir}/find_new_files.sh ./sums srav3h.ids.all.final.tsv ./paste all .tsv > find.run 2>&1 &

ls paste/all.??.??.manifest | perl -ne 'chomp; print "'${sdir}'/remove_extra_samples.sh $_\n";' > remove_extra_samples.jobs
/bin/bash -x remove_extra_samples.randoms100.jobs

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "/bin/bash -x '${sdir}'/group_sums.sh paste/all.$s.$r.manifest paste/all.$s.$r.grouped convert_to_int > runs/all.$s.$r.grouped.run 2>&1\n";' | sort -u > group_sums.u.jobs
/usr/bin/time -v parallel -j 33 < group_sums.u.jobs > group_sums.u.jobs.run 2>&1

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "/bin/bash -x '${sdir}'/paste_sums.sh paste/all.$s.$r.manifest paste/all.$s.$r.pasted 0 > runs/all.$s.$r.pasted.run 2>&1\n";' | sort -u > paste_sums_per_group.u.jobs
/usr/bin/time -v parallel -j33 < paste_sums_per_group.u.jobs > paste_sums_per_group.u.jobs.run 2>&1

cut -f 1-2 srav3h.ids.all.final.tsv | perl -ne 'chomp; ($study,$run)=split(/\t/,$_); $study=~/(..)$/; $s=$1; $run=~/(..)$/; $r=$1; print "$s\t$r\n";' | sort -u | perl -ne 'chomp;  ($s,$r)=split(/\t/,$_); print "ls paste/all.$s.??.pasted > paste/all.$s.pasted.files.list\n";' | sort -u > ls.paste.u.jobs
/bin/bash -x ls.paste.jobs.u > ls.paste.jobs.u.run 2>&1

ls paste/*.list | perl -ne 'chomp; $f=$_; $f=~/all\.(..)\.pasted\.files\.list/; $s=$1; print "/bin/bash -x '${sdir}'/paste_sums.sh $f paste/all.$s.pasted 0 dont_get_ids > runs/all.$s.pasted.run 2>&1\n";' | sort -u > paste_sums_per_study_group.u.jobs
/usr/bin/time -v parallel -j33 < paste_sums_per_study_group.u.jobs > paste_sums_per_study_group.u.jobs.run 2>&1

ls paste/all.??.pasted > paste/all.groups.pasted.files.list 

/usr/bin/time -v /bin/bash -x ${sdir}/paste_sums.sh paste/all.groups.pasted.files.list all.samples.pasted.gz 8 dont_get_ids ryten.random.bed
