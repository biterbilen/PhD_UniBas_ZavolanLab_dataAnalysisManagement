#!/bin/bash
SGE_TASK_ID=$1 #TODO doesn't work from clusters

soutdir=UniqueRawData_ALL_woContamination

sid_prot_f=~bilebi00/_ARE/data/protein_sampleid_list
mkdir -p $outdir/$soutdir; pushd $outdir/$soutdir;

id=`head -n $SGE_TASK_ID $sid_prot_f | tail -n 1 | cut -f 2`;
odir=DB_all_$id
if [ ! -e DB_all_$id ]; then
	rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$odir.tar.gz .

	tar -zxvf $odir.tar.gz
	rm -rf $odir.tar.gz
fi

f=$id.mutation_summary
cutoff=0
echo -e "mutation\tsum(copies)\tcount(*)\tcopies_cutoff" > $f
for i in `ls $odir | grep -v copies`; do echo "$i" | cut --bytes=-3; done | \
	sort | uniq | while read f; do
		less $odir/$f* | awk 'BEGIN{c=0;k=0;}{OFS="\t"; if ($4>cutoff) { c=$4*($3-$2)+c; k=k+1; } if (NR==nr) { print f,c,cutoff; }}' f=$f cutoff=$cutoff nr=`cat $odir/$f* | wc -l`;
	done | sort -k 2,2gr >> $f

echo $SGE_TASK_ID  $id is DONE
