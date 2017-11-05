#!/bin/bash
SGE_TASK_ID=$1 #TODO doesn't work from clusters

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
outdir=Analysis/UniqueRawData_ALL_woContamination

mkdir -p $outdir; pushd $outdir;

id=`head -n $SGE_TASK_ID $sid_prot_f | tail -n 1 | cut -f 2`;
odir=DB_all_$id
minCov=1
minCov=0

f=$id.mutation_summary

echo -e "mutation\tsum(copies)\tcount(*)" > $f
for i in `ls $odir | grep -v copies`; do echo "$i" | cut --bytes=-3; done | \
	sort | uniq | while read mut; do
		nr=`cat $odir/$mut* | wc -l`
		cat $odir/$mut* | \
			awk 'BEGIN{c=0;k=0;OFS="\t";} { if ( $4 > 0 ) { c=$4*($3-$2)+c; k=k+($3-$2); } if (NR==nr) { print f,c,k; }}' f=$mut nr=$nr;
	done | sort -k 2,2gr >> $f

echo $SGE_TASK_ID  $id is DONE
