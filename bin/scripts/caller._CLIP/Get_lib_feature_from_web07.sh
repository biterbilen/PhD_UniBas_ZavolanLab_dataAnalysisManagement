#!/bin/bash
SGE_TASK_ID=$1 #TODO doesn't work from clusters

#TODO set
wiggle=0
MTsummary=0

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
outdir=Analysis
soutdir=UniqueRawData_ALL_woContamination
#soutdir=UniqueRawData_mRNA_repeat_none

mkdir -p $outdir/$soutdir; pushd $outdir/$soutdir;

id=`head -n $SGE_TASK_ID $sid_prot_f | tail -n 1 | cut -f 2`;

odir=DB_all_$id
if [ ! -e $odir ] && [ $wiggle == 1 ]; then
	rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$odir.tar.gz .
	tar -zxvf $odir.tar.gz
	rm -rf $odir.tar.gz
else
	echo $odir.tar.gz or wiggle=$wiggle-0 exists
fi

#summary files
if [ "`echo $soutdir | grep ALL`" != "" ]; then
	rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$id.annotation_summary .
	rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$id.general_summary .
	#mutations statistics
	#sums the copies*length for each entry
	of=$id.mutation_summary
	if [ ! -e $of ] && [ $MTsummary == 1 ]; then
		cutoff=0
		echo -e "mutation\tsum(copies)\tcount(*)\tcopies_cutoff" > $of
		for i in `ls $odir | grep -v copies`; do echo "$i" | cut --bytes=-3; done | \
			sort | uniq | while read f; do
				less $odir/$f* | \
					awk 'BEGIN{c=0;k=0;}{OFS="\t"; if($4>cutoff) { c=$4*($3-$2)+c; k=k+($3-$2); } if (NR==nr) { print f,c,k,cutoff;}}' f=$f cutoff=$cutoff nr=`cat $odir/$f* | wc -l`;
			done | sort -k 2,2gr >> $of
	else
		echo $f exists;
	fi
else
	rsync -t mirz@bc2-web08.bc2.unibas.ch:~/BITER/$soutdir/$id.effective_sample_size .
fi

echo $SGE_TASK_ID  $id is DONE
