#!/bin/bash
#===============================================================================
#       USAGE: #run it max 3 at a time; main script submits to clusters
#       PURPOSE: 
#       REQUIREMENTS:  
#       AUTHOR:  Biter Bilen
#       COMPANY:  Biozentrum, University of Basel
#===============================================================================

#TODet
libs=(271 272 651 647 1037 1039 3258 3262 3263 3268 3027 3029 3031 3033 271 272 651 647 248 238 288 289 309 1598 1628 1629 299 445 446 308 240)


SGE_TASK_ID=$1 #set to run in parallel

if [ $SGE_TASK_ID -gt ${#libs[@]} ]; then
	echo "$SGE_TASK_ID > ${#libs[@]}"
	exit;
fi

#make directory
soutdir=Unmapped
mkdir -p $soutdir; pushd $soutdir;

#get unique aligmnets from the mysql database with the header lines
#for SGE_TASK_ID in `seq 1 ${#libs[@]}`; do
	id=${libs[$((SGE_TASK_ID-1))]}
	mapped=$id.tab

	if [ ! -e $mapped.gz ]; then
		echo Database dump for t_mapped_sequence_$id t_sequence_genome_alignment_$id t_chromosome_contig
		/import/bc2/soft/bin/mysql --quick -u clipzRead -h wnzmsql01 -D clipz --password=readClipz -A > $mapped <<END
select * from t_unmapped_sequence_$id;
quit
END
	gzip $mapped
	else
		echo Database dump done before
	fi
	echo $id DONE
#done


