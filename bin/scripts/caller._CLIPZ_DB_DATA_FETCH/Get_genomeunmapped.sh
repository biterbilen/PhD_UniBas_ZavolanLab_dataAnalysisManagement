#!/bin/bash

#===============================================================================
#       USAGE: #run it max 3 at a time; main script submits to clusters
#       PURPOSE: 
#       REQUIREMENTS:  
#       AUTHOR:  Biter Bilen
#       COMPANY:  Biozentrum, University of Basel
#===============================================================================

#TODO set
libs=(271 272 651 647 1037 1039)
t_organism_id=1 #human

SGE_TASK_ID=$1 #set to run in parallel

if [ $SGE_TASK_ID -gt ${#libs[@]} ]; then
	echo "$SGE_TASK_ID > ${#libs[@]}"
	exit;
fi

#make directory
soutdir=GenomeUnmapped
mkdir -p $soutdir; pushd $soutdir;

#get unique aligmnets from the mysql database with the header lines
#for SGE_TASK_ID in `seq 1 ${#libs[@]}`; do
	id=${libs[$((SGE_TASK_ID-1))]}
	mapped=$id.tab

	if [ ! -e $mapped.gz ]; then
		echo Database dump for t_mapped_sequence_$id t_sequence_genome_alignment_$id t_chromosome_contig
		date
		#/import/bc2/soft/bin/mysql --quick -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped <<END
		/import/bc2/soft/bin/mysql --quick -u clipzRead -h wnzmsql01 -D clipz --password=readClipz -A > $mapped <<END
SELECT id,sequence,copies
	FROM t_mapped_sequence_$id JOIN 
		(SELECT distinct t_sequence_anno_alignment_$id.t_sequence_id FROM t_sequence_anno_alignment_$id LEFT JOIN t_sequence_genome_alignment_$id
		ON t_sequence_anno_alignment_$id.t_sequence_id = t_sequence_genome_alignment_$id.t_sequence_id WHERE t_sequence_genome_alignment_$id.t_sequence_id IS NULL) as tsa
	WHERE id=tsa.t_sequence_id;
quit
END
date
#SELECT distinct t_mapped_sequence_$id.id,t_mapped_sequence_$id.sequence,t_mapped_sequence_$id.copies
#	FROM t_sequence_anno_alignment_$id LEFT JOIN t_sequence_genome_alignment_$id ON 
#		t_sequence_anno_alignment_$id.t_sequence_id = t_sequence_genome_alignment_$id.t_sequence_id,t_mapped_sequence_$id
#	WHERE t_sequence_genome_alignment_$id.t_sequence_id IS NULL and t_mapped_sequence_$id.id=t_sequence_anno_alignment_$id.t_sequence_id;
#quit




	gzip $mapped
	else
		echo Database dump done before
	fi
	echo $id DONE
#done


