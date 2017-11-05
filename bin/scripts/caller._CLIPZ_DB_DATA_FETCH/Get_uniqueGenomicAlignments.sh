#!/bin/bash

#===============================================================================
#       USAGE: #run it max 3 at a time; main script submits to clusters
#       PURPOSE: 
#       REQUIREMENTS:  
#       AUTHOR:  Biter Bilen
#       COMPANY:  Biozentrum, University of Basel
#===============================================================================

#TODO set
libs=(4188 4189 4190 4191 186 183 3353 3354 3355 3356 3027 3042 3029 3044 3031 3045 3033 3046 3015 2722 2123 2124 321 1215 1216 1217 1218 1219 1221 1431 1106 1107 1108 1122 1117 1120 1121 1635 1598  1639 1640 1637 1638 1633 1634 1628 1629  1218 1219 1344  271 272 651 647 1037 1039  211 212 213 1012 1019 1020 1021 1597 433 434 435 436 472 858 238 239 445 446 294 295 478 479 480 485 482 483 481 484 502 503 299 300 308 309 288 289 984 1341 1342 240 314 315)
t_organism_id=1 #human

SGE_TASK_ID=$1 #set to run in parallel

if [ $SGE_TASK_ID -gt ${#libs[@]} ]; then
	echo "$SGE_TASK_ID > ${#libs[@]}"
	exit;
fi

#make directory
soutdir=UniqueGenomicAlignments
mkdir -p $soutdir; pushd $soutdir;

#get unique aligmnets from the mysql database with the header lines
#for SGE_TASK_ID in `seq 1 ${#libs[@]}`; do
	id=${libs[$((SGE_TASK_ID-1))]}
	mapped=$id.bedplus

	if [ ! -e $mapped.gz ]; then
		echo Database dump for t_mapped_sequence_$id t_sequence_genome_alignment_$id t_chromosome_contig
		#/import/bc2/soft/bin/mysql --quick -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped <<END
		/import/bc2/soft/bin/mysql --quick -u clipzRead -h wnzmsql01 -D clipz --password=readClipz -A > $mapped <<END
select 
	name,chr_start-1,chr_end,t_mapped_sequence_$id.id,copies,strand,alignment,sequence,annotation,genome_error,anno_error,five_utr_fract,cds_fract,three_utr_fract,exon_fract,intron_fract,unknown_fract,seq_len 
	from t_mapped_sequence_$id,t_sequence_genome_alignment_$id,t_chromosome_contig 
	where t_mapped_sequence_$id.id=t_sequence_id and genome_count_total=1 and t_organism_id=$t_organism_id and t_chromosome_contig.id=t_chromosome_contig_id;
quit
END
	gzip $mapped
	else
		echo Database dump done before
	fi
	echo $id DONE
#done


