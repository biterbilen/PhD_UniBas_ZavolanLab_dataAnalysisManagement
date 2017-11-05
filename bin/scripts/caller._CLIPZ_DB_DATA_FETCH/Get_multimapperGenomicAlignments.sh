#!/bin/bash
#run it max 3 at a time; main script submits to clusters

annot=tRNA
t_organism_id=1 #human
libs=(1598 1635 1628 1629 1637 1638 1639 1640 308 309 271 272 651 647 1037 1039  1598 211 212 213 1012 1019 1020 1021 1597 433 434 435 436 472 858 238 239 445 446 294 295 478 479 480 485 482 483 481 484 502 503 299 288 289 984 1341 1342 240 314 315)
SGE_TASK_ID=$1 #set to run in parallel

if [ $SGE_TASK_ID -gt ${#libs[@]} ]; then
	echo "$SGE_TASK_ID > ${#libs[@]}"
	exit;
fi

#make directory
soutdir=MultimapperGenomicAlignments
mkdir -p $soutdir; pushd $soutdir;

#get unique aligmnets from the mysql database with the header lines
#for SGE_TASK_ID in `seq 1 ${#libs[@]}`; do
	id=${libs[$((SGE_TASK_ID-1))]}
	mapped=$id.$annot.bedplus

	if [ ! -e $mapped.gz ]; then
		echo Database dump for t_mapped_sequence_$id t_sequence_genome_alignment_$id t_chromosome_contig
		/import/bc2/soft/bin/mysql --quick -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped <<END
select 
	name,chr_start-1,chr_end,t_mapped_sequence_$id.id,copies/genome_count_total,strand,alignment,sequence,annotation,genome_error,anno_error,five_utr_fract,cds_fract,three_utr_fract,exon_fract,intron_fract,unknown_fract,seq_len,genome_count_total
	from t_mapped_sequence_$id,t_sequence_genome_alignment_$id,t_chromosome_contig
	where t_mapped_sequence_$id.id=t_sequence_id and t_organism_id=$t_organism_id and t_chromosome_contig.id=t_chromosome_contig_id and annotation="$annot";
quit
END
	gzip $mapped
	else
		echo Database dump done before
	fi
	echo $id DONE
#done


