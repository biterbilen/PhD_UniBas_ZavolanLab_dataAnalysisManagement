#!/bin/bash
#$ -S /bin/bash
#$ -q long_m8
#$ -pe MT8 4
#$ -e stderr_$TASK_ID
#$ -o stdout_$TASK_ID
#$ -j n
#$ -N MARA_prep_job-1
#$ -cwd
#$ -l mem_total=6000M
#$ -t 1-4

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1"
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

#Usage: 
#nohup ./Get_uniq_mRNA_mapped_seq_MARA_bedFormat_extended_with_copies.sh 63 &> 63.log &
SGE_TASK_ID=1;
sampleIds=(502 503)
id=${sampleIds[$((SGE_TASK_ID-1))]}
annot=mRNA
uniq=1

mkdir -p bed_files/$id; pushd bed_files/$id

mapped_sequences=~mirz/clipzServer/data/database/$id/t_mapped_sequence 
genome_mappings=~mirz/clipzServer/data/cluster/$id/genome.sorted

if [ ! -e $mapped_sequences ]; then
	echo Database dump for t_mapped_sequence_$id
	mapped_sequences=t_mapped_sequence
	mysql -u mirzadmin --password=wUP4gZMW -h wnzmsql01 -D clipz -A > $mapped_sequences <<END
        select * from t_mapped_sequence_$id;
        quit
END
fi


~/BITER/scripts/selectAnnot.pl $mapped_sequences $genome_mappings $annot $uniq | sort -k 1,1 > genome_mappings 2> sample_size

less $mapped_sequences | sort -k 3,3 > mapped_sequences 

for try in `seq 1 3`; do
	join -1 3 -2 1 -t "`echo -e '\t'`" genome_mappings mapped_sequences \
		| awk 'BEGIN{OFS="\t";} { print $3,$7,$8,$2,$16,$6; }' \
		| ~/BITER/scripts/extend_bed_fromScoreField.pl 0 \
		> joined 2> LOG
	if [ ! -s LOG ]; then #zero size
		echo SUCCESS
		break;
 	fi
	echo RUN $try failed. Trying again... 
	try=$((try+1))
done

rm mapped_sequences genome_mappings
if [ ! -e LOG ]; then
	echo FAILED
else
	popd
	mv bed_files/$id/joined bed_files/$id.bed
	rm -rf bed_files/$id
fi
