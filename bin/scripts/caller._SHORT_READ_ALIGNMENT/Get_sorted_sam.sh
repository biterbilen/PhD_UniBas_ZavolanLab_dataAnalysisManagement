#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file_$TASK_ID
#$ -t 1-runtimes

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> log_file_$SGE_TASK_ID
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

PATH=$HOME/bin:$PATH

samples=samples_file
outdir=output_directory
pushd $outdir;

id=TopHat_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`
if [ -e $id/*sam ]; then echo $id/*sam exists; exit; fi

echo "Indexing and Sorting && bam2sam converting $id/accepted_hits.bam"
samtools index $id/accepted_hits.bam
samtools sort $id/accepted_hits.bam $id/accepted_hits.sorted
samtools view -H $id/accepted_hits.sorted.bam > $id/accepted_hits.sorted.sam
samtools view $id/accepted_hits.sorted.bam >> $id/accepted_hits.sorted.sam
echo DONE
