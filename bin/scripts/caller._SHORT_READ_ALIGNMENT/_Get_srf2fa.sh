#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
# $ -q fs_short@@qc_nehalem
#$ -l h_vmem=22G
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -t 1-4
#$ -o calisiyormusun

#usage
# nohup extract_from_srf.sh 1 &> 1.log &";
# nohup extract_from_srf.sh 2 &> 2.log &";

SGE_TASK_ID=$1 #if the files are in the archive, submit one by one since the clusters cannot see the archive directory
samples=~bilebi00/_SHORT_READ_ALIGNMENT/data/Zavolan_samples_srf
inp=`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1`
of=`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`
od=~bilebi00/_SHORT_READ_ALIGNMENT/data/Zavolan/
mkdir -p $od

#-C masks out bad quality sequences
~bilebi00/bin/srf2fasta -C $inp > $od/$of.fasta
~bilebi00/bin/srf2fastq -C $inp > $od/$of.fastq

echo DONE
