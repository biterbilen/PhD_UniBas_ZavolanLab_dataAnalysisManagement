#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
# $ -q fs_short@@qc_nehalem
#$ -l h_vmem=22G
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -t 1-2
#$ -o calisiyormusun

#usage
# nohup extract_from_srf.sh 1 &> 1.log &";
# nohup extract_from_srf.sh 2 &> 2.log &";

SGE_TASK_ID=$1

#for fo in Sam68_KO Sam68_WT; do
for fo in Sam68_KO_brainStem Sam68_WT_brainStem; do
	#if [ $SGE_TASK_ID == 1 ] && [ "$fo" != "Sam68_KO" ]; then
	#	continue
	#elif [ $SGE_TASK_ID == 2 ] && [ "$fo" != "Sam68_WT" ]; then
	#	continue
	#fi
	f=~bilebi00/_NRXN/data/$fo/*srf
	#-C masks out bad quality sequences
	~bilebi00/bin/srf2fasta -C $f | gzip > ~bilebi00/_NRXN/data/$fo.fasta.gz
	~bilebi00/bin/srf2fastq -C $f | gzip > ~bilebi00/_NRXN/data/$fo.fastq.gz
done

echo DONE
