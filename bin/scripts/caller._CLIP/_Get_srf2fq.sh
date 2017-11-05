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
ids=(1628_1637 1629_1638 308_1639 309_1640)

id=${ids[$((SGE_TASK_ID-1))]}

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
srff=`less $sid_prot_f | awk -F "\t" '$2==id{ print $13}' id=$id`
project=`less $sid_prot_f | awk -F "\t" '$2==id{ print $9}' id=$id`

outfile=Analysis/Project_$project/rawdata
mkdir -p $outfile; pushd $outfile

echo $srff $id
##-C masks out bad quality sequences
#~bilebi00/bin/srf2fastq -C $srff > $id.fastq
~bilebi00/bin/srf_info $srff > $id.info.log
gzip -f $id.fastq

echo DONE
