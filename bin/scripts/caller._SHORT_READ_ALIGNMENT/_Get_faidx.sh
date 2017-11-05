#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
# $ -o LOG._make_faidx
#$ -o LOG._make_faidx_mm9

PATH=$HOME/bin:$PATH

bti_dir=~bilebi00/DATA/Bowtie_DB_INDEX
db_name=hg18 #hg_TR #TODO set here
db_name=mm9 #hg_TR #TODO set here

#in=/import/bc2/home/zavolan/GROUP/miRNA/RefSeq/DB_10-01-13/hg/hg_TR.fa
#in=/import/bc2/home/zavolan/bilebi00/DATA/hg18/chr19.fa
in=~bilebi00/DATA/$db_name/*.fa
out=$db_name.fa

mkdir -p $bti_dir; pushd $bti_dir;
cat $in > $out

cmd="samtools faidx $out"

echo $cmd;
$cmd

ls -latr $db_name*

