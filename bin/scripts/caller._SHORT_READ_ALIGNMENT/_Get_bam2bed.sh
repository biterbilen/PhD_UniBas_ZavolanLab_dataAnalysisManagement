#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_bam2bed.$TASK_ID
#$ -t 1-3

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH

outdir=Shivendra/MARA

mkdir -p $outdir; cd $outdir;

bamfs=(mRNAseq_siCTRL-A_HeLa mRNAseq_siEWSR1_HeLa totRNAseq-Hela-0)

bamf=${bamfs[$(($SGE_TASK_ID - 1))]}


#bamToBed -split -i $bamf.bam | gzip -c > $bamf.bed.gz
bamToBed -tag NH -i $bamf.bam | awk '$5==1{ print; }' | gzip -c > $bamf.uniq.bed.gz

echo DONE




