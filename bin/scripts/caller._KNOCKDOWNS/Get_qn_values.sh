#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1" >> log_file
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=~bilebi00/bin:$PATH

outdir=output_directory
samples="samples"

#TODO replace me later
outdir=Normalization
samples=~bilebi00/_KNOCKDOWNS/data/KnockDowns/*/transcript_expression

mkdir -p $outdir; cd $outdir

file=`~bilebi00/_KNOCKDOWNS/scripts/outerJoin.pl 1 0 0 3 $samples`;
R --no-save --args $file < ~bilebi00/_KNOCKDOWNS/scripts/qn.R
matlab -nojvm -r "addpath ~bilebi00/CHRISTIANE/scripts/2ndGen; matrixplotCorrelations('$file.raw.image.png','$file.raw.scatter.png','$file.raw','Raw Expr'); exit;";
matlab -nojvm -r "addpath ~bilebi00/CHRISTIANE/scripts/2ndGen; matrixplotCorrelations('$file.qn.image.png','$file.qn.scatter.png','$file.qn','QN Expr'); exit;";


