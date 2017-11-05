#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
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

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=output_directory
ftag=ftag

pushd $outdir

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all_clsWMs_$ftag.pdf $ftag/$ftag*cls_wm*.pdf;
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all_topWMwindows$ftag.pdf $ftag/*-$ftag.pdf;

echo DONE
