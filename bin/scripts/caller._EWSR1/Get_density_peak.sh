#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-2

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

outdir=output_directory
indir=indir
sid_prot_file=sid_prot_file
tag=tag
idtag=idtag

pushd $outdir;

id=`less $sid_prot_file | cut -f 2 | head -n $SGE_TASK_ID | tail -n 1`;

soutdir=$tag$id;

if [ -e $soutdir ]; then
	exit;
fi

mkdir -p $soutdir; pushd $soutdir;

infile=$indir/$id/inputdata.tag.sites.sorted
densities=sorted_T2C_densities

~bilebi00/_EWSR1/scripts/window2densityPeakFormat.pl $idtag $infile > $densities
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $densities 5 > reverse_cum_$densities

