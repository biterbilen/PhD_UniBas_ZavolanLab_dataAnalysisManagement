#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
# $ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_junction_stats.chr19.rep_A

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

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

bedfile=~bilebi00/_SHORT_READ_ALIGNMENT/rep_A.hg18_chr19/MapSplice/best_junction.bed
tittag=MapSplice
#bedfile=~bilebi00/_SHORT_READ_ALIGNMENT/rep_A.hg18_chr19/tophat/junctions.bed
#tittag=TopHat
#bedfile=~bilebi00/_SHORT_READ_ALIGNMENT/rep_A.hg18/MapSplice/best_junction.bed
#tittag=MapSplice
tag="^JUNC"
inp=$bedfile.inp
cls=$inp.cls
out=$inp.out
stat=$inp.stat

overlap_fracs="0.4 0.6 0.8 1"

less ~bilebi00/_SHORT_READ_ALIGNMENT/data/coordinates/INTRON/chr* > $inp
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/bed2clsformat.pl $bedfile >> $inp

~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls 
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/fine_tune_junction_overlaps.pl $cls "$tag" > $out

#TODO consider intron length too
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/summ_overlap_stats.pl $out "$overlap_fracs" > $stat

~bilebi00/bin/R --no-save --args $stat $tittag < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/plot_junction_stats.R 


