#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
# $ -q fs_long@@qc_nehalem 
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_junction_stats_of_GMAP_tophat

#TODO change log_file_name

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

#1.
sjtool=tophat;
infile=junctions.bed
#2.
#sjtool=MapSplice
#infile=best_junction.bed

bedfile=~bilebi00/_SHORT_READ_ALIGNMENT/rep_A.hg18_chr19/$sjtool/$infile
ftag=GMAP
tag="INTRON"
r=$RANDOM
inp=$r.inp
cls=$r.cls
out=$bedfile.$ftag.out
stat=$bedfile.$ftag.stat
overlapIndex=6;
coverageindex=5;
wHitRegCov=1;
wScore=1

overlap_fracs="1" #"0.4 0.6 0.8 1" Not much difference
re=1; #remove extension
intronsf=~bilebi00/_SHORT_READ_ALIGNMENT/data/coordinates/INTRON/chr*
expressed_trans_file=~bilebi00/_SHORT_READ_ALIGNMENT/data/expressed_transcripts

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/select_cls_format_lines_given_ids.pl "$intronsf" $expressed_trans_file > $inp  
#~bilebi00/_SHORT_READ_ALIGNMENT/scripts/bed2clsformat.pl $bedfile $re >> $inp
~bilebi00/_DIS3/scripts/track2clsFormat.pl $bedfile 0 0 0 $wScore >> $inp

~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls 
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/fine_tune_junction_overlaps.pl $cls "$tag" > $out

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/summ_overlap_stats.pl $out "$overlap_fracs" $overlapIndex $coverageindex $wHitRegCov > $stat

~bilebi00/bin/R --no-save --args $stat "$ftag SJs of Expressed Transcripts hit by $sjtool" < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/plot_junction_stats.R 

rm $inp $cls
