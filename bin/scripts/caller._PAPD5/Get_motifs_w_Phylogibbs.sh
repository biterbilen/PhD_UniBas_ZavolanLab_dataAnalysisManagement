#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-40

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

export PERL5LIB=$HOME/lib:$PERL5LIB
phylogibbs=/import/bc2/home/nimwegen/GROUP/PhyloGibbs/DistributionApr9_2006/src/phylogibbs

outdir=output_directory
tag=tag
infile=in_file
t2ccut=t2c_cut
libcut=lib_cut
reglen=region_len
numSitesFrac=num_sites_frac

mkdir -p $outdir/$tag; cd $outdir/$tag;
ftag=ftag
file=$ftag.fa
less ../$infile | perl -e '$t2ccut=shift; $libcut=shift; $reglen=shift; while(<>){ @t=split; if ($t[5]>=$t2ccut and $t[6]>=$libcut) { my $i = int((length($t[2]) - $reglen) / 2 ); print ">$_".uc(substr($t[2], $i, $reglen))."\n"; } }' $t2ccut $libcut $reglen > $file;

track_file=$ftag.track;
out_file=$ftag.out;
wm_file=$ftag.wm;
num_seq=`wc -l $file | awk '{ print $1/2;}'`;
run=`echo "($SGE_TASK_ID-1)%10+1" |bc`; 
window_size=`echo "($SGE_TASK_ID-1)/10+5" |bc`;
num_motifs=1

#find motifs
out_dir=$ftag-run$run-$window_size-$num_motifs;	
rm -rf $out_dir 2> /dev/null;
mkdir -p $out_dir; 
cd $out_dir;
#should be seen in at least numSitesFrac of the sites
num_sites=`echo "$num_seq * $numSitesFrac * $num_motifs" | bc | awk -F . '{ print $1;}'`;
#N=1 markov order 1; single base frequencies in the background file
#N=0 markov order 0; 0.25 for each base
# $phylogibbs -D 0 -F ../$file -f ../$file -m $window_size -N 1 -o $out_file -S 500 -t $track_file -y $num_sites -z $num_motifs -r
$phylogibbs -D 0 -F ../$file -f ../$file -m $window_size -N 0 -o $out_file -S 500 -t $track_file -y $num_sites -z $num_motifs -r
~bilebi00/PUM/bperlscripts/make_logo_from_phylogibbs_file.pl $track_file .ws$window_size.r$run > $wm_file
# ~bilebi00/PUM/bperlscripts/make_logos.pl $wm_file 2> /dev/null

