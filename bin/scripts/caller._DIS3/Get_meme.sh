#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
# $ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-runtimes

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

export PATH=$HOME/bin:/import/bc2/soft/bin/hmmer:$PATH

outdir=output_directory

nmotifs=nmotifs
minw=minw
maxw=maxw
topN=topN
shuffleN=shuffleN
mod=mod

minsites=$((topN / 4));
maxsites=$topN;
if [ $mod == "anr" ]; then maxsites=$((topN * 2)); fi

meme=~bilebi00/bin/meme
shuffle=/import/bc2/soft/bin/hmmer/shuffle
tag=tag

pushd $outdir/MEME;

if [ $SGE_TASK_ID -gt $shuffleN ]; then
	odir=$tag;
	fa=$odir/seqs.fa
else
	odir=$tag.$SGE_TASK_ID;
	fa=$odir/seqs.fa
	mkdir -p $odir;
	less $tag/shuffled.fa | perl -e '$sorder=shift; $ts=shift; $i=1; while(<>){ $d=$_; $s=<>; print "$d$s" if ($i == $sorder); $i=0 if ($i==$ts); $i++; }' $SGE_TASK_ID $shuffleN > $fa
fi

$meme $fa -sf $fa -dna -mod $mod -nmotifs 1 -minsites $minsites -maxsites $maxsites -minw $minw -maxw $maxw -oc $odir -nostatus > /dev/null 
