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

outdir=output_directory
prot=prot
tag=tag
topN=topN
shuffleN=shuffleN

shuffle=/import/bc2/soft/bin/hmmer/shuffle
fa=$tag.fa

pushd $outdir/MEME;

#FIXME will fail when submitted parallel
if [ ! -e ../$prot.annot.sorted ]; then 
	less -S ../$prot.annot | sort -k 9,9gr > ../$prot.annot.sorted
fi

if [ ! -e ../$prot.annot.sorted.max9repeat ]; then 
	less ../$prot.annot.sorted | perl -e 'while(<>){ @t=split; $c=($t[5] =~ s/[acgt]/X/g); if ($c<10) { print $_; }}' > ../$prot.annot.sorted.max9repeat
fi

#input
mkdir -p $tag
less ../$prot.annot.sorted.max9repeat | perl -e '$N=shift; $seqc=0; while(<>){ last if ($seqc == $N);  @t=split; if ($t[9]==2) { print ">seq$seqc\n$t[5]\n"; $seqc++; }}' $topN > $tag/seqs.fa

$shuffle -d -n $shuffleN $tag/seqs.fa > $tag/shuffled.fa

echo DONE

exit;

