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
annot=annot
libcut=libcut
shuffleN=shuffleN

shuffle=/import/bc2/soft/bin/hmmer/shuffle
fa=$tag.fa

mkdir -p $outdir/MEME/$tag
pushd $outdir/MEME/$tag;

less ../../$prot.annot.sorted | perl -e '$annot=shift; while(<>){ @t=split; next if ($t[0] ne $annot);  $c=($t[5] =~ s/[acgt]/X/g); if ($c<10) { print $_; }}' $annot > $annot.sorted.max9repeat

#input
less $annot.sorted.max9repeat | perl -e '$N=shift; $lc=shift; $seqc=0; while(<>){ last if ($seqc == $N);  @t=split; if ($t[9]>=$lc) { print ">seq$seqc\n$t[5]\n"; $seqc++; }}' $topN $libcut > seqs.fa

$shuffle -d -n $shuffleN seqs.fa > shuffled.fa

echo DONE

exit;

