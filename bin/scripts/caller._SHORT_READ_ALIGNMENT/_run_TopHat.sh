#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_very@@qc_nehalem
# $ -q fs_long@@high_mem_node
#$ -P project_zavolan
# $ -l h_vmem=60G
#$ -j y
#$ -cwd
#$ -o LOG._run_topHat.chr19.rep_A

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

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
chr=chr19
db_name=hg18_$chr #hg_TR #TODO set here
tag=B
outdir=rep_$tag.$db_name

thread_num=48
inp=inp.fa

mkdir -p $outdir; cd $outdir;

less ~bilebi00/CLIP_samples/mRNASeq_No_4SU_No_XL_rep_$tag/PROTS/mRNASeq_No_4SU_No_XL_rep_$tag.bwa_oligoAln | perl -e '$l=shift; $c=shift; while(<>) {($id)=split;$chr=<>;chomp $chr;$_=<>;($x,$y,$z,$str)=split; $seq=<>; chomp $seq; $_=<>; $_=<>; $_=<>; if ($chr =~ /chr/) { $seq=~tr/ACGT/TGCA/ if ($str eq "-"); $seq =~ s/\-//g; print ">$id\t$seq\n" if (length($seq) == $l and ($c eq "" or $c eq $chr)); }} ' $len $chr > $inp
less ~bilebi00/_SHORT_READ_ALIGNMENT/data/rep_${tag}_unmapped.seq |perl -e ' $_=<>;  while(<>) { @t=split; print ">sequd_$t[0]\t$t[1]\n" if (length($t[1]) == 37); }' >> $inp

less $inp | sort -u | perl -e 'while(<>){@t=split; print $t[0], "\n", $t[1], "\n"; }' > $inp.tmp
mv $inp.tmp $inp

echo tophat --segment-length 18 --segment-mismatches 1 --keep-tmp -p $thread_num -o tophat $bti_dir/$db_name $inp
tophat --segment-length 18 --segment-mismatches 1 --keep-tmp -p $thread_num -o tophat $bti_dir/$db_name $inp


