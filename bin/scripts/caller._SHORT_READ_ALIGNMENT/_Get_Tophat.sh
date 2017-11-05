#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
# $ -q fs_long@@high_mem_node
#$ -P project_zavolan
# $ -l h_vmem=60G
#$ -j y
#$ -cwd
# $ -t 1-12
# $ -t 1-2
#$ -t 5-10
# $ -o LOG._Get_Tophat_Guo_$TASK_ID
# $ -o LOG._Get_Tophat_Rep_chr19_$TASK_ID
#$ -o LOG._Get_Tophat_Zavolan_$TASK_ID

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

bti_dir=~bilebi00/DATA/Bowtie_DB_INDEX
thread_num=8

#TODO set
db_name=hg18 #hg_TR #TODO set here
#db_name=hg18_chr19 #hg_TR #TODO set here

#TODO set task number to
outdir=Guo #TODO set task number to
outdir=Rep_chr19
outdir=Zavolan

mkdir -p $outdir; cd $outdir;

tag=TopHat
#TODO param
samples=~bilebi00/_SHORT_READ_ALIGNMENT/data/Guo_et_al_2010_nature_GSE21992/samples
samples=~bilebi00/_SHORT_READ_ALIGNMENT/data/Rep_chr19_samples
samples=~bilebi00/_SHORT_READ_ALIGNMENT/data/${outdir}_samples

inp=`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1`

#TODO param
#1.
#chr=chr19
#len=37
#less ~bilebi00/CLIP_samples/mRNASeq_No_4SU_No_XL_$inp/PROTS/mRNASeq_No_4SU_No_XL_$inp.bwa_oligoAln | perl -e '$l=shift; $c=shift; while(<>) {($id)=split;$chr=<>;chomp $chr;$_=<>;($x,$y,$z,$str)=split; $seq=<>; chomp $seq; $_=<>; $_=<>; $_=<>; if ($chr =~ /chr/) { $seq=~tr/ACGT/TGCA/ if ($str eq "-"); $seq =~ s/\-//g; print ">$id\t$seq\n" if (length($seq) == $l and ($c eq "" or $c eq $chr)); }} ' $len $chr | sort -u | perl -e 'while(<>){@t=split; print $t[0], "\n", $t[1], "\n"; }' > $inp
#2.
#inp=~bilebi00/_SHORT_READ_ALIGNMENT/data/Guo_et_al_2010_nature_GSE21992/${inp}_filtered_sequence.txt
#3.
inp=~bilebi00/_SHORT_READ_ALIGNMENT/data/$outdir/$inp.fastq

of=${tag}_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`

cmd="tophat --segment-length 18 --segment-mismatches 1 --keep-tmp -p $thread_num -o $of $bti_dir/$db_name $inp"
echo $cmd
$cmd

rm -rf $of/tmp; #clean up space
