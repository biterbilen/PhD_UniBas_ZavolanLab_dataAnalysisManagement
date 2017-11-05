#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l h_vmem=60G
#$ -j y
#$ -cwd
#$ -o LOG._run_SpliceMap.repA

echo Cannot run <50nt sequences
exit

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
db_name=hg18_chr19 #hg_TR #TODO set here
tag=A
outdir=splicemap_rep_$tag.$db_name

inp=inp.seq
thread_num=8

mkdir -p $outdir; cd $outdir;
cp ~bilebi00/_SHORT_READ_ALIGNMENT/data/SpliceMap.cfg run.cfg;

less ~bilebi00/CLIP_samples/mRNASeq_No_4SU_No_XL_rep_$tag/PROTS/mRNASeq_No_4SU_No_XL_rep_$tag.bwa_oligoAln | perl -e 'while(<>) {($id)=split;$chr=<>;chomp $chr;$_=<>;($x,$y,$z,$str)=split; $seq=<>; chomp $seq; $_=<>; $_=<>; $_=<>; if ($chr =~ /chr19/) { $seq=~tr/ACGT/TGCA/ if ($str eq "-"); $seq =~ s/\-//g; print ">$id\t$seq\n" if (length($seq) == 37); }} ' > $inp

less ~bilebi00/_SHORT_READ_ALIGNMENT/data/rep_${tag}_unmapped.seq |perl -e ' $_=<>;  while(<>) { @t=split; print ">sequd_$t[0]\t$t[1]\n" if (length($t[1]) == 37); }' >> $inp

less $inp | sort -u | perl -e 'while(<>){@t=split; print uc($t[1]), "\n"; }' > $inp.tmp
mv $inp.tmp $inp

~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/SpliceMap3352_example_linux-64/bin/runSpliceMap run.cfg

