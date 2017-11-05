#!/bin/bash
#$ -S /bin/bash
# $ -l sjpn=1
# $ -q fs_very_long@@qc_nehalem
#$ -q fs_long@@high_mem_node
#$ -P project_zavolan
#$ -l h_vmem=60G
#$ -j y
#$ -cwd
#$ -o LOG._run_hmmSplicer.chr.rep_A

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
#db_name=hg18_chr19
db_name=hg18 #hg_TR #TODO set here
tag=A
outdir=rep_$tag.$db_name

inp=inp.fastq
thread_num=48

mkdir -p $outdir; cd $outdir;

less ~bilebi00/CLIP_samples/mRNASeq_No_4SU_No_XL_rep_$tag/PROTS/mRNASeq_No_4SU_No_XL_rep_$tag.bwa_oligoAln | perl -e 'while(<>) {($id)=split;$chr=<>;chomp $chr;$_=<>;($x,$y,$z,$str)=split; $seq=<>; chomp $seq; $_=<>; $_=<>; $_=<>; $seq=~tr/ACGT/TGCA/ if ($str eq "-"); $seq =~ s/\-//g; print ">$id\t$seq\n" if (length($seq) == 37 and $chr =~ /chr/); } ' > $inp

less ~bilebi00/_SHORT_READ_ALIGNMENT/data/rep_${tag}_unmapped.seq | perl -e ' $_=<>;  while(<>) { @t=split; print ">sequd_$t[0]\t$t[1]\n" if (length($seq) == 37); }' >> $inp

less $inp | sort -u | perl -e 'while(<>){@t=split; print $t[0], "\n", $t[1], "\n"; }' > $inp.tmp
fq_all2std.pl fa2std $inp.tmp > $inp
rm $inp.tmp

echo python ~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/hmmSplicer-0.9.5/runHMM.py -o hmmSplicer -i $inp -g $bti_dir/$db_name.fa -b $bti_dir/$db_name -j 50 -k 500000 -d True -r True -p $thread_num -x False
python ~bilebi00/_DOWNLOAD/_FOR_ALIGNMENT/hmmSplicer-0.9.5/runHMM.py -o hmmSplicer -i $inp -g $bti_dir/$db_name.fa -b $bti_dir/$db_name -j 50 -k 500000 -d True -r True -p $thread_num -x False


