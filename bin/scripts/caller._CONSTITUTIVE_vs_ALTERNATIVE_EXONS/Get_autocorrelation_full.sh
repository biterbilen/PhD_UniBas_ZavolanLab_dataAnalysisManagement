#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -l sjpn=1
#$ -l mem_total=30G
#$ -q fs_short@@qc_nehalem
#$ -l h_vmem=50G
# $ -l h_vmem=1M
# $ -notify
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-98

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

chrs=(chr14 chr20 chr6_random chr17_random chr17 chr21 chr19 chr21_random chr15 chr7_random chr7 chr13 chr16 chr2 chr6 chr5 chr19_random chr10 chr8 chr1 chr12 chrX chr11 chr22 chr4_random chr9 chrY chr3 chrM chr18 chr1_random chr4 chr3_random chr15_random chr2_random chr22_random chr5_random chrX_random chr13_random chr10_random chr8_random chr9_random chr16_random chr11_random chr22_h2_hap1 chr6_cox_hap1 chr18_random chr5_h2_hap1 chr6_qbl_hap2);
chrids=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49);

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=output_directory

#set these from sed script
dist=distance
t2cCut=t2c_cut
tag=tag;

mkdir -p $outdir; cd $outdir;

str="+";
if [ $SGE_TASK_ID -lt 50 ]; then str="-"; fi
chridi=$(($SGE_TASK_ID % 49));
chrid=${chrids[$chridi]};
chr=${chrs[$chridi]};

inp=$chr$str.$tag.distance$dist.inp
cls=$inp.cls
#select top T2Cs
less sorted_T2C_densities | grep -P "$chr\s\\$str" | perl -e 'my $dist=shift; my $t2c_cut=shift; while(<>){chomp; ($ss,@t)=split; last if ($t[$#t] < $t2c_cut); $d=0; $s=1; $s=-1 if($t[4] eq "-"); print join("\t","${ss}_0",@t,$d),"\n"; for (1..$dist) { $t[2]-=$s; $t[3]-=$s; print join("\t","SS_$_",@t,$d),"\n"; }}' $dist $t2cCut > $inp

#cluster them
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $cls

#extract T2C based on nucleotide shift of distance
tagg="$chr$str${tag}.distance";
rm -rf $tagg*.corr;
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/parse_cls_for_autocorr_shifts.pl $cls 1 $dist $tagg

rm $inp $cls


