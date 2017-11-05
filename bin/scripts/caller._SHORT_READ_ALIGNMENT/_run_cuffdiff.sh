#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._run_new_Cuffdiff.chr19.rep_A

PATH=$HOME/bin:$PATH

bti_dir=~/DATA/Bowtie_DB_INDEX
db_name=hg18_chr19 #hg_TR #TODO set here
tag=A
outdir=rep_$tag.$db_name

thread_num=8

mkdir -p $outdir; cd $outdir;

inp=assemblies.txt
#gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/refseq_genes.gtf
gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf
out=cuffdiff

cmd="cuffcompare -R -r $gtff -s $bti_dir/$db_name.fa -i $inp -o $out"

echo $cmd
$cmd

