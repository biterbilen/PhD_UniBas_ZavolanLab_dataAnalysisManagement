#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_spliceArray_cuffdiffSplicing.sh 

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
export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=Scheiffele/array_comparison_cerebellum; 

mkdir -p $outdir; pushd $outdir;
arrayf=~bilebi00/_SHORT_READ_ALIGNMENT/data/Splice_arrays/Harald_31_03_2011_mouse_cerebellum/splicing-change-Harald_31_03_2011_mouse_cerebellum-High_Low_with_duplicates.txt
seqf=~bilebi00/_SHORT_READ_ALIGNMENT/Scheiffele/Cuffdiff_cerebellum/splicing.diff


~bilebi00/bin/R --no-save --args $arrayf < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/fill_table_empty_fields.R > /dev/null

~bilebi00/_SHORT_READ_ALIGNMENT/scripts/aggregate.pl sqrtJS 2 9 $seqf > sqrtJS_per_gene.stat
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/aggregate.pl absIrank 9 24  $arrayf.inp > absIrank_per_gene.stat

~bilebi00/_PAPD5/scripts/innerJoin.pl sqrtJS_per_gene.stat absIrank_per_gene.stat 1 1 '1-4,6-8' > array-cuffdiff_splicing_per_gene.stat

~bilebi00/bin/R --no-save --args array-cuffdiff_splicing_per_gene.stat < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/spliceArray_cuffdiffSplicing.R > /dev/null

less $seqf | grep "yes$" | cut -f 3 | perl -e 'while(<>){@t=split/,/; print join("\n", @t);}' | sort -u  > sqrtJS.significant_genes

cut1=1
less $arrayf.inp | awk -F "\t" '$25 != "NA" && $25 > cut && NR>1 { print $10;}' cut=$cut1 | perl -e 'while(<>){next if ($_ =~ /^$/); @t=split/,/; print join("\n", @t);}' | sort -u > Irank.enriched_genes.$cut1
cut2=1
less $arrayf.inp | awk -F "\t" ' $25 != "NA" && $25 < -cut && NR>1 { print $10;}' cut=$cut2 | perl -e 'while(<>){next if ($_ =~ /^$/); @t=split/,/; print join("\n", @t);}' | sort -u > Irank.depleted_genes.$cut2

cat Irank.enriched_genes.$cut1 Irank.depleted_genes.$cut2 | sort -u > Irank.significant_genes.$cut1.$cut2

~bilebi00/_PAPD5/scripts/innerJoin.pl Irank.significant_genes.$cut1.$cut2 sqrtJS.significant_genes 1 1 '1' > Irank.$cut1.$cut2.sqrtSJ.significant_genes.intersection

