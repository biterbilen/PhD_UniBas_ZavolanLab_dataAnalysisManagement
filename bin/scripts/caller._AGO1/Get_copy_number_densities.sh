#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -l s_vmem=50G
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

outdir=output_directory
sampleId=sample_id
scalenorm=scale_norm #0
annot=annot  #ALL
excAnnots="exc_annots"  

chrs=(chr14 chr20 chr6_random chr17_random chr17 chr21 chr19 chr21_random chr15 chr7_random chr7 chr13 chr16 chr2 chr6 chr5 chr19_random chr10 chr8 chr1 chr12 chrX chr11 chr22 chr4_random chr9 chrY chr3 chrM chr18 chr1_random chr4 chr3_random chr15_random chr2_random chr22_random chr5_random chrX_random chr13_random chr10_random chr8_random chr9_random chr16_random chr11_random chr22_h2_hap1 chr6_cox_hap1 chr18_random chr5_h2_hap1 chr6_qbl_hap2);
chrids=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49);

cd $outdir

strand="+";
if [ $SGE_TASK_ID -lt 50 ]; then strand="-"; fi
chridi=$(($SGE_TASK_ID % 49));
chrid=${chrids[$chridi]};
chr=${chrs[$chridi]};
outf=$chr$strand;

if [ $strand == "-" ]; then sleep 40m; fi

#Get densities
odir=$chr$strand;
~bilebi00/GEORGES_POLYA/scripts/copy_allFeaturesToWiggle.pl $sampleId $chrid $strand $odir $annot $scalenorm "$excAnnots";
if [ -s $odir/copies ]; then
	#corrects for 0 index
	#TODO ask chris if he uses the same code for T2Cs I subtracted 1 in that case
	less $odir/copies | perl -e 'while(<>){ @t=split; next if($t[3]==0); $t[1]+=1; shift @t; print join("\t", @t),"\n"; }' > $odir.copies
fi

rm -rf $odir;

