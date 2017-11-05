#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-70

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

chrs=(chr1_random chrY_random chr18 chrX chr9_random chrUn_random chr16 chrY chr7 chr17 chrX_random chr6 chr5 chr15 chr12 chr9 chr1 chr19 chr13 chr10 chrM chr13_random chr11 chr8 chr2 chr14 chr4 chr3 chr8_random chr4_random chr7_random chr5_random chr17_random chr3_random chr16_random);
chrids=(50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84);

cd $outdir

strand="+";
if [ $SGE_TASK_ID -lt 36 ]; then strand="-"; fi
chridi=$(($SGE_TASK_ID % 35));
chrid=${chrids[$chridi]};
chr=${chrs[$chridi]};
outf=$chr$strand;

#Get densities
odir=$chr$strand;
~bilebi00/GEORGES_POLYA/scripts/copy_allFeaturesToWiggle.pl $sampleId $chrid $strand $odir mRNA;
if [ -s $odir/copies ]; then
	#corrects for 0 index
	#TODO ask chris if he uses the same code for T2Cs I subtracted 1 in that case
	less $odir/copies | perl -e 'while(<>){ @t=split; next if($t[3]==0); $t[1]+=1; shift @t; print join("\t", @t),"\n"; }' > $odir.copies
fi

rm -rf $odir;

