#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@mpi.PE8
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-98
# $ -t 1-1

chrs=(chr14 chr20 chr6_random chr17_random chr17 chr21 chr19 chr21_random chr15 chr7_random chr7 chr13 chr16 chr2 chr6 chr5 chr19_random chr10 chr8 chr1 chr12 chrX chr11 chr22 chr4_random chr9 chrY chr3 chrM chr18 chr1_random chr4 chr3_random chr15_random chr2_random chr22_random chr5_random chrX_random chr13_random chr10_random chr8_random chr9_random chr16_random chr11_random chr22_h2_hap1 chr6_cox_hap1 chr18_random chr5_h2_hap1 chr6_qbl_hap2);
chrids=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49);

cd output_directory

strand="+";
if [ $SGE_TASK_ID -lt 50 ]; then strand="-"; fi
chridi=$(($SGE_TASK_ID % 49));
chrid=${chrids[$chridi]};
chr=${chrs[$chridi]};
outf=$chr$strand;

#Get densities
~bilebi00/CDS_CLIP_ElMMo/scripts/generateGenomeAlignmets.pl sample_id $chrid $strand T2C mRNA > $outf;
#~bilebi00/CDS_CLIP_ElMMo/scripts/generateGenomeAlignmets.pl 216 $chrid $strand T2C mRNA > $outf;

~bilebi00/CDS_CLIP_ElMMo/scripts/generateCumulativeScores.pl $outf > $outf-density

sort $outf-density -n -k 3 -r | awk -v chr=$chr -v strand=$strand 'BEGIN{OFS="\t";} { print chr, strand, $0; }' > densitySorted_$outf
rm $outf $outf-density
