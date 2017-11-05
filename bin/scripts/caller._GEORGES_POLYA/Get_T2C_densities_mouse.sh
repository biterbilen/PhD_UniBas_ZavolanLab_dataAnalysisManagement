#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@mpi.PE8
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file
#$ -t 1-70
# $ -t 1-1

chrs=(chr1_random chrY_random chr18 chrX chr9_random chrUn_random chr16 chrY chr7 chr17 chrX_random chr6 chr5 chr15 chr12 chr9 chr1 chr19 chr13 chr10 chrM chr13_random chr11 chr8 chr2 chr14 chr4 chr3 chr8_random chr4_random chr7_random chr5_random chr17_random chr3_random chr16_random);                                                                                      
chrids=(50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84);

cd output_directory

strand="+";
if [ $SGE_TASK_ID -lt 36 ]; then strand="-"; fi
chridi=$(($SGE_TASK_ID % 35));
chrid=${chrids[$chridi]};
chr=${chrs[$chridi]};
outf=$chr$strand;

#Get densities
~bilebi00/CDS_CLIP_ElMMo/scripts/generateGenomeAlignmets.pl sample_id $chrid $strand T2C mRNA > $outf;
#~bilebi00/CDS_CLIP_ElMMo/scripts/generateGenomeAlignmets.pl 216 $chrid $strand T2C mRNA > $outf;

~bilebi00/CDS_CLIP_ElMMo/scripts/generateCumulativeScores.pl $outf > $outf-density

sort $outf-density -n -k 3 -r | awk -v chr=$chr -v strand=$strand 'BEGIN{OFS="\t";} { print chr, strand, $0; }' > densitySorted_$outf
rm $outf $outf-density
