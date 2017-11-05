#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file_interesting_crosstalks
#$ -t 1-rep_lib_pair_count

outdir=output_directory #Crosstalks
mkdir -p $outdir; cd $outdir

libs=(libs) #("TIA1_15_CR8 TIA1_15_CR9" "TIA1_40_CR10 TIA1_40_CR11" "hnRNPC_Clip9 hnRNPC_Clip10" "Ago2_365nm_CLIP31 Ago2Nu_Clip11")
distance=xtalk_dist #10;
replibcount=rep_lib_count

i=`echo "$(( $(($SGE_TASK_ID-1)) % $replibcount)) " | bc`;
j=`echo "$(( $(($SGE_TASK_ID-1)) / $replibcount)) " | bc`;

if [ $i -le $j ]; then
	exit;
fi

rep1=`echo ${libs[$i]} | sed -e 's/ /_/g'`; 
rep2=`echo ${libs[$j]} | sed -e 's/ /_/g'`; 

echo $SGE_TASK_ID $i $j $rep1 $rep2

rand=$RANDOM
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $rep1.trusted > $rep1.cls.$rand
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/summarize_clusters.pl $rep1.cls.$rand $rep1 > $rep1.inxtalks.$rand
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $rep2.trusted > $rep2.cls.$rand
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/summarize_clusters.pl $rep2.cls.$rand $rep2 > $rep2.inxtalks.$rand

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/summarize_crosstalks.pl $rep1.inxtalks.$rand $rep2.inxtalks.$rand $distance > $rep1\_$rep2.xtalk 

rm -rf $rep1.cls.$rand $rep2.cls.$rand $rep1.inxtalks.$rand $rep2.inxtalks.$rand
