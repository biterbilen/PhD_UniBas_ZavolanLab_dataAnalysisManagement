#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file_pairwise_distance_distribution
#$ -t 1-lib_pair_count


outdir=output_directory #Crosstalks
mkdir -p $outdir; cd $outdir

dir=input_directory #/import/bc2/home/zavolan/GROUP/RNA-BP_CLIP_Project/Methods/MRNAExpression
libs=(libs) #(TIA1_15_CR8 TIA1_15_CR9 TIA1_40_CR10 TIA1_40_CR11 hnRNPC_Clip9 hnRNPC_Clip10 Ago2_365nm_CLIP31 Ago2Nu_Clip11);
maxsites=max_sites #10000
rand=random_times #100
scrdir=scr_directory #~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods
trf=rep_trx_file #/import/bc2/home/zavolan/GROUP/RNA-BP_CLIP_Project/Methods/representative_transcripts
libCount=lib_count

i=`echo "$(( $(($SGE_TASK_ID-1)) % $libCount)) " | bc`;
j=`echo "$(( $(($SGE_TASK_ID-1)) / $libCount)) " | bc`;

if [ $i -le $j ]; then
	exit;
fi

pat1=${libs[$i]};
pat2=${libs[$j]};

echo $SGE_TASK_ID $i $j $pat1 $pat2

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/wrapper_intersite_distance.pl $dir $pat1 $pat2 $maxsites $rand $scrdir $trf
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/wrapper_site_summary.pl $dir $pat1 $pat2 $scrdir $maxsites $pat1.out $pat2.out
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_plot_distances.sh  /import/bc2/soft/app/matlab/current/Linux/ $pat1-$pat2.$maxsites.real $pat1-$pat2.$maxsites.random 

