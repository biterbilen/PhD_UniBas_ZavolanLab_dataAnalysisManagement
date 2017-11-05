#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -j y
#$ -cwd
#$ -o log_file_trusted_sites
#$ -t 1-lib_rep_count

outdir=output_directory #Crosstalks
mkdir -p $outdir; cd $outdir

dir=input_directory #/import/bc2/home/zavolan/GROUP/RNA-BP_CLIP_Project/Methods/MRNAExpression
maxsites=max_sites #10000
libs=(libs) #("TIA1_15_CR8 TIA1_15_CR9" "TIA1_40_CR10 TIA1_40_CR11" "hnRNPC_Clip9 hnRNPC_Clip10" "Ago2_365nm_CLIP31 Ago2Nu_Clip11")

reps="${libs[$(($SGE_TASK_ID-1))]}";

ftag=`echo $reps | sed -e 's/ /_/g'`;
folders=($reps);
minLibCov=${#folders[@]};
rm -rf $ftag.inp
for pat in $reps; do
	less $dir/$pat/windows.foreground.allf | perl -e 'while(<>) {@s = split(/\|/, $_); print $s[3]/$s[4], "\t$_";}' | sort -gr -k 1 | head -1000 | awk '{print $2}' | perl -e 'my $label=shift;while(<>) {$_ =~ s/\s//g; @s = split(/\|/, $_); print "$label-$_\t$s[0]\t+\t$s[1]\t$s[2]\t",$s[3]/$s[4], "\n";}' $pat >> $ftag.inp 
done	

~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $ftag.inp > $ftag.cls.tmp
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/Methods/extract_trusted_sites.pl $ftag.cls.tmp $minLibCov > $ftag.trusted 

rm $ftag.inp $ftag.cls.tmp


