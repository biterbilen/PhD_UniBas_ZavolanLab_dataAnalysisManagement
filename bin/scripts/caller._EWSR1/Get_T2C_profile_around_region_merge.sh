#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o log_file

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
tag=tag
prot_libs_libcuts_f=protein_libs_libcuts_file
region_type=region_type
regions_f=regions_f

pushd $outdir;

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/prep_profiles_forplotting.pl $tag.profile $tag.profile.stat "./" $tag .stat;

less $prot_libs_libcuts_f | cut -f 1 | while read i; do 
	echo $i; 
	R --no-save --args $tag.profile $i < ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/plot_profile_of_event.R; 
done; 

grep -P "\t$region_type" $regions_f | cut -f 1 | while read i; do 
	echo $i; R --no-save --args $tag.profile $i < ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/plot_profile_around_region.R; 
done; 

R --no-save --args $tag.profile.stat < ~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/plot_regional_pereference.R;
	
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$tag.pdf $tag.profile_*pdf $tag.profile.stat.pdf;

echo DONE
