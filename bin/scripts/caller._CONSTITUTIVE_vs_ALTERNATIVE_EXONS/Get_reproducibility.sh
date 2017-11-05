#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
# $ -q fs_long@@qc_nehalem
#$ -q fs_short@@qc_nehalem
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

pwd=`pwd`;
outdir=output_directory
mkdir -p $outdir; cd $outdir;

folders="folders";
outfile=outfile
maxExt=max_extension
nmer=central_nmer

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/calc_positional_reproducibility.pl $pwd $folders > $outfile 

less $outfile | perl -e 'my $me=shift; while(<>) { @t=split; $s=$t[1];$t[2]-=$me; $t[3]+=$me; splice(@t,1,1); splice(@t,3,0,$s); print join("\t", @t), "\n"; }' $maxExt | ~/PUM/bperlscripts/extractSeqsFromCoords.pl > $outfile.fa 

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/count_nc_content.pl $nmer $outfile.fa 8 > $outfile.rep

rm $outfile.fa $outfile 

exit;

