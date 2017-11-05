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

outdir=output_directory
pushd $outdir;
pwd=`pwd`;

soutdir=Reproducibility
mkdir -p $soutdir; pushd $soutdir;

folders="folders";
protein=protein
outfile=outfile
maxExt=max_extension
nmer=central_nmer

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/calc_positional_reproducibility.pl $pwd $folders > $outfile 
db=~bilebi00/DATA/hg19

less $outfile | perl -e 'my $me=shift; while(<>) { @t=split; $s=$t[1];$t[2]-=$me; $t[3]+=$me; splice(@t,1,1); splice(@t,3,0,$s); print join("\t", @t), "\n"; }' $maxExt | ~/PUM/bperlscripts/extractSeqsFromCoords.pl $db > $outfile.fa 

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/count_nc_content.pl $nmer $outfile.fa 8 1 > $outfile.rep

rm $outfile.fa $outfile 

popd
ln -s $soutdir/$outfile.rep $protein

echo DONE
exit;

