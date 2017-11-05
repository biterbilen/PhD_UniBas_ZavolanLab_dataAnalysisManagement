#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@high_mem_node
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

#set these from sed script
outdir=output_directory
dist=distance
tag=tag;
t2cCut=t2c_cut

cat $outdir/chr*$tag.distance$dist.corr > $outdir/$tag.distance$dist.t2cCut${t2cCut}.corr; 
rm $outdir/chr*$tag.distance*.corr

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_autocorrelation_full.sh  /import/bc2/soft/app/matlab/current/Linux/ $outdir $dist $tag unnormalized_autocorrelation_full_merged
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_autocorrelation_full_pairwise.sh  /import/bc2/soft/app/matlab/current/Linux/ $outdir $dist $tag pairwise_autocorrelation_full_merged

