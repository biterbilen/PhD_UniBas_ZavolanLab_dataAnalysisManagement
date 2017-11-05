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
nhood=nhood
t2cCut=t2c_cut

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/run_call_plot_profile_error.sh  /import/bc2/soft/app/matlab/current/Linux/ $outdir $nhood $t2cCut

