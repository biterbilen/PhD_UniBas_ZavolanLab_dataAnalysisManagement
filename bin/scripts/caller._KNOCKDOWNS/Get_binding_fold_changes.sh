#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -cwd
#$ -j y
#$ -o log_file
#$ -t 1-jobs_count

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

export PATH=~bilebi00/bin:$PATH

outdir=output_directory
rnaSeqNormTable=rna_seq_norm_table
fgroups=(f_groups)
bgroups=(b_groups)
cgroups="c_groups"
cfolder=c_folder
trxi=trx_intersect
nobt=non_over_bin_type

mkdir -p $outdir; cd $outdir

f=${fgroups[$(($(($SGE_TASK_ID-1))/${#bgroups[@]}))]};
b=${bgroups[$(($(($SGE_TASK_ID-1))%${#bgroups[@]}))]};
if [ $f != $b ]; then 
	echo Getting Fold changes for foreground group $f and background group $b;
	for c in $cgroups; do
		~bilebi00/_KNOCKDOWNS/scripts/getFoldChange_forBinnedIdMatch.pl $cfolder/$c/windows.foreground.allf.fa $c $f $b $rnaSeqNormTable $trxi $nobt;
	done
	~bilebi00/_KNOCKDOWNS/scripts/run_call_wrapper_barweb.sh /import/bc2/soft/app/matlab/current/Linux/ $f
fi
