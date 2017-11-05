#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -t 1-runtimes
#$ -o log_file_$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1"
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

PATH=$HOME/bin:$PATH

samples=samples_file
db_name=db_name
bti_dir=bti_dir
thread_num=thread_num
outdir=output_directory
segment_len=segment_len
segment_mismatches=segment_mis

pushd $outdir;

inp=`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 1`
of=TopHat_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`

if [ -e $of/accepted_hits.bam ]; then echo $of/accepted_hits.bam exists; exit; fi #skip if folder exists

#cmd="tophat --segment-length 18 --segment-mismatches 1 -p $thread_num -o $of $bti_dir/$db_name $inp"
cmd="tophat --segment-length $segment_len --segment-mismatches $segment_mismatches -p $thread_num -o $of $bti_dir/$db_name $inp"
echo $cmd
$cmd
echo DONE

