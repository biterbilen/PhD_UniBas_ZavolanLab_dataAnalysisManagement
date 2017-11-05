#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file

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

bti_dir=bti_dir
db_name=db_name
thread_num=thread_num
guide_gtff=guide_gtff
indir_tag="indir_tag"

outdir=output_directory
pushd $outdir;

of=Cuffmerge_$indir_tag
if [ -e $of/merged.gtf ]; then echo $of/merged.gtf exists; exit; fi #skip if folder exists

mkdir -p $of; pushd $of

inp=_assemblies.txt
ls ../Cufflinks_*$indir_tag*/*gtf > $inp

cmd="cuffmerge -g $guide_gtff -s $bti_dir/$db_name.fa -p $thread_num --keep-tmp $inp"

echo $cmd
$cmd

mv merged_asm/* .; rm -rf merged_asm; #TODO

echo DONE
