#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
# $ -q fs_long@@qc_nehalem
#$ -l sjpn=1
# $ -q fs_long@@high_mem_node
#$ -P project_zavolan
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
ffile=ffile
bfile=bfile
minValue=min_value

mkdir -p $outdir; cd $outdir;
#The name should include DB
filetag=`perl -e '{$s=shift; $s=~/DB\w*_(\d+)/; print $1;}' $ffile`_`perl -e '{$s=shift; $s=~/DB\w*_(\d+)/; print $1;}' $bfile`;
echo $filetag

#~bilebi00/_PAPD5/scripts/tablize.pl 5 $filetag.inp 0 $ffile $bfile 
~bilebi00/_PAPD5/scripts/trim_min_table_values.pl $filetag.inp $minValue > $filetag.inp.trimmed
~bilebi00/_PAPD5/scripts/zScore.pl $filetag.inp.trimmed "1 0" 0 > $filetag.zscore 

rm $filetag.inp*;

