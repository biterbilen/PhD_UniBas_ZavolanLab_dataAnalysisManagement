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

db_name=db_name
bti_dir=bti_dir
thread_num=thread_num
guide_gtff=guide_gtff
indir_tag="indir_tag"

outdir=output_directory
pushd $outdir;

of=Cuffcompare_$indir_tag
if [ -e $of/*.combined.gtf ]; then echo $of/*.combined.gtf exists; exit; fi #skip if folder exists

mkdir -p $of; pushd $of

inp=_assemblies.txt
ls ../Cufflinks_*$indir_tag*/*gtf > $inp

cmd="cuffcompare -R -r $guide_gtff -s $bti_dir/$db_name.fa -i $inp -o $of"

echo $cmd
$cmd

#select j and = class transcripts which are assigned a TSS (are they called primary transcripts?)
minfpkm=1
minfraccut=0.6; #TODO set in the script 
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/select_trusted.pl $minfpkm $minfraccut *.combined.gtf *.tracking _assemblies.txt > Cuffcompare_$indir_tag.trusted.gtf
#less Cuffcompare_$indir_tag.combined.gtf |perl -e 'while (<>) { if ($_ =~ /class_code "[=j]"/ and $_=~/TSS/) { print $_;}  }' > Cuffcompare_$indir_tag.trusted.gtf 

echo DONE

