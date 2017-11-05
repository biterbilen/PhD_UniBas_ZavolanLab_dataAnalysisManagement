#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file_$TASK_ID
#$ -t 1-2

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
outdir=output_directory
rmsk_gtff=rmsk
rnaGene_gtff="rnaGene"
indir_tag="indir_tag"

pushd $outdir;

assemblies=(combined trusted);
assembly=${assemblies[$(($SGE_TASK_ID-1))]}
if=Cuffcompare_$indir_tag/Cuffcompare_$indir_tag.$assembly.gtf
#if=Cuffcompare_$indir_tag/Cuffcompare_$indir_tag.combined.gtf
#if=Cuffmerge_$indir_tag/merged.gtf
od=Cuffdiff_${assembly}_$indir_tag
if [ -e $od ]; then echo $od exists; exit; fi #skip if folder exists

#This part is redundant
maskfile=maskfile.$RANDOM
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rich' $rmsk_gtff > $maskfile
if [ $rnaGene_gtff != rnaGene ]; then
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rRNA,Mitochondrial' $rnaGene_gtff >> $maskfile
fi

samf="";
labels="";
for i in TopHat*$indir_tag*/*sam; do
	samf="$samf $i"
	l=`perl -e '$_=shift; $_=~/TopHat_(\S+)\//; print $1;' $i`;
	if [ "$labels" != "" ]; then
		labels="$labels,$l";
	else
		labels="$l";
	fi
done

cmd="cuffdiff -L $labels --no-update-check -q -M $maskfile -b $bti_dir/$db_name.fa -N -u -p $thread_num -o $od $if $samf"

echo $cmd
$cmd
mv $maskfile $od/_maskfile
echo DONE
