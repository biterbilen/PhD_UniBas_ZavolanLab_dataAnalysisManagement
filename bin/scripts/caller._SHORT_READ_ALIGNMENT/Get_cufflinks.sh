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
guide_gtff=guide_gtff
rmsk_gtff=rmsk
rnaGene_gtff="rnaGene"

pushd $outdir;

id=TopHat_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`
od=Cufflinks_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`

if [ -e $od/transcripts.gtf ]; then echo $od/transcripts.gtf exists; exit; fi #skip if folder exists

maskfile=maskfile.$RANDOM
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rich' $rmsk_gtff > $maskfile
if [ $rnaGene_gtff != rnaGene ]; then
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rRNA,Mitochondrial' $rnaGene_gtff >> $maskfile
fi

#-j --pre-mrna-fraction <0.0-1.0>
m=8
cmd="cufflinks --no-update-check -q -M $maskfile -b $bti_dir/$db_name.fa -N -u -p $thread_num -o $od -g $guide_gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 $id/accepted_hits.bam"

echo $cmd
$cmd
mv $maskfile $od/_maskfile
echo DONE
