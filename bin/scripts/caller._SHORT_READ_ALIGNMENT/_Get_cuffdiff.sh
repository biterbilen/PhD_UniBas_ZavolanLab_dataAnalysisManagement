#!/bin/bash
#$ -S /bin/bash
#$ -l sjpn=1
#$ -q fs_long@@qc_nehalem
# $ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -t 1-12
#$ -o LOG._Get_cufflinks_Guo_$TASK_ID

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

PATH=$HOME/bin:$PATH

bti_dir=~bilebi00/DATA/Bowtie_DB_INDEX
db_name=hg18 #hg_TR #TODO set here
thread_num=8
outdir=Guo #TODO set task number to

mkdir -p $outdir; cd $outdir;

#gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/refseq_genes.gtf
gtff=~bilebi00/_SHORT_READ_ALIGNMENT/data/GMAP_EXONS.gtf
maskfile=maskfile.$RANDOM
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rich' ~bilebi00/DATA/hg18_ucsc_tracks/rmsk.gtf.gz > $maskfile
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/unambiguate_gff_transcript_ids.pl 'rRNA,Mitochondrial' ~bilebi00/DATA/hg18_ucsc_tracks/rRNAGene.gtf.gz >> $maskfile

samples=~bilebi00/_SHORT_READ_ALIGNMENT/data/Guo_et_al_2010_nature_GSE21992/samples
id=TopHat_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`
od=Cufflinks_`less $samples | head -n $SGE_TASK_ID | tail -n 1 | cut -f 2`

#-j --pre-mrna-fraction <0.0-1.0>
m=8
cmd="cufflinks --no-update-check -q -M $maskfile -b $bti_dir/$db_name.fa -N -u -p $thread_num -o $od -g $gtff --min-frags-per-transfrag $m --min-intron-length 50 -I 500000 $id/accepted_hits.bam"

echo $cmd
$cmd
mv $maskfile $od/maskfile
echo DONE
