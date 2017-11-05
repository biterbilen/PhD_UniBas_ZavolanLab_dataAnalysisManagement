#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q mpi8_very_long
#$ -l sjpn=1
#$ -t 1-4
#$ -j y
#$ -cwd
#$ -o LOG._Get_cufflinks_suite.$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1";
}
trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH

#TODO set
project=ARE
project=DIS3L2
project=MethodValidation
project=EWSR1

themedir=/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
fraction=1 #sampling fraction for backward compatibility
indir=~bilebi00/_CLIP/data/clipz_UniqueGenomicAlignments

#ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Cufflinks" { print $2;}' project=$project`)
files=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Cufflinks" { print $13;}' project=$project`)

if [ $SGE_TASK_ID -gt ${#files[@]} ]; then
	exit;
fi

id=${ids[$((SGE_TASK_ID-1))]}
f=${files[$((SGE_TASK_ID-1))]}

echo Doing $id $f

outdir=Analysis/Project_$project/Cufflinks/
mkdir -p $outdir; pushd $outdir;

#tophat
BOWTIE2_INDEXES=~/DATA/Bowtie2_DB_INDEX/hg19
GTF=~/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_exon.gtf.gz

if [ $SGE_TASK_ID == 1 ]; then
	#select ncRNA and mRNA genes
	if [ ! -e refGene.gtf ]; then
		zless $GTF | awk -F "\t" 'BEGIN{OFS="\t"} { if ($3 == "mRNA" || $3 == "ncRNA") { $3="exon"; print }}' > refGene.gtf
		touch prep_GTF
	else
		echo refGene.gtf exists
	fi
fi

#wait GTF file prep
while [ ! -e prep_GTF ]; do
	echo waiting for prep_GTF
	sleep 1m
done

#spliced alignment
if [ ! -e out_tophat_$SGE_TASK_ID ]; then
	tophat --no-novel-juncs -a 5 -G refGene.gtf -p 8 -o out_tophat_$SGE_TASK_ID $BOWTIE2_INDEXES $f
	#tophat --segment-length 18 --segment-mismatches 1 -p $thread_num -o $of $bti_dir/$db_name $inp
else
	echo out_tophat_$SGE_TASK_ID exists
fi

#transcript assembly
if [ ! -e out_cufflinks_$SGE_TASK_ID ]; then
	cufflinks --no-update-check -u -b $BOWTIE2_INDEXES.fa -p 8 -o out_cufflinks_$SGE_TASK_ID out_tophat_$SGE_TASK_ID/accepted_hits.bam
else
	echo out_cufflinks_$SGE_TASK_ID exists
fi

if [ $SGE_TASK_ID == 1 ]; then

	#merge assemblies
	if [ ! -e merged_asm ]; then
		find . -name transcripts.gtf > assemblies.txt
		cuffmerge -g refGene.gtf -s $BOWTIE2_INDEXES.fa -p 8 assemblies.txt
	else
		echo merged_asm exists
	fi

	#differential expression
	cuffdiff -o out_cuffdiff -u -b $BOWTIE2_INDEXES.fa -p 8 -L siCTRL,siEWSR1 merged_asm/merged.gtf \
		./out_tophat_1/accepted_hits.bam,./out_tophat_3/accepted_hits.bam \
		./out_tophat_2/accepted_hits.bam,./out_tophat_4/accepted_hits.bam
fi

echo DONE

