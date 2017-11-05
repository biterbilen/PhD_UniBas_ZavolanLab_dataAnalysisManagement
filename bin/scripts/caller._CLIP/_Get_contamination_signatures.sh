#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-3
#$ -o LOG._Get_contamination_signatures.$TASK_ID

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
#-----------------
#SGE_TASK_ID=$1 #for rsync
#TODO set
tag=EWSR1_FlpIn
tag=EWSR1

cell=""; #HEK293
indir=~bilebi00/_EWSR1/data/clipz_CLIP_tags_mRNA_repeat_none/
db=hg19
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list

outdir=Analysis/CLIP/contaminationSignature
mkdir -p $outdir; pushd $outdir;

ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort`)

if [ $SGE_TASK_ID -le ${#ids[@]} ]; then

	id=${ids[$((SGE_TASK_ID-1))]}
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	if [ ! -e $indir/$id.other.bed.gz ]; then
		mkdir -p $indir
		rsync -t mirz@web08:~/BITER/CLIP/scratch_unique_other/bed_files_hg19/$id.bed.gz $indir/$id.other.bed.gz
		exit;
	else
		echo $indir/$id.other.bed.gz exists
	fi
	#get start base signature
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	#TODO change directory so that it does not need xlinkEnrichment to be run first
	#TODO 2. add ~bilebi00/_CLIP/Analysis/UniqueRawData_mRNA_repeat_none  for mutation ratio and rewrite xlinkEnrichment
	#TODO for 2. bedplus with contamination copy count T2C count, check also wiggle files used in UniqueRawData_mRNA_repeat_none (it should be the same when the copies are extracted from db in clipz_CLIP_tags
	#check regions quickly before doing anything else; check nuc content; do motif search in top
	#TODO less EWSR1.signature.289 | awk -F "\t" 'NR>1{ r=$11/$10; if (r>b) { print r"\t"$0}}' b=0.00334973930141287037227 | sort -k 7,7gr | less
	ifile=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/EWSR1.$id.$id.T.xlink_count.gtf.gz
	ifile=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/xlinkEnrichment/EWSR1.$id.$id.T.xlink_mu_refmax.wo_singleton.bed.gz
	ifile=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/xlinkEnrichment/EWSR1.$id.$id.T.xlink_mu_refmax.bed.gz
	bedtools annotate -s -counts -i $ifile -files $indir/$id.bed.gz $indir/$id.other.bed.gz \
		> $tag.signature.$id
	#TODO ekle contamination count: $count>2 for filtering out for not losing low copy regions
	median=`less $tag.signature.$id | sort -k 5,5g | awk '$8>0{print "tag\t"($8/($7+$8))"\t"$1":"$2"-"$3"\t"$0}' | bedtools groupby -g 1 -c 2 -o median -i stdin | cut -f 2`
	echo $median
	less $tag.signature.$id | sort -k 5,5g | awk 'BEGIN{OFS="\t"}{ if (($8/($7+$8))<m){ x="*"} else {x="";} print x,($8/($7+$8)),$1":"$2"-"($3+50),$0;}' m=$median > $id.4ucsc

	touch $tag.runfinished$SGE_TASK_ID
fi


#TODO change
#if [ $SGE_TASK_ID == 1 ]; then
if [ $SGE_TASK_ID == 1 ]; then
	#wait
	while [ "`ls $tag.runfinished* | wc -l`" -lt ${#ids[@]} ]; do
		echo waiting for other runs
		sleep 1m;
	done

	#TODO
#	#merge
#	statf=$tag.contamination_signature.stats
#	echo -e "tag\tA\tC\tG\tT" > $statf
#	cat $tag.signature* >> $statf

	#plot
	#TODO

	#clean
	rm $tag.signature* $tag.runfinished*
fi

echo DONE
