#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-3
#$ -o LOG._Get_T1nuclease_signature.$TASK_ID

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
#SGE_TASK_ID=$1
#TODO set
tag=EWSR1_FlpIn
tag=EWSR1

cell=""; #HEK293
indir=~bilebi00/_EWSR1/data/clipz_CLIP_tags_mRNA_repeat_none/
db=hg19
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list

outdir=Analysis/CLIP/T1signature
mkdir -p $outdir; pushd $outdir;

ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort`)

if [ $SGE_TASK_ID -le ${#ids[@]} ]; then

	id=${ids[$((SGE_TASK_ID-1))]}
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	if [ ! -e $indir/$id.bed.gz ]; then
		mkdir -p $indir
		#TODO ths part does not work from the clusters
		rsync -t mirz@web08:~/BITER/CLIP/scratch/bed_files_hg19/$id.bed.gz $indir/.
	else
		echo $indir/$id.bed.gz exists
	fi
	#get start base signature
	#T1
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	zless $indir/$id.bed.gz | \
		bedtools flank -s -l 1 -r 0 -i stdin -g ~/aux/human.$db.genome | \
		bedtools nuc -s "-fi" ~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa -bed stdin | \
		awk 'BEGIN{OFS="\t"}{ if ($1~/chr/) { print idn,$9,$10,$11,$12}}' idn="$idn" | \
		bedtools groupby  -g 1 -c 2,3,4,5 -o sum,sum,sum,sum > $tag.signature.$id
	#length signatures
	zless $indir/$id.bed.gz | \
		awk 'BEGIN{OFS="\t"}{ print idn,$3-$2; }' idn="$idn" | \
		bedtools groupby  -g 1 -c 2,2,2 -o mean,stdev,count -i stdin > $tag.signature_length.$id
	touch $tag.runfinished$SGE_TASK_ID
fi


if [ $SGE_TASK_ID == 1 ]; then
	#wait
	while [ "`ls $tag.runfinished* | wc -l`" -lt ${#ids[@]} ]; do
		echo waiting for other runs
		sleep 2m;
	done

	#merge
	#T1
	statf=$tag.T1nuclease_signature.stats
	echo -e "tag\tA\tC\tG\tT" > $statf
	cat $tag.signature\.* >> $statf
	#length
	statf=$tag.length_signature.stats
	echo -e "tag\tmean\tstdev\tN" > $statf
	cat $tag.signature_length\.* >> $statf

	#plot
	#TODO

	#clean
	rm $tag.signature* $tag.runfinished*
fi

echo DONE
