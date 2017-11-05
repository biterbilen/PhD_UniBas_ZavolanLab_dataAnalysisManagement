#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-10
#$ -o LOG._Get_region_exp.$TASK_ID

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

date
export PATH=$HOME/bin:$PATH
#-----------------
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
subregion=intron
indir=~bilebi00/_EWSR1/data/clipz_RNAseq_tags/
regionf=~bilebi00/_ARE/Analysis/regions/hg19_GMAP_GENE${subregion}s_refseq.gtf.gz
geneStructuref=~bilebi00/_ARE/Analysis/regions/hg19_GMAP_GENE_refseq.gtf.gz

outdir=Analysis/Expression/Erfc_$subregion
mkdir -p $outdir; pushd $outdir;

ids=(`less $sid_prot_f | grep -w Erfc | cut -f 2 | sort`)

if [ $SGE_TASK_ID -le ${#ids[@]} ]; then

	id=${ids[$((SGE_TASK_ID-1))]}
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	if [ ! -e $indir/$id.bed.gz ]; then
		rsync -t mirz@web08:~/BITER/RNAseq/scratch/bed_files_hg19/$id.bed.gz $indir/.
	else
		echo $indir/$id.bed.gz exists
	fi

	#use reads that full length in the region
	if [ ! -e $id.coverage.gz ]; then
		bedtools intersect -s -a $indir/$id.bed.gz -b $regionf -f 1 | \
			bedtools coverage -a - -b $regionf -counts -s | gzip -c > $id.coverage.gz
		echo coverageBed finished 
	else
		echo $id.coverage.gz exists
	fi

	~bilebi00/bin/R --no-save --args $id.coverage.gz $geneStructuref "'$idn'" count < ~bilebi00/_EWSR1/scripts/plot_intronExpScore_wrtIntron_fromGtfPlus.R > /dev/null

	echo Plot done
	less $id.txt | awk 'NR>1{print "tag\t"$5}' | bedtools groupby -g 1 -c 2,2,2 -o mean,stdev,count  | cut -f 2- | \
		awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn'" type="'Erfc.$subregion'" > freqstat.$SGE_TASK_ID

	touch runfinished.$SGE_TASK_ID 
fi

if [ $SGE_TASK_ID == 1 ]; then
	#wait
	while [ "`ls runfinished* | wc -l`" -lt ${#ids[@]} ]; do
		echo waiting for other runs
		sleep 2m;
	done

	#merge
	statf=all.stats
	echo -e "tag\ttype\tmu\tsd\tN" > $statf
	cat freqstat* >> $statf

	#plot
	~bilebi00/bin/R --no-save --args $statf "Erfc" "" "" < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

	#clean
	rm freqstat* runfinished*
fi

date
echo DONE
