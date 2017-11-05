#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q mpi8_very_long
# $ -q fs_long@@qc_nehalem
# $ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-10
#$ -o LOG._Get_refseqPA_expressionEstimation_4nonFusionGenes.sh.EWSR1.$TASK_ID

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

#TODO get data if not done before
#rsync -t mirz@web08:~/BITER/UniqueGenomicAlignments/$id.bedplus.gz $indir/.
#exit;

export PATH=$HOME/bin:$PATH
project=DIS3L2
project=EWSR1
project=ARE
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
indir=~bilebi00/_CLIP/data/clipz_UniqueGenomicAlignments
RefGene_3UTR_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_3UTR.gtf.gz
RefGene_exon_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_exon.gtf.gz
RefGene_gene_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_gene.gtf.gz
minReadLen=25

ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Expr" { print $2;}' project=$project`)
annots=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 ~ "Expr" { print $10;}' project=$project`)
annot=${annots[0]}

if [ $SGE_TASK_ID -gt ${#ids[@]} ]; then
	exit;
fi

if [ `echo ${annots[@]} | sed 's/ /\n/g' | sort | uniq | wc -l` -gt 1 ]; then
	echo inconsistent annotation categories for project $project ids:${ids[@]};
	exit;
fi

annotStr=" -wP $annot"
if [ $annot == ALL ]; then
	annotStr=" -wvP bacterial|fungus|vector|Markers_Adaptors|viral"
fi

id=${ids[$((SGE_TASK_ID-1))]}

echo Doing $id $indir for $annotStr

outdir=Analysis/Project_$project/RNAseq
mkdir -p $outdir; pushd $outdir

#select representative transcripts according to the project name
if [ $SGE_TASK_ID == 1 ]; then
	#precedence for representative transcript selection is done by 
	# 1. NM_ or NR_ / NM_ if project=ARE
	# 2. RefSeq status(Reviewed or Validated or Inferred or Provisional or Predicted) 
	# 3. length of gene locus (takes longest)
	# 4. length of exonic region (takes longest) / length of 3'UTR region if project=ARE (takes longest)
	# for mRNA and ncRNA
	# omitting fusion genes (i.e. have a dash in the gene name)
	if [ ! -e representative.id.gid ]; then
		if [ $project == ARE ]; then
			zless $RefGene_3UTR_regions
		else 
			zless $RefGene_exon_regions
		fi | \
			awk -F "\t" 'BEGIN{OFS="\t"}{ print $1, $9, $5-$4+1}' | \
			sort -k 1,2 | bedtools groupby -g 1,2 -c 3 -o sum > length		

		~/_PAPD5/scripts/outerJoin.pl $RefGene_gene_regions length 1,9 1,2 '1-9,12' NULL | \
			#corrected on 20121110 to include those non-fusion genes with a dash (-) in gene symbol
#			awk 'BEGIN{OFS="\t";}{ if ($6 < 6 && ($3 == "mRNA" || $3 == "ncRNA") && $10 !~ /-/) { k=1; if ($0 ~ /NR_/) { k=2; } print $10,$12,k,$6,($5-$4+1),$13} } ' | \
		awk 'BEGIN{OFS="\t";}
			{ 
				if (project == "ARE") {
					k=1;
					if ($6 < 6 && $3 == "mRNA" && ($10 ~ /-[1-9]/ || $10 !~ /-/)) {
						print $10,$12,k,$6,($5-$4+1),$13;
					} 
				} else {
					if ($6 < 6 && ($3 == "mRNA" || $3 == "ncRNA") && ($10 ~ /-[1-9]/ || $10 !~ /-/)) {
						k=1; if ($0 ~ /NR_/) { k=2; } print $10,$12,k,$6,($5-$4+1),$13;
					} 
				}
			}' project=$project | \
			sed 's/";*//g' | sort -k 1,1 -k 3,3g -k 4,4g -k 5,5gr -k 6,6gr > refgene.id.gid 
			less refgene.id.gid	| \
			perl -e '$g=""; while(<>){ @t=split; print $t[1],"\t",$t[0],"\n" if ($g ne $t[0]); $g=$t[0];}' \
			> representative.id.gid

		rm length
	else
		echo representative.id.gid exists
	fi

	#get exonic regions
	if [ ! -e representative.bed.gz ]; then
		~/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_exon_regions representative.id.gid 0 0 transcript_id | \
			~/_CLIP/scripts/gtf2bed.pl | \
			sort -k 1,1 -k 2,2g | \
			gzip -c > representative.bed.gz
	else
		echo representative.bed.gz exists
	fi

	touch prep
fi

while [ ! -e prep ]; do
	echo waiting for prep
	sleep 1m
done

#select reads in mRNAseq based on minReadLen and annotStr
if [ ! -e $id.bed.gz ]; then
	zless $indir/$id.bedplus.gz | \
		awk '$NF > minReadLen{ print }' minReadLen=$minReadLen | \
		grep $annotStr | cut -f 1-6 | \
		sort -k 1,1 -k 2,2g | \
		gzip -c > $id.bed.gz
else
	echo $id.bed.gz exists
fi

#map mRNAseq reads onto representative genes exons
if [ ! -e $id.exp ]; then
	zless representative.bed.gz | \
		bedtools map -null 0 -s -o sum -a stdin -b $id.bed.gz | \
		awk 'BEGIN{OFS="\t"}{print $0,$3-$2}' | \
		sort -k 4,4 | bedtools groupby -g 4 -c 8,7 -o sum,sum \
		> $id.exp
else
	echo $id.exp exists
fi

touch run.$SGE_TASK_ID

if [ $SGE_TASK_ID == 1 ]; then
	while [ `ls run* | wc -l` -lt ${#ids[@]} ]; do
		echo waiting for runs;
		sleep 1m
	done
	rm -rf prep run*

	echo running bayseq...
	#TODO
	#merge counts
	#DIS3L2
	if [ $project == DIS3L2 ]; then
		ls *exp | while read i; do 
			id=`basename $i .exp`; less $sid_prot_f | \
		#	awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /oeDIS/ { print $2"\t"$3}' id=$id project=$project; 
			awk -F "\t" '$2==id && $9==project && $3 !~ /HeLa/ && $3 !~ /siDIS/ { print $2"\t"$3}' id=$id project=$project; 
		done
	else
		ls *exp | while read i; do 
			id=`basename $i .exp`; less $sid_prot_f | \
			awk -F "\t" '$2==id && $9==project && $3 ~ /^si/ { print $2"\t"$3}' id=$id project=$project; 
		done
	fi | sort -k 2,2 > selected

	less selected | cut -f 2 | while read id; do
		echo -n " $id.id $id.len $id.exp"
	done | sed 's/ /\t/g' | sed 's/\t//' > input_4_bayseq
	echo >> input_4_bayseq
	files=`less selected | cut -f 1 | sed 's/$/.exp/' `
	paste $files >> input_4_bayseq 

	gidf=representative.id.gid
	nullData=T
	bootStraps=100
	ncl=8
	date
	R --no-save --args input_4_bayseq $gidf $nullData $bootStraps $ncl < ~/_CLIP/scripts/baySeq_dave.R  &> log.bayseq
	date

	R --no-save --args input_4_bayseq_DE.txt .pernucExp < ~/_CLIP/scripts/brandNewPrettyLookingScatter_4baySeq_wColumnSelection.R > /dev/null

fi

echo DONE
