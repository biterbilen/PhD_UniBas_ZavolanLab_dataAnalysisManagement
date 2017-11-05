#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 3-5
#$ -o LOG._Get_xlinkEnrichment_annotation_4top.$TASK_ID

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
#XXX set
exptag=Untreated #HEK293 #CTRL
cell=HEK293 #""
topN=10000
c=xlinkEnrichment
regionsf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/regions
clipz_RNAseq_f=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_$exptag/RNAseq.w_geneSymbol
updownext=2000
db=hg19
genomelen=~bilebi00/aux/human.$db.genome
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") {print $2;}' tag=$tag cell="$cell" | sort`)

#-----------------
#link files
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_${c}_topN${topN}_cell${cell}_w_downGSR
mkdir -p $outdir; pushd $outdir;
ln -sf $regionsf .
ln -sf $clipz_RNAseq_f clipz_RNAseq_f

#-----------------
#extract introns exons and upstream and downstream of expressed genes for doing the annotations
less clipz_RNAseq_f  | awk '$1 ~ /NM_/{ print $NF; }' | sort | uniq > exp_genes
annotf=exp_gene_subregions.gtf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENE_refseq.gtf.gz exp_genes | \
	bedtools flank -s -l $updownext -r 0 -i stdin -g $genomelen | awk -F "\t" 'BEGIN{OFS="\t";}{$3=c; print;}' c=upGSR \
		> $annotf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENE_refseq.gtf.gz exp_genes | \
	bedtools flank -s -l 0 -r $updownext -i stdin -g $genomelen | awk -F "\t" 'BEGIN{OFS="\t";}{$3=c; print;}' c=downGER \
		>> $annotf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENE_refseq.gtf.gz exp_genes | \
	bedtools flank -s -l 0 -r $updownext -i stdin -g $genomelen | awk -F "\t" 'BEGIN{OFS="\t";}{$3=c; print;}' c=downGSR \
		>> $annotf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENE_refseq.gtf.gz exp_genes | \
	bedtools flank -s -l $updownext -r 0 -i stdin -g $genomelen | awk -F "\t" 'BEGIN{OFS="\t";}{$3=c; print;}' c=upGER \
		>> $annotf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENEexons_refseq.gtf.gz exp_genes >> $annotf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl regions/hg19_GMAP_GENEintrons_refseq.gtf.gz exp_genes >> $annotf


gzip $annotf
annotf=$annotf.gz

#-----------------
#get how many transcripts from each gene type exists
gene_counts=gene.counts
echo -e "type\tgeneCount" > $gene_counts
zless $annotf | cut -f 3,9 | sort | uniq | bedtools groupby -i - -g 1 -c 2 -o count >> $gene_counts
types=(`less $gene_counts | grep -v "^type" | cut -f 1`)


#-----------------
#gets topN of mu_refmax
statf=position.stats
echo -e "lib\ttype\tvalue\tN" > $statf
substatf=position_4subregion.stats
echo -e "lib\ttype\tvalue\tN" > $substatf
ftags=()
idns=()
for id in ${idsf[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	idns=("${idns[@]}" "'$idn'")
	#input file
	file=../$c/$tag.$id.$id.T.xlink_mu_refmax.bed.gz
	ln -sf $file .

	#output file for topN
	bedgz=$id.bed.w_type_genename.gz

	#attach types and names of intersecting transcripts to top $topN regions;
	zless $file | head -n $topN | \
		bedtools intersect -wao -s -a stdin -b $annotf | \
		awk 'BEGIN{OFS="\t";} { if($7 !~ /chr/) { $9="other"; } print $1,$2,$3,$4,$5,$6,$9,$16; }' | \
		sort -T ~bilebi00/scratch/ | uniq | sed 's/";*//g' | \
		sort -k 5,5g -T ~bilebi00/scratch/ | gzip -c > $bedgz

	#separate transcript names with type
	files4int=()
	for type in ${types[@]}; do
		zless $bedgz | awk '{ if ($7 == type) { print; } }' type=$type | cut -f 8 | sort -T ~bilebi00/scratch/ | \
			bedtools groupby -i - -g 1 -c 1 -o count | sort -k 2,2gr | \
			sort -T ~bilebi00/scratch/ -k 2,2gr | sed 's/";*//g' > $id.gene_$type
		files4int=(${files4int[@]} $id.gene_$type) 
	done
	R --no-save --args ${#files4int[@]} ${files4int[@]} $id.gene_regions_intersectDiagram $id.gene_ < ~bilebi00/_EWSR1/scripts/intersectDiagram_woUniverse.R > /dev/null
	echo intersectDiagram_woUniverse done

	#get distinct xlinked gene counts
	echo -e "type\txlinkedGeneCount" > $id.xlinkedGene.counts
	zless $bedgz | cut -f 7-8 | sort -T ~bilebi00/scratch/ | uniq | \
		bedtools groupby -g 1 -c 2 -o count | sort -T ~bilebi00/scratch/ -k 2,2gr >> $id.xlinkedGene.counts

	#get distinct xlink counts
	total=`zless $bedgz | cut -f 1-6 | sort | uniq | wc -l`;
	echo -e "type\txlinkCount\ttotalXlinkCount" > $id.xlink.counts
	zless $bedgz | sort -T ~bilebi00/scratch/ -k 7,7 | \
		bedtools groupby -g 7 -c 7 -o count | sort -T ~bilebi00/scratch/ -k 2,2gr | awk '{ print $1"\t"$2"\t"total;}' total=$total >> $id.xlink.counts

	R --no-save --args $id.xlink.counts $id.xlinkedGene.counts $gene_counts "$idn" $id. < ~bilebi00/_DIS3L2/scripts/annot2D.R > /dev/null

	#write in stat file 
	Nt=`zless $file | wc -l`;
	lib="$idn"
	N=`zless $bedgz | cut -f 1-6 | sort | uniq | wc -l`
	echo | awk 'BEGIN{OFS="\t";}{ print lib,type,n/N*100,N;}' lib="$lib" type="All Positions" n=$N N=$Nt >> $statf;
	N=`zless $bedgz | awk '$7!="other"{print;}' | cut -f 1-6 | sort | uniq | wc -l`
	echo | awk 'BEGIN{OFS="\t";}{ print lib,type,n/N*100,N;}' lib="$lib" type="mRNA Positions" n=$N N=$Nt >> $statf;
	zless $bedgz | awk '$7!="other" { print $7; }' | sort | bedtools groupby -i - -g 1 -c 1 -o count | \
		awk 'BEGIN{OFS="\t";}{ print lib,$1,$2/N*100,N;}' lib="$lib" N=$N >> $substatf 
	echo statf done
done

echo $statf and $substatf generated

#-----------------
topN=100
stack=0
R --no-save --args $statf $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
R --no-save --args $substatf $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
if [ $tag == EWSR1 ] && [ "`less $statf | grep 90kDa`" != "" ]; then
	less $statf | grep -v "90kDa" > $statf.288.289
	R --no-save --args $statf.288.289 $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
	less $substatf | grep -v "90kDa" > $substatf.288.289
	R --no-save --args $substatf.288.289 $topN $stack < ~bilebi00/_EWSR1/scripts/barchart_summary.R > /dev/null
fi

echo $statf and $substatf barcharts generated

#-----------------
for type in ${types[@]}; do
	ifiles=""
	for id in ${idsf[@]}; do
		ifiles="$ifiles,$id.gene_$type"
	done
	ifiles=`echo $ifiles | cut --bytes 2-`;
	R --no-save --args "$ifiles" exp_genes venn.$type "${idns[@]}" < ~bilebi00/_DIS3/scripts/venn.R > /dev/null
done

echo $statf and $substatf venns generated

echo DONE

