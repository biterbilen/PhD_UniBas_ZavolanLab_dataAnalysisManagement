tag=molcel_4184_mmc2_all_GC_skewed_regions
indir=~bilebi00/_EWSR1/data/Ginno2012/
over=~bilebi00/_EWSR1/data/over/hg18ToHg19.over.chain.gz
db=hg19
genomelen=~bilebi00/aux/human.hg19.genome

outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/annotation_xlinkEnrichment_topN10000_cellHEK293/Rloop
mkdir -p $outdir; pushd $outdir;

less $indir/$tag.txt  | grep -v "^#" | awk 'BEGIN{c=0;OFS="\t";}{c=c+1; print $1, $2-1, $3, c"_"$4"skew", 0, "+"; }' > ${tag}_hg18.bed
~bilebi00/www/MARA/bin/liftOver ${tag}_hg18.bed $over $tag.bed $tag.unmapped

#XXX
refseqf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/regions/hg19_GMAP_mRNAs.refseq.gtf.gz
intronf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/regions/hg19_GMAP_GENEintrons_refseq.gtf.gz
clipf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/annotation_xlinkEnrichment_topN10000_cellHEK293/intron/28*.bed.w_common_subregion_genes

rm -rf stats;
l=500
r=0
#get expressed genes extended -500 to the promoter region
#less $refseqf | grep NM_ | \
less $refseqf | \
	bedtools slop -s -l $l -r 0 -i stdin -g $genomelen > ext_genes.gtf

~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl ext_genes.gtf $expf 0 0 > ext_genes.exp
bedtools intersect -u -a molcel_4184_mmc2_all_GC_skewed_regions.bed -b ext_genes.exp > genic_skewed

#get promoter regions of refseq transcripts of expressed genes
#less $refseqf | grep NM_ | \
less $refseqf | \
	awk -F "\t" 'BEGIN{OFS="\t";}{if ($7=="+") {$5=$4;} else {$4=$5;} print $0; }' | \
	bedtools flank -s -l $l -r $r -i stdin -g $genomelen > core_prom.gtf
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl core_prom.gtf $expf 0 0 > core_prom.exp
echo `wc -l core_prom.exp` >> stats
bedtools intersect -u -a molcel_4184_mmc2_all_GC_skewed_regions.bed -b core_prom.exp > core_prom_skewed
echo "core_prom_skewed_over_genic(%)" `less core_prom_skewed | wc -l` `less genic_skewed | wc -l` | awk 'BEGIN{OFS="\t";}{ print $1,100*$2/$3,$3;}' >> stats

#get expressed refseq introns
~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl $intronf $expf 0 1 > introns.exp
all=`zless introns.exp | wc -l`;

#get which Ginno regions are intronic -unstranded
#bedtools intersect -wa -a introns.exp -b $tag.bed | sort | uniq > introns.skewed
#bedtools intersect -wa -a introns.exp -b genic_skewed | sort | uniq > introns.skewed
bedtools intersect -wa -a introns.exp -b core_prom_skewed | sort | uniq > introns.skewed
skewed=`less introns.skewed | wc -l`;

#get clipped introns -stranded
introns_clipped="";
for f in $clipf; do
	of=`basename $f`;
	cat $f | \
#		bedtools intersect -s -wa -a introns.exp -b stdin | sort | uniq > introns.clipped.$of
		bedtools intersect -s -u -a introns.exp -b stdin | sort | uniq > introns.clipped.$of
	introns_clipped="$introns_clipped introns.clipped.$of";
done
cat $introns_clipped | sort | uniq -D | uniq > introns.clipped
rm -rf $introns_clipped
clipped=`less introns.clipped | wc -l`;

#get both skewed and clipped -stranded
bedtools intersect -s -a introns.skewed -b introns.clipped > introns.clipped.skewed

both=`less introns.clipped.skewed | wc -l`;

#echo "testname alternative common skewed clipped expressed p-value oddsRatio" | sed -e 's/ /\t/g' > fisher.stats
echo "testname alternative common core_prom_skewed jointly_clipped expressed p-value oddsRatio" | sed -e 's/ /\t/g' > fisher.stats
R --no-save --args $both $skewed $clipped $all < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> fisher.stats;

#SKEWED TODO
#all=`wc -l core_prom_skewed`
#echo "testname alternative common core_prom_skewed jointly_clipped expressed p-value oddsRatio" | sed -e 's/ /\t/g' > fisher.stats
#R --no-save --args $both $skewed $clipped $all 1 < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> fisher.stats;
#less fisher.stats  | awk 'NR>1{ print "common_over_jointly_clipped_introns(%)", 100*$3/$5,$5; }' >> stats

echo DONE
