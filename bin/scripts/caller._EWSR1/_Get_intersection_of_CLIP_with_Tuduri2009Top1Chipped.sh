tag=shTop1
inf=~bilebi00/_EWSR1/data/Tuduri2009/GSM437526_shTop1.bed.gz
over=~bilebi00/_EWSR1/data/over/hg18ToHg19.over.chain.gz

outdir=Analysis/Top1
mkdir -p $outdir; pushd $outdir;

cp $inf .
gzip -d `basename $inf`
~bilebi00/www/MARA/bin/liftOver `basename $inf .gz` $over $tag.bed $tag.unmapped

#XXX
intronf=../Rloop/introns.exp
clipf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/annotation_xlinkEnrichment_topN10000_exptagCTRL/intronHEK293/28*.bed.w_common_subregion_genes

#get introns of expressed refseq genes restricted to chr6 and chr1
less $intronf | awk '$1=="chr1" || $1 == "chr6" {print}' > introns.exp
all=`zless introns.exp | wc -l`;

#get which regions ie overlapping with Tuduri2009
bedtools intersect -u -a introns.exp -b $tag.bed | sort | uniq > introns.top1
top1=`less introns.top1 | wc -l`;

#get clipped introns
#get clipped introns -stranded
introns_clipped="";
for f in $clipf; do
	of=`basename $f`;
	cat $f | \
		#   bedtools intersect -s -wa -a introns.exp -b stdin | sort | uniq > introns.clipped.$of
	bedtools intersect -s -u -a introns.exp -b stdin | sort | uniq > introns.clipped.$of
	introns_clipped="$introns_clipped introns.clipped.$of";
done
cat $introns_clipped | sort | uniq -D | uniq > introns.clipped
rm -rf $introns_clipped
clipped=`less introns.clipped | wc -l`;

#get both skewed and clipped
bedtools intersect -s -a introns.top1 -b introns.clipped > introns.clipped.top1
both=`less introns.clipped.top1 | wc -l`;

echo "testname alternative common Top1chipped clipped expressedChr1Chr6 p-value oddsRatio" | sed -e 's/ /\t/g' > fisher.stats
R --no-save --args $both $top1 $clipped $all < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> fisher.stats;
echo DONE

