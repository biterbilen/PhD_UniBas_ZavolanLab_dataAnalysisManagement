
outdir=~bilebi00/_EWSR1/Analysis/GO_4clipzRNAseqFCup_down
pushd $outdir

#XXX apply shiv.m to the RNAdata and produce 
#in /Volumes/Data/IBM_import/_Zavolan_Group_Meetings/EWSR1/RNAseq_paper_material/HW3
#RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock.w_Gaussian
file=RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock.w_Gaussian
ifile=$file.w_outlier.w_Kishore
~bilebi00/_EWSR1/scripts/label_outlier_regulation.pl $file 1.61 > $ifile
R --no-save --args $ifile < ~bilebi00/_EWSR1/scripts/gaussian.R > /dev/null
less $ifile | awk 'NR>1{print}'| sort -k 8,8 -k 6,6 -k 9,9 | bedtools groupby -g 8,6,9 -c 5,5 -o mean,stdev -i stdin | less | sort -k 4,4g > $ifile.stats

echo "globalsetname set1name set2name testname alternative #ofcommonGenes #ofgenesInKishoreCombinedSet #ofgenesInHancock #ofgenesinglobal p-value oddsRatio" | \
	sed 's/ /\t/g' > $ifile.Fisher.stats
for ol in Outlier notOutlier; do
	all=`less $ifile | grep -w $ol | wc -l`; 
	for kr in kishoreUp kishoreDown; do
		s1=`less $ifile | grep -w $ol | grep -w $kr | wc -l`; 
		for hc in up down; do
			s2=`less $ifile | grep -w $ol | grep -w $hc | wc -l`; 
			b=`less $ifile | grep -w $ol | grep -w $kr | grep -w $hc | wc -l`; 
			echo $ol $kr $hc `R --no-save --args $b $s1 $s2 $all < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"`
		done
	done
done | sort -k 10,10g | sed -e 's/ /\t/g' >> $ifile.Fisher.stats
#------------

less RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock | awk 'NR>1&&$6!=exc{print;}' exc=None | sort -k 5,5g | \
 	awk 'BEGIN{c=0;k=0;OFS="\t";}{c=c+1;if($6==type) { k=k+1;} print $0, k, c, k/c; }' type=down | head -n 120 | grep -w down > mz_down
less RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock | awk 'NR>1&&$6!=exc{print;}' exc=None | sort -k 5,5gr | \
 	awk 'BEGIN{c=0;k=0;OFS="\t";}{c=c+1;if($6==type) { k=k+1;} print $0, k, c, k/c; }' type=up | head -n 15 | grep -w up > mz_up
#------------

exit;
echo "sample1 sample2 testName testAlternative #ofGenesInBoth #ofGenesInSample1 #ofGenesInSample2 #ofGenesExpressed_in_HeLa p-value oddsRatio" | sed -e 's/ /\t/g' > Hancock2008.top.stats

for i in batch*/studyset; do 
	ii=`perl -e '$_=shift; $_ =~ /([^\/]*)/ ; print $1' $i`;
	is=`less $i |wc -l`; 
	for jj in ~/_EWSR1/data/HancockCC7/*; do 
		j=`basename $jj`;
		#innerJoin.pl $jj population 1 1 1 | head -n $is > $j
		innerJoin.pl $jj population 1 1 1 > $j
		js=`less $j | wc -l`; 
		isjs="$ii-$j"
		innerJoin.pl $i $j 1 1 1 > $isjs;
		b=`less $isjs | wc -l`;
		echo -en "$i\t$j\t"
		R --no-save --args $b $is $js `less population | wc -l` < ~/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"; 
	done; 
done | sort -k 9,9g | sed 's/\/studyset//g' >> Hancock2008.top.stats

less Hancock2008.top.stats | awk 'NR>1{print $1"-"$2}' | while read i; do  echo $i; cat $i; echo; rm $i; done > Hancock2008.top.stats.common.genes

less ST7.EWS-FLI_Upregulated_signature | awk '{print $1"\tup"}' > ST7-8
less ST8.EWS-FLI_Downregulated_signature | awk '{print $1"\tdown"}' >> ST7-8
~/_PAPD5/scripts/outerJoin.pl RNAseq.raw.qn.w_geneSymbol.FC ST7-8 2 1 '' None | \
	awk 'BEGIN{OFS="\t";}{ if (NR>1) {$5=log(2^$3+2^$4)/log(2)-1; } else { $5="averageFC";} print; }'  > RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock
sed -i 's/averageFC\tNone/averageLog2FC\tHancock2008Regulation/g' RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock;
rm ST7-8

echo -e "type\ttag\tmu\tsd\tN" > Hancock2008.FC.stats
awk 'NR>1{ print }' RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock | sort -k 6,6 | bedtools groupby -g 6 -c 3,3,3 -o mean,stdev,count -i stdin | awk '{print"batch1\t"$0}' >> Hancock2008.FC.stats
awk 'NR>1{ print }' RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock | sort -k 6,6 | bedtools groupby -g 6 -c 4,4,4 -o mean,stdev,count -i stdin |  awk '{print"batch2\t"$0}' >> Hancock2008.FC.stats
awk 'NR>1{ print }' RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock | sort -k 6,6 | bedtools groupby -g 6 -c 5,5,5 -o mean,stdev,count -i stdin |  awk '{print"batch12\t"$0}' >> Hancock2008.FC.stats

R --no-save --args Hancock2008.FC.stats "Fold Change (log2)" "" < ~/_EWSR1/scripts/CI.R > /dev/null
