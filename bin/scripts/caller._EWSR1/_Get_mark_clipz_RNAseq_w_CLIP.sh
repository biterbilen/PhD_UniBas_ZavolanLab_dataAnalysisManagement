

outdir=Analysis/clipz_RNAseq/FoldChange_woQN
outdir=Analysis/clipz_RNAseq/FoldChange_wSN_woMM_w_clip_w_chip
mkdir -p $outdir; pushd $outdir;

#XXX set
infile=../_RNAseq.raw.qn.w_geneSymbol
infile=../_RNAseq.raw.w_geneSymbol
infile=../_RNAseq_common_scale.raw.w_geneSymbol
file=`basename $infile`
gi=8
i1s=(3 5 7)
i2s=(2 4 6)
clipfiles="../../xlinkEnrichment/annotation/EWSR1/pernuc_merged_XL.gene_intron ../../xlinkEnrichment/annotation/EWSR1/pernuc_merged_XL.gene_mRNA"
chipfile="../../ChIPseq2_severin/outzvals.z3.5.annot"

#pseudocount for log ratio
ps=0.2
r=$RANDOM
#---
#CLIP intersection
cat $clipfiles |  sed 's/^ *//g' | sed 's/";*//g' | sed 's/ /\t/g' | sort -k 2,2 | \
	bedtools groupby -i stdin -g 2 -c 1 -o sum | awk '{print $2"\t"$1;}' > clip
echo `head -n 1 $infile` xlinkCount xlinkStatus | sed 's/ /\t/g' > $file.w_clip 
~bilebi00/_PAPD5/scripts/leftJoin.pl $infile clip $gi 2 '' 0 | \
	awk 'BEGIN{OFS="\t";}{if($(NF-1)>0){$NF="CLIPPED";} else{$NF="NOTCLIPPED";} if(NR>1){print}}' >> $file.w_clip

#---
#ChIP intersection
cat $chipfile | cut -f 15 |  sed 's/^ *//g' | sed 's/";*//g' | sed 's/ /\t/g' | sort -k 2,2 | \
	bedtools groupby -i stdin -g 2 -c 2 -o count | awk '{print $2"\t"$1;}' > chip
echo `head -n 1 $file.w_clip` chipCount chipStatus | sed 's/ /\t/g' > $file.w_clip.w_chip
~bilebi00/_PAPD5/scripts/leftJoin.pl $file.w_clip chip $gi 2 '' 0 | \
	awk 'BEGIN{OFS="\t";}{if($(NF-1)>0){$NF="CHIPPED";} else{$NF="NOTCHIPPED";} if(NR>1){print}}' >> $file.w_clip.w_chip

#calculate log2Ratio
echo -e "tag\ttype\tmu\tsd\tN" > log2Ratio_clip_chip.stats
echo -e "tag\ttype\tmu\tsd\tN" > log2Ratio_clip.stats
echo -e "tag\ttype\tmu\tsd\tN" > log2Ratio_chip.stats
indices='1,2,5,6'
for i in `seq 0 $((${#i1s[@]}-1))`; do
	i1=${i1s[$i]}
	i2=${i2s[$i]}
	indices=$indices","$(($((i+1))*7))

	tag=`less $file.w_clip.w_chip | head -n 1 | cut -f $i1`_over_`less $file.w_clip.w_chip | head -n 1 | cut -f $i2`
	echo $tag
	less $file.w_clip.w_chip | \
		awk 'BEGIN{OFS="\t";}{if(NR==1){print $1,$gi,$i1,$i2,$(NF-2),$NF,t;}else{print $1,$gi,$i1,$i2,$(NF-2),$NF,(log($i1+ps)-log($i2+ps))/log(2);}}' i1=$i1 i2=$i2 gi=$gi ps=$ps t=$tag > \
		$r$i

	less $r$i | awk 'BEGIN{OFS="\t";} { if (NR>1){print;}}' | sort -k 5,5 | \
		bedtools groupby -g 5 -i stdin -c 7,7,7 -o mean,stdev,count | awk 'BEGIN{OFS="\t";}{print $1,type,$2,$3,$4; }' type=$tag >> log2Ratio_clip.stats

	less $r$i | awk 'BEGIN{OFS="\t";} { if (NR>1){print;}}' | sort -k 6,6 | \
		bedtools groupby -g 6 -i stdin -c 7,7,7 -o mean,stdev,count | awk 'BEGIN{OFS="\t";}{print $1,type,$2,$3,$4; }' type=$tag >> log2Ratio_chip.stats

	less $r$i | awk 'BEGIN{OFS="\t";} { if (NR>1){ if ($5=="CLIPPED")  { if ($6=="CHIPPED") { $5="BOTH"; } else { $5=tag1; } } else { if ($6=="CHIPPED") { $5=tag2; } else { $5=tag3;  }}  print; }}' tag1=CLIPPED tag2=CHIPPED tag3=NONE | \
		sort -k 5,5 | \
		bedtools groupby -g 5 -i stdin -c 7,7,7 -o mean,stdev,count | awk 'BEGIN{OFS="\t";}{print $1,type,$2,$3,$4; }' type=$tag >> log2Ratio_clip_chip.stats

done

~bilebi00/bin/R --no-save --args log2Ratio_clip_chip.stats "Expression Fold Change (log2)" "CLIP and CHIP Statistics" "" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
~bilebi00/bin/R --no-save --args log2Ratio_clip.stats "Expression Fold Change (log2)" "CLIP Statistics" "" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
~bilebi00/bin/R --no-save --args log2Ratio_chip.stats "Expression Fold Change (log2)" "ChIP Statistics" "" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

#merge files plot and merge also plots
files=`ls $r* | sort -g`;
echo $files
paste $files | cut -f $indices > $file.w_clip.w_chip.w_log2Ratio
less $file.w_clip.w_chip.w_log2Ratio | cut -f 1,5- > $file.w_log2Ratio
~bilebi00/_DIS3/scripts/matrix2dataframe.pl $file.w_log2Ratio 1-$((${#i1s[@]})) exp > $file.w_log2Ratio.df
R --no-save --args $file.w_log2Ratio $file.w_log2Ratio.df 1 "Expression Fold Change (log2)" < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/splom_density_exp.R > /dev/null

rm $file.w_clip $r* clip chip $file.w_log2Ratio $file.w_log2Ratio.df

#get overlap of top N for fold change
#---
topN=600
infile=$file.w_clip.w_chip.w_log2Ratio
header=(`less $infile | head -n 1`)
rm -rf $infile.common
count=${#i1s[@]}
gi=2
less $infile | awk '$gi!="gid" && $gi!="NA"{print $gi}' gi=$gi > universe.txt
for t in UP DOWN; do
	for i in `seq 0 $((count-1))`; do
		k=$((i+5))
		if [ $t == UP ]; then
			less $infile | sort -k $k,${k}gr > $infile.sorted
		else
			less $infile | sort -k $k,${k}g > $infile.sorted
		fi 
		less $infile.sorted | awk '$gi!="gid" && $gi!="NA"{print $gi}' gi=$gi | head -n $topN > ${header[$((k-1))]}$t;
		rm -rf $infile.sorted
	done
	if [ $count -lt 1 ]; then #deprecated
		R --no-save --args `ls *$t` "" < ~/_DIS3/scripts/venn$count.R > /dev/null;
		mv venn$count.txt $t.tab
		mv venn$count.pdf $infile.$t.pdf
		~bilebi00/_PAPD5/scripts/leftJoin.pl $t.txt $infile 1 1 '' | \
			perl -e '$type=shift; $count=shift; while(<>){@t=split; $sum=0; for (1..$count) { $sum+=$t[$_]; } print $type,"\t",$t[0],"\t",$t[$count+2],"\n" if ($sum == $count); }' $t $count >> $infile.common
	else
		R --no-save --args $t $t universe.txt `ls *$t` < ~bilebi00/_DIS3/scripts/intersectDiagram.R > /dev/null
	fi
	rm *$t
done

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all$file.pdf $file.*.pdf log2*pdf UP.pdf DOWN.pdf 
rm $file.*pdf log2*pdf universe.txt UP.pdf DOWN.pdf

echo DONE
exit

#---
#z-value
#s1=`less $file | awk 'NR>1{ print "aa\t"$0; }' | bedtools groupby -i stdin -g 1 -c $((1+i1)) -o stdev | cut -f 2`
#s2=`less $file | awk 'NR>1{ print "aa\t"$0; }' | bedtools groupby -i stdin -g 1 -c $((1+i2)) -o stdev | cut -f 2`
#n=`wc -l $file | awk '{ print $1-1}'`
#awk 'BEGIN{OFS="\t";}{if (NR==1){print $1,$gi,$i1,$i2,"z-value";}else{print $1,$gi,$i1,$i2,($i1-$i2)/sqrt((s1*s1+s2*s2)/2);}}' i1=$i1 i2=$i2 gi=$gi s1=$s1 s2=$s2 n=$n $file > $file.w_zvalue
#grep -P -w "t-value|DIS3L2|HSPA1B|MAZ|ACTG1|BIRC5|SCD|EEF2|TMBIM6|TMEM107|TUBA1B|ATP2A2|GNB2L1|ATP1A1|CD164|HSPA1A|KDELR1|KDELR1|PMAIP"

