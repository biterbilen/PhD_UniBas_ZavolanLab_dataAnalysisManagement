#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
# $ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-3000
# $ -o LOG._Get_xlinkEnrichment_position.DIS3L2
# $ -o LOG._Get_xlinkEnrichment_position.HUR
#$ -o LOG._Get_xlinkEnrichment_position.AGO2

#XXX FIXME (x2) positions with the following comment
#xlink count or xlinks score
#xlinkscore has very low resolution when only xlinked pernucs are used; Analysis/xlinkEnrichment/xlinkPositionEnrichment_tRNA
#xlinkcount is in Analysis/xlinkEnrichment/xlinkPositionEnrichment_tRNA.back_w_xlinkcount

#merge
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=tRNA.pos_enrichment.pdf *pos_enrichment*/*tRNAs*pdf
#-----------------
#XXX TODO set
tag=AGO2
#tag=HuR
#tag=DIS3L2
type=tRNA #there are other parameters
L=72
#-----------------
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#---
export PATH=$HOME/bin:$PATH
#c=clusters
c=Nucleotide
#-----------------
outdir=Analysis/xlinkEnrichment/xlinkPositionEnrichment_$type/$tag
mkdir -p $outdir; pushd $outdir;

pref=xlinked.$L.tRNAs_genename
nuc=T
otag=_${pref}.gtf
d=pernuc_merged_XL.w_type_genename_$type

#-----------------
if [ $SGE_TASK_ID == 1 ]; then
	#get $type annotated pernucs and $type
	#get the substitution value for 0 as the smallest -log10(smallestnonzero_xlinkScore)
	subst=`less ../../annotation/$tag/pernuc_merged_XL.w_type_genename | \
		awk '$7==type{ print; }' type=$type | sort -k 5,5g | \
		perl -e 'while(<>){ @t=split; next if ($t[4] == 0); print -log($t[4])/log(10); exit; }'`;
	less ../../annotation/$tag/pernuc_merged_XL.w_type_genename | \
		awk 'BEGIN{OFS="\t";}{ if ($7==type) { if ($5 == 0) { $5=subst; } else { $5=-log($5)/log(10); } print; }}' type=$type subst=$subs > $d
	less $d | cut -f 8 | sort | uniq | sed 's/";*//g' > xlinked.$type

	#get L=$L clipped tRNAs in $pref.gtf.gz
	zless ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs_genename.gtf.gz | awk '{if ($5-$4+1 == L) { print; } }' L=$L > $L.${type}_genename.gtf
	grep -w -f xlinked.$type $L.${type}_genename.gtf > $pref.gtf
	gzip $pref.gtf

	#get original bins for L=$L tRNAs in bins
	less $d | cut -f 1-6  | \
		bedtools intersect -s -a stdin -b $pref.gtf.gz -wb | sort -k 5,5g | \
		awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*BIN); if (l==L) { if (bin<0) { bin=-bin;} print $1,$2,$3,$4,$5,$6,bin; }}' BIN=$L L=$L | \
		#xlink count or xlinks score
#		cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | \
	sort -k 7,7g | bedtools groupby -i stdin -g 7 -c 5 -o mean | awk 'BEGIN{c=0;}{c=c+$2; print $2,$1,c;}' | \
		tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 1,1g | \
		awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c,L,tag}' tag=original L=$L > bins$L

	#get nuc=$nuc positions along these clipped L=$L tRNAs to be used in shuffle afterwards
	bedtools nuc "-fi" $genomefa -bed $pref.gtf.gz -s -seq | cut -f 1-9,19 | awk 'NR>1{ print; }' | \
		gzip -c > ${pref}.gtf.nucbed.gz
	~bilebi00/_EWSR1/scripts/split_gtf_nucbed.pl ${pref}.gtf.nucbed.gz "$nuc" $otag 0
	gzip -f $nuc$otag;

#XXX
#	R --no-save --args bins < ~bilebi00/_DIS3L2/scripts/bins.R > /dev/null
#	evince bins.txt.pdf
	touch prepfinished
else
	while [ ! -e prepfinished ]; do
		echo Waiting for prep for 1min!
		sleep 1m
	done
fi

#get random bins for L=$L tRNAs -20120411 bug corrected:the shuffling should be done with $nuc positions not all along tRNA
#less $d | bedtools shuffle -i stdin -incl $pref.gtf.gz -g $genomelen > $SGE_TASK_ID.int
less $d | cut -f 1-6 | bedtools shuffle -i stdin -incl $nuc$otag.gz -g $genomelen > $SGE_TASK_ID.int
less $SGE_TASK_ID.int | cut -f 1-6  | \
	bedtools intersect -s -a stdin -b $pref.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*BIN); if (l==L) { if (bin<0) { bin=-bin;} print $1,$2,$3,$4,$5,$6,bin; }}' BIN=$L L=$L | \
#xlink count or xlinks score
#	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | \
sort -k 7,7g | bedtools groupby -i stdin -g 7 -c 5 -o mean | awk 'BEGIN{c=0;}{c=c+$2; print $2,$1,c;}' | \
	tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 1,1g | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c,L,tag}' tag=random.$SGE_TASK_ID L=$L > $SGE_TASK_ID.bins

rm -rf $SGE_TASK_ID.int;

touch $SGE_TASK_ID.runfinished

if [ $SGE_TASK_ID == 1 ]; then
	while [ "`ls *.runfinished | wc -l`" -lt 3000 ]; do
		echo Waiting 1 minute for the other tasks!
		sleep 1m
	done

	rm -rf *runfinished prepfinished

	#~bilebi00/_DIS3L2/scripts/merge_distances.pl "*\.bins" > bins.random
	~bilebi00/_DIS3L2/scripts/merge_distances_confInt.pl "*\.bins" 0.9 > bins$L.random

	rm -rf *\.bins

	~bilebi00/_PAPD5/scripts/outerJoin.pl bins$L bins$L.random 1 1 '1-2,8-10' 0 | sort -k 1,1g > $L.tRNAs.$tag 
	R --no-save --args $L.tRNAs.$tag "Average Crosslink Score (log10)" < ~bilebi00/_DIS3L2/scripts/closest_confInt.R > /dev/null

fi

popd
echo DONE
exit;
