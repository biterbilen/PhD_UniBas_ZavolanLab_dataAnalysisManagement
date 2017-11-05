#!/bin/bash

outdir=Guo
tag=mrna

outdir=Scheiffele
tag=brainStem

pushd $outdir

tmp=$RANDOM
j=0;
for i in final*/*gtf; do 
	echo $i;
#	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/gtf2clsformat.pl $i transcript_id,seqname,strand,start,end,exon_number,gene_id,FPKM,frac,conf_lo,conf_hi > $i.cls_format
#	~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $i.cls_format > $i.cls;
	awk '$6 == "." {print;}' $i.cls_format > $i.trx
	if [ $j == 0 ]; then
		echo if 
		indices=""
		header=""
		cp $i.trx $tmp
	else
		echo else
		~/_PAPD5/scripts/innerJoin.pl $tmp $i.trx 1,6 1,6 '' > $tmp.1 
		mv $tmp.1 $tmp
	fi
	header="$header $i";
	indices="$indices,"$((8+$(($j*11))))
	j=$((j+1))
	echo $header
	echo $indices
done

echo "trxId $header" | sed -e 's/ /\t/g' > final.FPKM
cut -f 1$indices $tmp >> final.FPKM
rm $tmp

mkdir SI; pushd SI;
for of in brainStem cerebellum; do
	mkdir -p $of;
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/extract_exon_inclusion_levels.pl ../finalCufflinks_Sam68_WT_$of/transcripts.gtf.cls > $of/WT_$of
	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/extract_exon_inclusion_levels.pl ../finalCufflinks_Sam68_KO_$of/transcripts.gtf.cls > $of/KO_$of
	echo id gid strand begin end A_SEQ_RAW A_GENE_RAW B_SEQ_RAW B_GENE_RAW | sed 's/ /\t/g' > $of/$of.txt 
	~bilebi00/_PAPD5/scripts/innerJoin.pl $of/KO_$of $of/WT_$of 1 1 '1-7,13-14' >> $of/$of.txt 
	~bilebi00/bin/R --no-save --args KO_vs_WT_$of $of/$of.txt $of EI < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/calculateDifferentialSplicing_SVG.R > /dev/null
#	/usr/bin/convert $of/*sgv $of.
done

of=cerebellum
f1=cerebellum/KO_vs_WT_cerebellum_SI_Values.txt
f2=~bilebi00/_SHORT_READ_ALIGNMENT/data/Splice_arrays/Harald_31_03_2011_mouse_cerebellum/splicing-change-Harald_31_03_2011_mouse_cerebellum-High_Low_with_duplicates.txt.inp
echo -e `less $f1 | head -n 1` `less $f2 | head -n 1` | sed 's/\s/\t/g' | cut -f 1-50 > $of/merged
~bilebi00/_PAPD5/scripts/innerJoin.pl $f1 $f2 1 8 '' | cut -f 1-50 >> $of/merged 
#TODO organize R script
~bilebi00/bin/R --no-save --args $of/merged <  $of/a.R > /dev/null





#-------
#graph equivalence class
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/gtf2clsformat.pl Cuffcompare_$tag/Cuffcompare_$tag.combined.gtf transcript_id,seqname,strand,start,end,exon_number,gene_name,nearest_ref,class_code,gene_id,oId,source,tss_id > combined.$tag.exon
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/get_trx_from_exonClsFormatFile.pl combined.$tag.exon transcript > combined.$tag.trx
i=KO; j=WT;
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/refine_combined.pl Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.refmap Cufflinks_Sam68_${j}_$tag/Cuffcompare_$tag.transcripts.gtf.refmap combined.$tag.trx 
cat combined.$tag.exon combined.$tag.trx.refined > combined.$tag.inp
for i in KO WT; do
	less Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap | awk '$7>0 {print;}' > Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap.refined
#	~bilebi00/_PAPD5/scripts/innerJoin.pl combined.$tag.trx Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap 7,8,9 1,2,3 '1-13,17-26' | awk '{sub(/transcript/,pat); print;}' pat=$i >> combined.$tag.inp
	~bilebi00/_PAPD5/scripts/innerJoin.pl combined.$tag.trx.refined Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap.refined 7,8,9 1,2,3 '1-13,14-26' | awk '{sub(/transcript/,pat); if ($11 == $18) print;}' pat=$i >> combined.$tag.inp
done
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel combined.$tag.inp > combined.$tag.cls

#-----
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/gtf2clsformat.pl Cuffcompare_$tag/Cuffcompare_$tag.combined.gtf transcript_id,seqname,strand,start,end,oId,exon_number,gene_id,nearest_ref,gene_name,class_code,source,tss_id > combined.$tag.exon
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/get_trx_from_exonClsFormatFile.pl combined.$tag.exon transcript > combined.$tag.trx
cat combined.$tag.exon combined.$tag.trx > combined.$tag.inp
for i in KO WT; do
	~bilebi00/_PAPD5/scripts/innerJoin.pl combined.$tag.trx Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap 6 5 '' | awk '{sub(/transcript/,pat); print;}' pat=$i >> combined.$tag.inp
done
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel combined.$tag.inp > combined.$tag.cls
~bilebi00/_SHORT_READ_ALIGNMENT/scripts/parse_cls_for_exon_and_gene_exp.pl combined.$tag.cls > combined.$tag.parsed

#TODO devam et


#r=$RANDOM
for i in KO WT; do
	echo $i
#	less Cufflinks_Sam68_${i}_$tag/Cuffcompare_Sam68.transcripts.gtf.refmap | awk -F "[,\t]" '$5 ~! // {print $1; }' | sort -u > ${i}_$tag.ambiguous.genes
#	~bilebi00/_SHORT_READ_ALIGNMENT/scripts/filter_lines.pl 7 $tag.combined.cls ${i}_$tag.ambiguous.genes > $r.$i 2> /dev/null
#	~bilebi00/_PAPD5/scripts/innerJoin.pl $r.$i Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap 8,7,9 1,2,3 '1-13,17-26' > combined.$tag.$i.tmap

#	~bilebi00/_PAPD5/scripts/innerJoin.pl $tag.combined.cls Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap 8,7,9 1,2,3 '1-13,17-26' > combined.$tag.$i.tmap
	~bilebi00/_PAPD5/scripts/innerJoin.pl $tag.combined.cls Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap 8,11,9 1,5,3 '1-13,17-26' > combined.$tag.$i.tmap
done
#rm $r.*;

~/_SHORT_READ_ALIGNMENT/scripts/gtf2clsformat.pl Cuffmerge_$tag/merged.gtf seqname,strand,start,end,gene_id,transcript_id,nearest_ref,gene_name,class_code,exon_number,oId,source,tss_id > $tag.merged.cls

mergedtotalTCONS=`less $tag.merged.cls | cut -f 6 | sort -u | wc -l`
mergedtotalXLOC=`less $tag.merged.cls | cut -f 5 | sort -u | wc -l`
#TCONS in merged: 61013
#XLOC in merged: 31182

combinedtotalTCONS=`less $tag.combined.cls | cut -f 6 | sort -u | wc -l`
combinedtotalXLOC=`less $tag.combined.cls | cut -f 5 | sort -u | wc -l`
#TCONS in combined: 98392
#XLOC in combined: 57085

guidetotalGenes=`less ~/_SAM68/data/GMAP_EXONS.gtf | awk '{ print $10; }' | sort -u | wc -l`
guidetotalTranscripts=`less ~/_SAM68/data/GMAP_EXONS.gtf | awk '{ print $12; }' | sort -u | wc -l`
#genes in guide: 29214
#transcripts in guide: 122104

echo mergedtotalTCONS: $mergedtotalTCONS
echo mergedtotalXLOC: $mergedtotalXLOC

echo combinedtotalTCONS: $combinedtotalTCONS
echo combinedtotalXLOC: $combinedtotalXLOC

echo guidetotalGenes: $guidetotalGenes
echo guidetotalTranscripts: $guidetotalTranscripts

less $tag.merged.cls  | grep -P "\t1\s\S+$" | cut -f 9 | sort | uniq -c |awk '{ printf("%s\t%.3f\t%d\n", $2,$1/total*100,$1); }' total=$mergedtotalTCONS | sort -k 2,2gr  
less $tag.combined.cls  | grep -P "\t1\s\S+" | cut -f 9 | sort | uniq -c |awk '{ printf("%s\t%.3f\t%d\n", $2,$1/total*100,$1); }' total=$combinedtotalTCONS | sort -k 2,2gr


#get FMI inconsistent transcripts
for i in KO WT; do
	echo $i
	less Cufflinks_Sam68_${i}_$tag/Cuffcompare_$tag.transcripts.gtf.tmap  | cut -f 2,5,6,7,12 | perl -e 'while(<>){@t=split; if ($t[1] eq $t[4] and ($t[2] != 100 && $t[3] > 0.000001) ) { print $_;} }' | sort -u
done

~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $tag.combined.cls 


#class distinguised transcripts
less Cufflinks_Sam68_KO_brainStem/Cuffcompare_brainStem.transcripts.gtf.refmap  | awk 'BEGIN{OFS="\t";} { print $1,$4,$3,$2;}' | sort -k 4,4gr | uniq -f 3 -c | grep "  2 " | awk '{ print $5; }' | while read i; do echo $i; grep $i Cufflinks_Sam68_KO_brainStem/Cuffcompare_brainStem.transcripts.gtf.refmap; echo; done | less


#0 expression level transcripts
for t in KO WT; do echo $t; less Cufflinks_Sam68_${t}_brainStem/*tmap | awk '{ if ($7 == 0 ) { print $2; } }' | while read i; echo $i; do grep -w $i Cufflinks_Sam68_${t}_brainStem/*refmap; echo; done; done | less -S
