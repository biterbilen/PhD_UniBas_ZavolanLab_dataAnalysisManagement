#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
# $ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-24
#$ -o LOG._Get_xlinkEnrichment_pernuc.sh._tagorder.$TASK_ID

#HOWTO for submission TODO set
#c=""; s=_Get_xlinkEnrichment_pernuc.sh; for i in `seq 1 2`; do sed -e "s/_tagorder/$i/" $s > $s$i.$c; qsub $s$i.$c; echo $i; done

#HOWTO check runs in a batch?
#lines=TODO; grep DONE LOG._Get_xlinkEnrichment_pernuc.sh* | awk  -F ":" '{ print $1; }' | while read i; do wc -l $i; done | grep "^$lines " | awk '{ print $2; }' | while read i; do rm $i; done

#XXX hard coded for 2 parclip as foreground and 4 rnaseq as background
#for i in LOG._Get_xlinkEnrichment.sh.289.pernuc.*; do echo $i `less $i | grep -w waiting -v | wc -l`; done | grep -v -w 7
#-----------------
rm -rf runfinished.$SGE_TASK_ID
export PATH=/import/bc2/home/zavolan/bilebi00/bin:$HOME/bin:$PATH
#-----------------
#XXX set
tagorder=_tagorder; #1..8; 4.EWSR1 8.TIA1
tagorder=2
annot=mRNA #"": ALL
covCut=1
fpath=/import/bc2/home/zavolan/bilebi00/DATA/hg18_ucsc_tracks
tracks=($fpath/representative_transcripts.exons.bed $fpath/representative_transcripts.introns.bed $fpath/mir.gtf.gz $fpath/tRNAs.gtf.gz)
annots=(mRNA intron miRNA tRNA)

indir=/import/bc2/home/zavolan/bilebi00/_xlinkEnrichment/Analysis/UniqueRawData_ALL_woContamination/
sid_prot_f=/import/bc2/home/zavolan/bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((tagorder-1))]}
c=mutatedInCLIP; 

#-----------------
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk '$1==tag { print $2; }' tag=$tag | sort` \
 	`less $sid_prot_f | grep -w mRNAseq | grep -w xlinkBackground | cut -f 2 | sort`)
nucs=(A C G T); 
#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_$c/$tag
mkdir -p $outdir; pushd $outdir;
#-----------------
i=`echo $((SGE_TASK_ID-1)) / ${#nucs[@]} % ${#ids[@]}  | bc`;
id=${ids[$i]};

i=`echo $((SGE_TASK_ID-1)) / 1           % ${#nucs[@]} | bc`;
nuc=${nucs[$i]};

regionfile=${nuc}_${tag}_clusters.gtf.gz
regionfile=${nuc}_${tag}_clusters.gtfplus.gz

#-----------------
#get index of the annot
track=ALL
for i in `seq 0 $((${#annots[@]}-1))`; do
	if [ "$annot" == ${annots[$i]} ]; then
		track=${tracks[$i]}
		break 
	fi 
done

#TODO
#Questions: 
#Do we care about other mutations?
#Is the distribution of mutations different in different regions? what is the best fitted beta dictribution? 
#use wiggle files to reduce model from a zero inflated model to a poisson beta mixture
#search for convergence multiple EM steps 
#how to merge the replicates?

#define regions to make the calculations
if [ $SGE_TASK_ID -le $((2*${#nucs[@]})) ]; then
	r=$id.$nuc
	for s in "-" "+"; do
		echo Doing $r $s strand
		#extract min tag covered regions
		less $indir/DB_all_$id/copies_$s | \
			awk 'BEGIN{OFS="\t"; c=0; }{ if ($4>minCov) { c=c+1; id=c strand; print $1,src,fea,$2+1,$3,$4,strand,".","gene_id \""id"\"; transcript_id \""id"\";"; } }' src=$tag fea=$c strand=$s minCov=$covCut > $r.$s.gtf
		#select mutated CLIP regions from those with original mutation positions
		cat $indir/DB_all_$id/M$nuc?_$s | \
			awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s > $r.$s.bed.tmp
		#annotate with track 
		bedtools intersect -wb -s -a $r.$s.gtf -b $r.$s.bed.tmp | cut -f 1-9 | \
			bedtools intersect -wa -s -a stdin -b $track | sort | uniq > $r.$s.gtf.tmp2 
		mv $r.$s.gtf.tmp2 $r.$s.gtf
		gzip $r.$s.gtf
		echo $r.$s.gtf.gz generated
	done
	touch $r.regionprep
else 
	while [ ! -e ${ids[0]}.T.regionprep ]; do
		echo Waiting for regionprep files 
		sleep 30;
	done
	while	[ "`ls *regionprep | wc -l`" -lt $((2*${#nucs[@]})) ]; do
		echo Waiting for regionprep files 
		sleep 30;
	done	
fi

#-----------------
#calculate the copies_count and mut_count on the regions defined by the CLIP samples separately
for s in "-" "+"; do
	for r in ${ids[0]}.$nuc ${ids[1]}.$nuc; do
		otag=$id.$r
		echo Doing $r $s $id strand
		#project the copies over the region
		cat $indir/DB_all_$id/copies_$s | \
			awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
			bedtools intersect -wao -s -b stdin -a $r.$s.gtf.gz | sort -k1,9 -T ~bilebi00/scratch | \
			awk -F "\t" '$10!="."{ print }' | \
			bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
			awk -F "\t" 'BEGIN { OFS="\t"; } {$9 = $9 " " type " \"" $10 "\";"; print $0; }' type=copies_count | \
			cut -f 1-9 > $otag.copies_count.$s.gtf
		echo Copies done 

		cat $indir/DB_all_$id/M$nuc?_$s | \
			awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
			bedtools intersect -wao -s -b stdin -a $r.$s.gtf.gz | sort -k1,9 -T ~bilebi00/scratch | \
			bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
			#		sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {$9 = $9 " " type " \"" $10 "\";"; print $0; }' type=mut_count | \
			sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {if ($10 >= mx) { $9 = $9 " " type " \"" $10 "\";"; print $0; }}' type=mut_count mx=0 | \
			cut -f 1-9 > $otag.mut_count.$s.gtf
		echo Mutations done

		if [ $nuc == "T" ]; then
#			cat $indir/DB_all_$id/MTC_$s $indir/DB_all_$id/DT_$s | \
			cat $indir/DB_all_$id/MTC_$s | \
				awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
				bedtools intersect -wao -s -b stdin -a $r.$s.gtf.gz | sort -k1,9 -T ~bilebi00/scratch | \
				bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
				#			sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } { $9 = $9 " " type " \"" $10 "\";"; print $0; }' type=mut_count | \
				sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {if ($10 >= mx) { $9 = $9 " " type " \"" $10 "\";"; print $0; }}' type=mut_count mx=0 | \
				cut -f 1-9 > $otag.xlink_count.$s.gtf
			echo xlink done
		else
			echo xlink skipped
		fi
	done
done

#merge copies_count and mut_count
for r in ${ids[0]}.$nuc ${ids[1]}.$nuc; do
	otag=$id.$r
	echo Merging $r $id
	cat $otag.copies_count.[-+].gtf | gzip -c > $otag.copies_count.gtf.gz
	cat $otag.mut_count.[-+].gtf    | gzip -c > $otag.mut_count.gtf.gz
	if [ $nuc == "T" ]; then
		cat $otag.xlink_count.[-+].gtf  | gzip -c > $otag.xlink_count.gtf.gz
	fi
	echo Merge for strands done

	bedtools intersect -wb -f 1 -s -a $otag.copies_count.gtf.gz -b $otag.mut_count.gtf.gz | \
		grep -v "copies_count \"0\";" | cut -f 1-9,18 | sed 's/;\tgene.* mut/; mut/' | gzip -c > \
		$tag.$otag.mut_count.gtf.gz
	echo Copies and mut intersection done
	zless $tag.$otag.mut_count.gtf.gz | \
		perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $tag.`basename $otag .$nuc` $nuc.mut | \
		bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count > $tag.$otag.mut_count.freqstat
	echo Copies and mut statistics done

	if [ $nuc == "T" ]; then
		bedtools intersect -wb -f 1 -s -a $otag.copies_count.gtf.gz -b $otag.xlink_count.gtf.gz | \
			grep -v "copies_count \"0\";" | cut -f 1-9,18 | sed 's/;\tgene.* mut/; mut/' | gzip -c > \
			$tag.$otag.xlink_count.gtf.gz
		echo Copies and xlink intersection done
		zless $tag.$otag.xlink_count.gtf.gz | \
			perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $tag.`basename $otag .$nuc` $nuc.xlink | \
			bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count > $tag.$otag.xlink_count.freqstat
		echo Copies and xlink statistics done
	else
		echo Copies and xlink intersection skipped 
		echo Copies and xlink statistics skipped
	fi
done

touch runfinished.$SGE_TASK_ID

#clean
while [ "`ls -latr runfinished* | wc -l`" -lt 24 ]; do
	echo Waiting for all runs to finish
	sleep 20
done

rm -rf $id*

echo Cleaning done
echo DONE


