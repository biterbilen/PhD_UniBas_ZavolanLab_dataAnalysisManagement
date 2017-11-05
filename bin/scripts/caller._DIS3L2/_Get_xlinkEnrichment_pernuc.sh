#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-24
#$ -o LOG._Get_xlinkEnrichment_pernuc.sh._tagorder.$TASK_ID

#HOWTO for submission TODO set
#c=""; s=_Get_xlinkEnrichment_pernuc.sh; for i in `seq 1 1`; do sed -e "s/_tagorder/$i/" $s > $s$i.$c; qsub $s$i.$c; echo $i; done

#HOWTO check runs in a batch?
#lines=TODO; grep DONE LOG._Get_xlinkEnrichment_pernuc.sh* | awk  -F ":" '{ print $1; }' | while read i; do wc -l $i; done | grep "^$lines " | awk '{ print $2; }' | while read i; do rm $i; done


#XXX hard coded for 2 parclip as foreground and 4 rnaseq as background
#for i in LOG._Get_xlinkEnrichment.sh.289.pernuc.*; do echo $i `less $i | grep -w waiting -v | wc -l`; done | grep -v -w 7
#-----------------
export PATH=/import/bc2/home/zavolan/bilebi00/bin:$HOME/bin:$PATH
#-----------------
#set
tagorder=_tagorder; #1..8; 4.EWSR1 8.TIA1

indir=/import/bc2/home/zavolan/bilebi00/_DIS3L2/Analysis/UniqueRawData_ALL_woContamination/
sid_prot_f=/import/bc2/home/zavolan/bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((tagorder-1))]}
c=clusters; 

#-----------------
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk '$1==tag { print $2; }' tag=$tag | sort` \
 	`less $sid_prot_f | grep -w mRNAseq | grep -w background | cut -f 2 | sort`)
nucs=(A C G T); 

#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide/$tag
mkdir -p $outdir; pushd $outdir;
#-----------------
i=`echo $((SGE_TASK_ID-1)) / ${#ids[@]} % ${#nucs[@]} | bc`;
nuc=${nucs[$i]};

i=`echo $((SGE_TASK_ID-1)) / 1          % ${#ids[@]}  | bc`;
id=${ids[$i]};
#-----------------
r=$id.$c.$nuc
regionfile=${nuc}_${tag}_$c.gtf.gz
#-----------------
for s in "-" "+"; do
	echo `date` Doing $r $s strand
	zless ../../regions/$regionfile | grep -w "$s" > $r.$s.gtf
	echo $r.gtf generated
	cat $indir/DB_all_$id/copies_$s | \
		awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
		bedtools intersect -wao -s -b stdin -a $r.$s.gtf | sort -k1,9 -T ~bilebi00/scratch | \
		bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
		sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {$9 = $9 " " type " \"" $10 "\";"; print $0; }' type=copies_count | \
		cut -f 1-9 > $r.copies_count.$s.gtf
	echo `date` Copies done 

	cat $indir/DB_all_$id/M$nuc?_$s | \
		awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
		bedtools intersect -wao -s -b stdin -a $r.$s.gtf | sort -k1,9 -T ~bilebi00/scratch | \
		bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
		sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {$9 = $9 " " type " \"" $10 "\";"; print $0; }' type=mut_count | \
		cut -f 1-9 > $r.mut_count.$s.gtf
	echo `date` Mutations done
	if [ $nuc == "T" ]; then
		cat $indir/DB_all_$id/MTC_$s $indir/DB_all_$id/DT_$s | \
			awk 'BEGIN{OFS="\t"; c=0; } { if ($4 > 0) { c=c+1; print $1,$2,$3,c,$4,strand; } }' strand=$s | \
			bedtools intersect -wao -s -b stdin -a $r.$s.gtf | sort -k1,9 -T ~bilebi00/scratch | \
			bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 14 -o sum | \
			sed 's/\t-1$/\t0/' | awk -F "\t" 'BEGIN { OFS="\t"; } {$9 = $9 " " type " \"" $10 "\";"; print $0; }' type=mut_count | \
			cut -f 1-9 >$r.xlink_count.$s.gtf
		echo `date` xlink done
	else
		echo `date` xlink skipped
	fi
done

cat $r.copies_count.[-+].gtf | gzip -c > $r.copies_count.gtf.gz
cat $r.mut_count.[-+].gtf    | gzip -c > $r.mut_count.gtf.gz
if [ $nuc == "T" ]; then
	cat $r.xlink_count.[-+].gtf  | gzip -c > $r.xlink_count.gtf.gz
fi
echo `date` Merge for strands done

bedtools intersect -wb -f 1 -s -a $r.copies_count.gtf.gz -b $r.mut_count.gtf.gz | \
	grep -v "copies_count \"0\";" | cut -f 1-9,18 | sed 's/;\tgene.* mut/; mut/' | gzip -c > \
	$tag.$r.mut_count.gtf.gz
echo `date` Copies and mut intersection done
zless $tag.$r.mut_count.gtf.gz | \
	perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $tag.$r.mut freq | \
	bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count > $tag.$r.mut_count.freqstat
echo `date` Copies and mut statistics done

if [ $nuc == "T" ]; then
	bedtools intersect -wb -f 1 -s -a $r.copies_count.gtf.gz -b $r.xlink_count.gtf.gz | \
		grep -v "copies_count \"0\";" | cut -f 1-9,18 | sed 's/;\tgene.* mut/; mut/' | gzip -c > \
		$tag.$r.xlink_count.gtf.gz
	echo `date` Copies and xlink intersection done
	zless $tag.$r.xlink_count.gtf.gz | \
		perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $tag.$r.xlink freq | \
		bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count > $tag.$r.xlink_count.freqstat
	echo `date` Copies and xlink statistics done
else
	echo `date` Copies and xlink intersection skipped 
	echo `date` Copies and xlink statistics skipped
fi

rm -rf $r*
echo `date` Cleaning done
echo DONE


