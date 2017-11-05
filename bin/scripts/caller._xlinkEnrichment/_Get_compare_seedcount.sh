#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-6
#$ -o LOG._Get_compare_seedcount.sh.$TASK_ID

#Note:
#reduced xlinkEnrichment to mRNA moltype using representative transcripts exons

export PATH=$HOME/bin:$PATH
#-----------------
#c=clusters
#c=Nucleotide
c=mutatedInCLIP
annot=mRNA
db=hg18
#-----------------
#XXX define mirs
mirs="hsa-let-7a hsa-miR-103 hsa-miR-106a hsa-miR-10a hsa-miR-19a hsa-miR-25 hsa-miR-19a hsa-miR-30a hsa-miR-320a hsa-miR-7"
mirhsa=~bilebi00/CDS_CLIP_ElMMo/data/hsa.fa
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#-----------------
sid_prot_f=~bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
#tag=${tags[$((SGE_TASK_ID-1))]}
tag=AGO2
sids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_$c/$tag/revCompSeedCount
mkdir -p $outdir; pushd $outdir;

methods=( xlinkEnrichment xlinkCount tagCount tagEnrichment mRNAtagCount mRNAtagEnrichment )
method=${methods[$((SGE_TASK_ID-1))]}

jobcount=${#methods[@]}

#set these accordingly
buffer=8000
bk=8
binsize=1000
topNs=(1000 2000 3000 4000 -4000)

#------------------------
#prepare reverse complement of 2-7 seed regions
if [ $SGE_TASK_ID == 1 ]; then
	for i in $mirs; do 
		echo `grep -w $i$ $mirhsa -A 1 | tail -n 1`; 
	done | perl -e 'while(<>){ chomp; $a=substr($_,1,7); $a =~ tr/ACGU/TGCA/ ; $b=reverse $a; print $b,"|";}' > mirseeds_rev_comp
	sed -i 's/|$/\n/g' mirseeds_rev_comp
fi

#------------------------
#prepare input data
if [ $method == mRNAtagCount ]; then
	for N in ${sids[@]}; do
		sn=`less $sid_prot_f | grep -w $N | cut -f 3`;
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			zless ~bilebi00/_xlinkEnrichment/data/MRNAExpression/$sn/windows.foreground.maxf.fa.gz | \
				perl -e 'while(<>){chomp; @t=split/\|/; $_=<>; print join("\t",$t[3],$_); }' | \
				sort -k 1,1gr | head -n $buffer | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == mRNAtagEnrichment ]; then
	for N in ${sids[@]}; do
		sn=`less $sid_prot_f | grep -w $N | cut -f 3`;
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			zless ~bilebi00/_xlinkEnrichment/data/MRNAExpression/$sn/windows.foreground.maxf.fa.gz | \
				perl -e 'while(<>){chomp; @t=split/\|/; $_=<>; print join("\t",$t[3]/$t[4],$_); }' | \
				sort -k 1,1gr | head -n $buffer | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == tagEnrichment ]; then
	for N in ${sids[@]}; do
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			#extend 10 for 20 nucleotide overlap
			less ../copies_count/$tag.$N.$N.T.xlink.bed | head -n $(($buffer*$bk)) | \
				bedtools slop -b 10 -i stdin -g $genomelen > tmp$N.$method.slop 
			#get nonoverlapping windows
			bedtools intersect -s -a tmp$N.$method.slop -b tmp$N.$method.slop -wao | \
				~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl | \
			 	bedtools slop -b 10 -i stdin -g $genomelen > tmp$N.$method.nonover
			#get background copies in these windows
			bfile=../`less ../pbeta.params | grep mu_ref_max | grep ".$N " | awk '{ print $1"."$2;}'`_count.gtf.gz
			bedtools intersect -s -a tmp$N.$method.nonover -b $bfile -wao | \
				perl -e 'while(<>){($a)=$_=~/copies_count "(\d+)"/; print (($a||0),"\t",$_); }' | \
				bedtools groupby -g 2,3,4,5,6,7 -c 1 -o sum -i stdin | \
				awk 'BEGIN{OFS="\t";}{ $5=$7; print; }' | cut -f 1-6 > tmp$N.$method.background
			#get ratio and sequence for the top
			ps=5
			bedtools intersect -wao -f 1 -s -a tmp$N.$method.nonover -b tmp$N.$method.background | \
				awk 'BEGIN{OFS="\t";}{ $5=(($5+ps)/($11+ps)); print; }' ps=$ps | cut -f 1-6 | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > tmp$N.$method.nucbed
			#sort the scores and get buffer size
			awk 'NR>1{print $5"\t"$16;}' tmp$N.$method.nucbed | sort -k 1,1gr | head -n $buffer | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == tagCount ]; then
	for N in ${sids[@]}; do
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			less ../copies_count/AGO2.$N.$N.T.xlink.bed | head -n $(($buffer*$bk)) | \
				bedtools slop -b 10 -i stdin -g $genomelen > tmp$N.$method.slop 
			bedtools intersect -s -a tmp$N.$method.slop -b tmp$N.$method.slop -wao | \
				~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl | head -n $buffer | \
				bedtools slop -b 10 -i stdin -g $genomelen | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > tmp$N.$method.nucbed
			awk 'NR>1{print $5"\t"$16;}' tmp$N.$method.nucbed | sort -k 1,1gr | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == xlinkCount ]; then
	for N in ${sids[@]}; do
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			less ../mut_count/AGO2.$N.$N.T.xlink.bed | head -n $(($buffer*$bk)) | \
				bedtools slop -b 10 -i stdin -g $genomelen > tmp$N.$method.slop 
			bedtools intersect -s -a tmp$N.$method.slop -b tmp$N.$method.slop -wao | \
				~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl | head -n $buffer | \
				bedtools slop -b 10 -i stdin -g $genomelen | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > tmp$N.$method.nucbed
			awk 'NR>1{print $5"\t"$16;}' tmp$N.$method.nucbed | sort -k 1,1gr | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == xlinkEnrichment ]; then
	for N in ${sids[@]}; do
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			less ../xlinkEnrichment/AGO2.$N.$N.T.xlink_mutratemax.bed | head -n $(($buffer*$bk)) | \
				bedtools slop -b 10 -i stdin -g $genomelen > tmp$N.$method.slop 
			bedtools intersect -s -a tmp$N.$method.slop -b tmp$N.$method.slop -wao | \
				~bilebi00/_xlinkEnrichment/scripts/get_nonoverlapping_from_top.pl | head -n $buffer | \
				bedtools slop -b 10 -i stdin -g $genomelen | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > tmp$N.$method.nucbed
			awk 'NR>1{print $5"\t"$16;}' tmp$N.$method.nucbed | sort -k 1,1g | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
else
	echo Unknown method; 
	exit;
fi

#calculate seed match count in each bin
echo -e "tag\tmu\tsd\tN\ttype" > $method.stats
for topN in ${topNs[@]}; do
	for N in ${sids[@]}; do
		if [ $topN -lt 0 ]; then 
			zless tmp$N.$method.gz | tac | head -n $topN | gzip -c > tmp$N.$method.slop.nucbed.gz
		else 
			zless tmp$N.$method.gz | head -n $topN | tail -n $binsize | gzip -c > tmp$N.$method.slop.nucbed.gz
		fi
		~bilebi00/_xlinkEnrichment/scripts/calc_multiplePatternCount_from_nucbed.pl tmp$N.$method.slop.nucbed.gz "`cat mirseeds_rev_comp`" 1 | \
			awk '{print tag"\t"$3;}' tag=$topN | \
			bedtools groupby -g 1 -c 2,2,2 -o mean,stdev,count -i stdin | \
			awk '{ print $0"\t"type;}' type=$method-$N >> $method.stats
	done
done

~bilebi00/bin/R --no-save --args $method.stats "Seed Reverse Complement Count" "" "Top N" < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

touch $SGE_TASK_ID.runfinished

if [ $SGE_TASK_ID == 1 ]; then
	if [ "`ls *runfinished | wc -l`" -lt $jobcount ]; then
		sleep 1m
	fi

	cat *stats | sort -r | uniq > $tag.seedCount_comparison
	~bilebi00/bin/R --no-save --args $tag.seedCount_comparison "Seed Reverse Complement Count" "Seed Reverse Complement Count Comparison" "Top N" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

	rm -rf *runfinished
fi

rm -rf tmp*\.$method\.*
echo DONE

exit;

