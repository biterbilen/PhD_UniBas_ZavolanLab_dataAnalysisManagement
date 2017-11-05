#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-7
#$ -o LOG._Get_compare_affinity.Nucleotide.sh.$TASK_ID
# $ -o LOG._Get_compare_affinity.clusters.sh.$TASK_ID

export PATH=$HOME/bin:$PATH
#SGE_TASK_ID=2
#-----------------
#c=clusters
#c=Nucleotide
c=mutatedInCLIP
annot=mRNA
db=hg18
#-----------------
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#-----------------
sid_prot_f=~bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
#tag=${tags[$((SGE_TASK_ID-1))]}
tag=HuR
sids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
snas=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 3 | sort`)
#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_$c/$tag/affinity
mkdir -p $outdir; pushd $outdir;

methods=( xlinkEnrichment xlinkCount tagEnrichment tagCount mRNAtagCount mRNAxlinkCount mRNAtagEnrichment )
method=${methods[$((SGE_TASK_ID-1))]}
jobcount=${#methods[@]}

buffer=20000
bk=13
binsize=2500
topNs=(2500 5000 7500 -7500)

if [ $method == mRNAtagCount ]; then
	for N in ${sids[@]}; do
		sn=`less $sid_prot_f | grep -w $N | cut -f 3`;
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			zless ~bilebi00/_xlinkEnrichment/data/MRNAExpression/$sn/windows.foreground.maxf.tab.gz | \
			awk 'NR>1{print $4"\t"$10;}' | \
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
			zless ~bilebi00/_xlinkEnrichment/data/MRNAExpression/$sn/windows.foreground.maxf.tab.gz | \
				awk 'BEGIN{OFS="\t";}{if(NR>1){$4=$4/$5;print $4"\t"$10;}}' | \
				sort -k 1,1gr | head -n $buffer | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
elif [ $method == mRNAxlinkCount ]; then
	for N in ${sids[@]}; do
		sn=`less $sid_prot_f | grep -w $N | cut -f 3`;
		if [ ! -e $N.$method.gz ]; then
			echo Generating $N.$method.gz
			zless ~bilebi00/_xlinkEnrichment/data/MRNAMutationsInit/$sn/mRNA_sites.top20000.tab.gz | \
				awk 'NR>1{print $5"\t"$6;}' | \
				sort -k 1,1gr | head -n $buffer | gzip -c > $N.$method.gz
		else
			echo $N.$method.gz exists
		fi
		ln -sf $N.$method.gz tmp$N.$method.gz
	done
#	for N in 1 2; do
#	zless ~bilebi00/_xlinkEnrichment/data/MRNAMutationsInit/${sids[$((N-1))]}/mRNA_sites.top20000.tab.gz | \
#		awk 'NR>1{print $5"\t"$6;}' | sort -k 1,1gr | gzip -c > tmp$N.$method.gz
#	done
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
			less ../copies_count/HuR.$N.$N.T.xlink.bed | head -n $(($buffer*$bk)) | \
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
			less ../mut_count/HuR.$N.$N.T.xlink.bed | head -n $(($buffer*$bk)) | \
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
			less ../xlinkEnrichment/HuR.$N.$N.T.xlink_mutratemax.bed | head -n $(($buffer*$bk)) | \
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

#calculate mean and stdev of affinities in each bin
echo -e "tag\tmu\tsd\tN\ttype" > $method.stats
for topN in ${topNs[@]}; do
	for N in 003 004; do
		if [ $topN -lt 0 ]; then 
			zless tmp$N.$method.gz | tac | head -n $topN | gzip -c > tmp$N.$method.slop.nucbed.gz
		else 
			zless tmp$N.$method.gz | head -n $topN | tail -n $binsize | gzip -c > tmp$N.$method.slop.nucbed.gz
		fi
		~bilebi00/_xlinkEnrichment/scripts/calc_affinity_from_nucbed.pl ~bilebi00/_xlinkEnrichment/data/$tag.affinity tmp$N.$method.slop.nucbed.gz 1 | \
			awk '{print tag"\t"$3;}' tag=$topN | \
			bedtools groupby -g 1 -c 2,2,2 -o mean,stdev,count -i stdin | \
			awk '{ print $0"\t"type;}' type=$method-$N >> $method.stats
	done
done

~bilebi00/bin/R --no-save --args $method.stats "Affinity" "" "Top N" < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

touch $SGE_TASK_ID.runfinished

if [ $SGE_TASK_ID == 1 ]; then
	if [ "`ls *runfinished | wc -l`" -lt $jobcount ]; then
		sleep 1m
	fi

	cat *stats | sort -r | uniq > $tag.affinity_comparison
	~bilebi00/bin/R --no-save --args $tag.affinity_comparison "Affinity" "Affinity Comparison -ALL" "Top N" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null

	rm -rf *runfinished
fi


rm -rf tmp*\.$method\.*
echo DONE

exit;

