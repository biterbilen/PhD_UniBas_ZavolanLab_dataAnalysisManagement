#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-6
#$ -o LOG._Get_xlinkEnrichment_sub_freqstats.HuR.sh.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_sub_freqstats.AGO2.sh.$TASK_ID

#-----------------
export PATH=/import/bc2/home/zavolan/bilebi00/bin:$HOME/bin:$PATH
#-----------------
#XXX set
tag=HuR
#tag=AGO2
jobcount=6
c=Nucleotide;
#-----------------
#XXX set annot tracks and generalize for hg19 and others
fpath=/import/bc2/home/zavolan/bilebi00/DATA/hg18_ucsc_tracks
tracks="$fpath/representative_transcripts.exons.bed $fpath/representative_transcripts.introns.bed $fpath/mir.gtf.gz $fpath/tRNAs.gtf.gz"
types="mRNA intron miRNA tRNA"
#-----------------

outdir=Analysis/xlinkEnrichment/basedOn$c/$tag
mkdir -p $outdir; pushd $outdir;
#-----------------

gzs=(`ls *xlink_count.gtf.gz | sort`)
gz=${gzs[$((SGE_TASK_ID-1))]}

atag=`basename $gz _count.gtf.gz`
bedtools annotate -s -i $gz -files $tracks -names $types > tmp.$gz
i=10

rm -rf $gz.subregionfreqstats
for type in $types; do
	less tmp.$gz | awk -F "\t" 'NR>1 && $i == 1 { print;}' i=$i | cut -f 1-9 | \
		perl -e '$tag=shift; $type=shift; while(<>){ $_ =~ /copies_count "(\d+).*mut_count "(\d+)";$/; print join("\t", $tag, $type, $2/$1),"\n"; }' $atag $type | \
		bedtools groupby -g 1,2 -c 3,3,3 -o mean,stdev,count >> $gz.subregionfreqstats
	i=$((i+1))
done 

touch $gz.subregionfreqstatsfinished

if [ $SGE_TASK_ID == 1 ]; then
	while [ "`ls *subregionfreqstatsfinished | wc -l`" -lt $jobcount ]; do
		echo waiting subregionfreqstats;
		sleep 2m
	done
	#merge and plot
	echo -e "tag\ttype\tmu\tsd\tN" > subregion.$tag.freqstats
	cat *subregionfreqstats >> subregion.$tag.freqstats
	~bilebi00/bin/R --no-save --args subregion.$tag.freqstats "Mutation Rate" "Comparison of Region Type for Mutation Rate" "" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
	#clean
	rm *subregionfreqstats *subregionfreqstatsfinished
fi

rm tmp.$gz*
echo DONE

