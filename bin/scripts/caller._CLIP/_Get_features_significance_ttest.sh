#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-200
#$ -o LOG._Get_features_significance.$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1";
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV
export PATH=$HOME/bin:$PATH
#-----------------
verbose=1
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
over=~bilebi00/_EWSR1/data/over/hg18ToHg19.over.chain.gz
subregion=intron
tag=EWSR1
cell=HEK293

date
#TODO set which type of feature
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_xlinkEnrichment_topN10000_cell$cell/features_ttest_$subregion
#outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/motif_topN_non_overlapping/features_significance
mkdir -p $outdir; pushd $outdir;

idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort` joint)

statf=feature.stats.$SGE_TASK_ID
features=(SNP ChIPseq DNaseClustered HeLatopxlink FETtopxlink HeLa FET RegionExpression Histone TfbsClustered RHM GC_pct G_pct Conservation Rloop Top1)
featureIs=(11 11 11 11 11 11 11 11 11 11 20 11 14 11 11 11)
featureJobCounts=(1 5 1 1 10 1 10 10 30 200 10 1 1 1 1 1)

#feature specific files
#----------------------
SNPf=~bilebi00/DATA/Ensembl/Annotation/snp135Common.gtf.gz
ChIPseqf=~bilebi00/_EWSR1/Analysis/ChIPseq*_severin*/outzvals.z3
HeLatopxlinkf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/motif_topN_non_overlapping/foreground984.nucbed.gz
FETtopxlinkf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1_FlpIn/motif_topN_non_overlapping/
HeLaf=~bilebi00/_CLIP/Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/EWSR1/motif_topN_non_overlapping/foreground984.nucbed.gz
FETf=~bilebi00/DATA/Hoell_et_al_2011_SRA025082/nsmb.2163-S2.txt
TfbsClusteredf=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClustered.bed.gz
DNaseClusteredf=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/wgEncodeRegDnaseClustered/wgEncodeRegDnaseClustered.bed.gz
RegionExpressionf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/RegionExpression_intron
Histonef=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/Histone_tracks.4download
Rloopf=~bilebi00/_EWSR1/data/Ginno2012//molcel_4184_mmc2_all_GC_skewed_regions.txt
Top1f=~bilebi00/_EWSR1/data/Tuduri2009/GSM437526_shTop1.bed.gz
Conservationf=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/regions/hg19_GMAP_GENE${subregion}s_refseq_wPhastCons_opmean_ws1.gtf.gz

#get positive and negative sets
if [ $SGE_TASK_ID == 1 ]; then
	for id in ${idsf[@]}; do
		if [ ! -e $id.CLIPPED ] || [ ! -e $id.NOTCLIPPED ]; then
			zless ../$subregion/$id.CLIPPED | \
				awk -F "\t" 'BEGIN{OFS="\t"}{ $2=src; print $0; } ' src=$id.CLIPPED
			zless ../$subregion/$id.NOTCLIPPED | \
				awk -F "\t" 'BEGIN{OFS="\t"}{ $2=src; print $0; } ' src=$id.NOTCLIPPED
		fi
	done > input.gtf
	touch dataprep
	echo dataprep done 
fi

while [ ! -e dataprep ]; do
	echo waiting for dataprep;
	sleep 1m;
done

#-----------------------------------
#subfeature arrays
RHM=(CC.CC.T..CC.C CCTCCCT CCCCACCCC G.GG..A.GG.GG AGGGAGG GGGGTGGGG .GGG.GG TGGGTGG GGGGAGGGG)
G_pct=(G_pct)
GC_pct=(GC_pct)
Rloop=($Rloopf)
Top1=($Top1f)
DNaseClustered=($DNaseClusteredf)
HeLatopxlink=($HeLatopxlinkf)
HeLa=($HeLaf)
Conservation=($Conservationf)
FETtopxlink=(`ls $FETtopxlinkf/foreground* `)
RegionExpression=(`ls $RegionExpressionf/*txt`)
ChIPseq=(`ls $ChIPseqf`)
SNP=($SNPf)

dir=`dirname $Histonef`;
Histone=(`less $Histonef | grep -P "H3K4me1|H3K27ac|H3K9me3|H3K27me|H3K4me3|H3K36me3" | awk '{print dir"/"$1"/"$2;}' dir=$dir`)

#takes time so only evaluated when set as a feature 
if [ "` echo ${features[@]} | grep -w FET`" != "" ]; then
	FET=(`zless $FETf | awk '$1~/chr/{ print $4}' | sort | uniq`)
fi

if [ "` echo ${features[@]} | grep -w TfbsClustered`" != "" ]; then
#	TfbsClustered=(`zless $TfbsClusteredf | head -n 1 | cut -f 4 | sort | uniq`)
	TfbsClustered=(`zless $TfbsClusteredf | cut -f 4 | sort | uniq`)
fi

#-------------------------------------
if [ ! -e $statf ]; then
#	echo -e "lib\ttag\ttype\tvalue" > $statf
	touch $statf
fi

#Intersection with the features
for i in `seq 1 ${#features[@]}`; do
	f=${features[$((i-1))]};
	fii=${featureIs[$((i-1))]};
	#"It's worth to note, that even an index will be substituted at time the variable is evaluated:"
	#Bash indirect array addressing
	subfeatureRef=$f[$((SGE_TASK_ID-1))]
	subfeature=${!subfeatureRef}
	#	if [ $SGE_TASK_ID -gt ${featureJobCounts[$((i-1))]} ] || [ "$subfeature" == "" ]; then
	if [ "$subfeature" == "" ]; then
		continue;
	fi
	xx=`basename $subfeature | awk -F "." '{print $1}'`
	normFeature=0
	stranded="-s"
	featureScoreSign=1
	featureScore=-1
	featureName=0
	featureStartBase=0
	addFeatureStrand=0
	featureLiftOver=0
	featureConvertFromGtf=0
	featureNuc=0
	if [ $f == ChIPseq ]; then
		stranded=0
		xx="`echo $subfeature | perl -e '$_=<>; $_=~/(ChIPseq[^_]*)/; print $1; '`"
	elif [ $f == RegionExpression ]; then
		featureStartBase=1
		id=`basename $subfeature .txt`
		xx="`less $sid_prot_f | awk -F "\t" '$2==id{print $6}' id=$id`";
	elif [ $f == Histone ]; then
		stranded=0
		xx="`perl -e '$a=shift; $a=~/Histone([\w\d]*)\./; print $1,"\n"' $subfeature`";
	elif [ $f == FETtopxlink ] || [ $f == HeLatopxlink ]; then
		featureScoreSign=-1
		xx="`perl -e '$a=shift; $a=~/foreground([\d]*)\./; print $1,"\n"' $subfeature | while read id; do less $sid_prot_f | awk -F "\t" '$2==id{print $6;}' id=$id; done`";
	elif [ $f == HeLa ]; then
		featureScore=1000
		xx="`perl -e '$a=shift; $a=~/foreground([\d]*)\./; print $1,"\n"' $subfeature | while read id; do less $sid_prot_f | awk -F "\t" '$2==id{print $6;}' id=$id; done`";
	elif [ $f == FET ]; then
		featureName=$subfeature
		fileref=${f}f
		subfeature=${!fileref}
	elif [ $f == DNaseClustered ]; then
		stranded=0
		addFeatureStrand=1
	elif [ $f == TfbsClustered ]; then
		stranded=0
		featureName=$subfeature
		fileref=${f}f
		subfeature=${!fileref}
	elif [ $f == Rloop ]; then
		featureScore=1000
		addFeatureStrand=1
		stranded=0
		featureLiftOver=$over
	elif [ $f == Top1 ]; then #XXX introns are reduced to chr1 and chr6 
		stranded=0
		featureScore=1000
		featureLiftOver=$over
	elif [ $f == SNP ]; then
		featureConvertFromGtf=6
		featureScore=1000
		featureStartBase=1
	elif [ $f == Conservation ]; then
		featureConvertFromGtf=10
		featureStartBase=1
	elif [ $f == RHM ]; then
		xx="$subfeature"
		featureName=$subfeature
		subfeature=tmp.$f.$SGE_TASK_ID
		less input.gtf | bedtools nuc -seq -s "-fi" $genomefa -bed stdin | \
			~bilebi00/_CLIP/scripts/calc_motif_freqs_for_nucbed.pl $featureName 19 | \
			awk -F "\t" 'BEGIN{OFS="\t";}{print $1,$4-1,$5,featureName,$20,$7;}' featureName=$featureName > $subfeature
	elif [ $f == G_pct ]; then
		featureNuc=$genomefa
		featureConvertFromGtf=14
		featureStartBase=1
		normFeature=1; 
	elif [ $f == GC_pct ]; then
		featureNuc=$genomefa
		featureConvertFromGtf=11
		featureStartBase=1
	else
		echo unknown feature: $f
		continue;
	fi
	echo ~bilebi00/_CLIP/scripts/get_feature_density.pl input.gtf $subfeature $normFeature $stranded $featureScoreSign $featureScore $featureName $featureStartBase $addFeatureStrand $featureLiftOver $featureConvertFromGtf $featureNuc "'$f.$xx'" $verbose ">>" $statf
	~bilebi00/_CLIP/scripts/get_feature_density.pl input.gtf $subfeature $normFeature $stranded $featureScoreSign $featureScore $featureName $featureStartBase $addFeatureStrand $featureLiftOver $featureConvertFromGtf $featureNuc "'$f.$xx'" $verbose >> $statf
	echo $f $subfeature DONE
done

if [ "`less $statf | wc -l`" -gt 0 ]; then
	R --no-save --args $statf  < ~bilebi00/_CLIP/scripts/get_ttest_stats_for_gtfplus.R | grep -P "^method|^Two" | sort | uniq > significance.$statf
fi

rel=free #xaxis scale
if [ $SGE_TASK_ID == 200 ]; then
	while [ "`ls feature* | wc -l`" -lt 200 ]; do
		echo waiting for other jobs
		sleep 1m
	done

	#merge
	#1
	allstatf=all.feature.stats
	cat feature.stats* > $allstatf
	#2
	allsigstatf=all.significance.stats
	cat significance* | grep "^method" | sort | uniq > $allsigstatf
	cat significance* | grep -v "^method" | grep -v NaN | sort -t $'\t' -k 3,3g -k 7,7g >> $allsigstatf

	#clean
	rm feature.stats* dataprep significance*

	#plot
	omitfeatures="";
	for i in `seq 1 ${#features[@]}`; do
		f=${features[$((i-1))]};
		fjc=${featureJobCounts[$((i-1))]};
		if [ $fjc -gt 1 ]; then
			less $allstatf | grep -w -P "$f" | sed "s/$f.//g" > $f.feature.stats
			less $allsigstatf | grep -w -P "method|$f" | sed "s/$f.//g" > $f.significance.stats
			#			~bilebi00/bin/R --no-save --args $f.feature.stats "$f" "" "" 0 $rel < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
			omitfeatures="$omitfeatures|$f"
		fi
		echo $f done
	done
	omitfeatures="`echo $omitfeatures | cut --bytes 2-`";
	less $allstatf | grep -v -w -P "$omitfeatures" > other.feature.stats
	less $allsigstatf | grep -v -w -P "$omitfeatures" > other.significance.stats
	if [ "`less other.feature.stats | wc -l`" -gt 1 ]; then
		echo mmmmmmmmmmmmmmmmmmmm
		#		~bilebi00/bin/R --no-save --args other.feature.stats "Feature" "" "" 0 $rel < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
	fi

	f=RegionExpression
	if [ "`echo ${features[@]} | grep $f`" != "" ]; then
		less $f.feature.stats | grep -w -P "CLIPPED"  | grep -P "Rep" | \
			awk -F "\t" 'BEGIN{OFS="\t";}{tag="CTRL"; if ($3 !~ /CTRL/) { tag="KD"; } gsub("CLIPPED",tag,$2); gsub("CTRL|EWSR1","",$3); print $0; }' \
			> $f.CTRL_KD.feature.stats
		head -n 1 $allsigstatf > $f.CTRL_KD.significance.stats
		R --no-save --args $f.CTRL_KD.feature.stats 1 < ~bilebi00/_CLIP/scripts/get_ttest_stats_for_gtfplus.R | \
			grep "^Two" | grep -v NaN | sort -t $'\t' -k 3,3g -k 7,7g >> $f.CTRL_KD.significance.stats
	fi

	#TrxExpression
	f=TrxExpression;
	if [ "`echo ${features[@]} | grep $f`" != "" ]; then
		~bilebi00/_xlinkEnrichment/scripts/grep_f_w_like_w_column.pl subregion/common_subregion_genes original.$f 7 1  > $f.CLIPPED

		less $f.CLIPPED | ~bilebi00/_CLIP/scripts/reshape.pl | \
			awk -F "\t" 'BEGIN{OFS="\t";}{tag="CTRL"; if ($3 !~ /CTRL/) { tag="KD"; } if (NR==1) { print $0;} else {print $1,tag,$3,$4;}}' | \
			sed "s/EWSR1_//" | sed "s/CTRL_//" > $f.CTRL_KD.feature.stats
		head -n 1 $allsigstatf > $f.significance.CTRL_KD.feature.stats
		R --no-save --args $f.CTRL_KD.feature.stats 1 < ~bilebi00/_CLIP/scripts/get_ttest_stats_for_lib_type.R | \
			grep "^Two" | grep -v NaN | sort -t $'\t' -k 3,3g -k 7,7g >> $f.significance.CTRL_KD.stats
	fi
fi

date
echo DONE
