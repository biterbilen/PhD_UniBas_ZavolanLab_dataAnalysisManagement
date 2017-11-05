#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-10
#$ -o LOG._Get_features_significance_fisher.$TASK_ID

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
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
subregion=intron
tag=EWSR1
cell=HEK293

date
#TODO set which type of feature
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_xlinkEnrichment_topN10000_cell$cell/features_fisher_$subregion
originals=../features_ttest_intron
mkdir -p $outdir; pushd $outdir;

idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort`)

statf=feature.stats.$SGE_TASK_ID
features=(FET RegionExpression Rloop Top1 HeLa)
featureJobCounts=(10 10 1 1 1)

if [ "`echo ${features[@]} | grep -w FET`" != "" ]; then
	fets=(`zless $originals/original.FET | grep "^chr" | awk '{ print $4}' | sort | uniq`)
fi

#FIXME strange parametrization
RegionExpressiond=~bilebi00/_CLIP/Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorm_Untreated/RegionExpression_intron
if [ "`echo ${features[@]} | grep -w RegionExpression`" != "" ]; then
	expressionLibs=(`ls $originals/original.RegionExpression/*txt`)
fi

#get positive and negative sets
if [ $SGE_TASK_ID == 1 ]; then
	for id in ${idsf[@]}; do
		idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
		if [ ! -e $id.CLIPPED ] || [ ! -e $id.NOTCLIPPED ]; then
			#patch for nucbed
			if [ "`echo $outdir | grep motif`" != "" ]; then
				zless ../foreground$id.nucbed.gz | \
					awk 'BEGIN{OFS="\t"}{ if (NR>1) { print $1, "GMAP", tag, $2+1, $3, ".", $6, ".", "gene_id \""$8"\"; transcript_id \""$8"\";"} } ' tag=foreground > $id.CLIPPED
				zless ../background$id.nucbed.gz | \
					awk 'BEGIN{OFS="\t"}{ if (NR>1) { print $1, "GMAP", tag, $2+1, $3, ".", $6, ".", "gene_id \""$8"\"; transcript_id \""$8"\";"} } ' tag=background > $id.NOTCLIPPED
			else
				ln -sf ../$subregion/$id.CLIPPED .
				ln -sf ../$subregion/$id.NOTCLIPPED .
			fi
			echo negative and positive sets prep
		else
			echo negative and positive sets exist
		fi
	done
	for suff in CLIPPED NOTCLIPPED; do
		files="";
		for id in ${idsf[@]}; do
			files="$files $id.$suff"							
		done
		if [ ! -e joint.$suff ]; then
			if [ "`echo $outdir | grep motif`" != "" ]; then
				cat $files | sort | uniq > joint.$suff
			else
				ln -sf ../$subregion/joint.$suff .
			fi
		else
			echo joint.$suff exists
		fi
	done
	touch dataprep
fi

while [ ! -e dataprep ]; do
	echo waiting for dataprep;
	sleep 30;
done

idsf=(${idsf[@]} joint)

if [ ! -e $statf ]; then
	echo -e "feature\tsample\tclipped\tnotclipped\tmethod\talternative\tfeatureClipped\tNotFeatureClipped\tFeatureNotClipped\tNotFeatureNotClipped\tpValue\toddsRatio" > $statf
fi
#features
for id in ${idsf[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	if [ $id == joint ]; then
		idn=$id
	fi
	Nclipped=`less $id.CLIPPED | wc -l`
	Nnotclipped=`less $id.NOTCLIPPED | wc -l`
	for i in `seq 1 ${#features[@]}`; do
		f=${features[$((i-1))]};
		fjc=${featureJobCounts[$((i-1))]};
		if [ $fjc -lt $SGE_TASK_ID ]; then
			continue;
		fi
		if [ $f == Rloop ] || [ $f == Top1 ]; then #Rloop regions are wo strand
			echo -en "$f\t$idn\t$Nclipped\t$Nnotclipped\t" >> $statf
			c_feat_clipped=`bedtools intersect -u -a $id.CLIPPED -b $originals/$f.bed | wc -l`
			c_feat_notclipped=`bedtools intersect -u -a $id.NOTCLIPPED -b $originals/$f.bed | wc -l`
			c_notfeat_clipped=`bedtools intersect -v -a $id.CLIPPED -b $originals/$f.bed | wc -l`
			c_notfeat_notclipped=`bedtools intersect -v -a $id.NOTCLIPPED -b $originals/$f.bed | wc -l`
		elif [ $f == HeLa ]; then
			echo -en "$f\t$idn\t$Nclipped\t$Nnotclipped\t" >> $statf
			c_feat_clipped=`zless $originals/original.$f | awk 'NR>1{print}' | cut -f 1-6 | bedtools intersect -s -u -a $id.CLIPPED -b stdin | wc -l`
			c_feat_notclipped=`zless $originals/original.$f | awk 'NR>1{print}' | cut -f 1-6 | bedtools intersect -s -u -a $id.NOTCLIPPED -b stdin | wc -l`
			c_notfeat_clipped=`zless $originals/original.$f | awk 'NR>1{print}' | cut -f 1-6 | bedtools intersect -s -v -a $id.CLIPPED -b stdin | wc -l`
			c_notfeat_notclipped=`zless $originals/original.$f | awk 'NR>1{print}' | cut -f 1-6 | bedtools intersect -s -v -a $id.NOTCLIPPED -b stdin | wc -l`
		elif [ $f == RegionExpression ]; then
			#XXX TODO
			if [ $SGE_TASK_ID -gt ${#expressionLibs[@]} ]; then
				echo Skip $id $f
				continue;
			fi
			x=${expressionLibs[$((SGE_TASK_ID-1))]}
			rid=`basename $x .txt`
			xx=`less $sid_prot_f | awk -F "\t" '$2==id{print $6}' id=$rid`;
			N=`less $x | wc -l | awk '{print int($1/20)}'`; #5%
			echo -en "$f($xx)\t$idn\t$Nclipped\t$Nnotclipped (Nfeature:$N)\t" >> $statf
			c_feat_clipped=`zless $x | cut -f 1-6 | awk 'BEGIN{OFS="\t"} { if (NR>1) {$2=$2-1; print}}' | head -n $N | bedtools intersect -s -u -a $id.CLIPPED -b stdin | wc -l` 
			c_notfeat_clipped=`zless $x | cut -f 1-6 | awk 'BEGIN{OFS="\t"} { if (NR>1) {$2=$2-1; print}}' | tail -n $N | bedtools intersect -s -u -a $id.CLIPPED -b stdin | wc -l` 
			c_feat_notclipped=`zless $x | cut -f 1-6 | awk 'BEGIN{OFS="\t"} { if (NR>1) {$2=$2-1; print}}' | head -n $N | bedtools intersect -s -u -a $id.NOTCLIPPED -b stdin | wc -l` 
			c_notfeat_notclipped=`zless $x | cut -f 1-6 | awk 'BEGIN{OFS="\t"} { if (NR>1) {$2=$2-1; print}}' | tail -n $N | bedtools intersect -s -u -a $id.NOTCLIPPED -b stdin | wc -l` 
		elif [ $f == FET ]; then
			if [ $SGE_TASK_ID -gt ${#fets[@]} ]; then
				echo Skip $id $f
				continue;
			fi
			x=${fets[$((SGE_TASK_ID-1))]}
			echo -en "$f($x)\t$idn\t$Nclipped\t$Nnotclipped\t" >> $statf
			c_feat_clipped=`zless $originals/original.$f | awk 'BEGIN{OFS="\t"} { if($4==x) {print $1,$2,$3,$4,$5,$6;}}' x=$x | bedtools intersect -s -u -a $id.CLIPPED -b stdin | wc -l` 
			c_notfeat_clipped=`zless $originals/original.$f | awk 'BEGIN{OFS="\t"} { if($4==x) {print $1,$2,$3,$4,$5,$6;}}' x=$x | bedtools intersect -s -v -a $id.CLIPPED -b stdin | wc -l` 
			c_feat_notclipped=`zless $originals/original.$f | awk 'BEGIN{OFS="\t"} { if($4==x) {print $1,$2,$3,$4,$5,$6;}}' x=$x | bedtools intersect -s -u -a $id.NOTCLIPPED -b stdin | wc -l` 
			c_notfeat_notclipped=`zless $originals/original.$f | awk 'BEGIN{OFS="\t"} { if($4==x) {print $1,$2,$3,$4,$5,$6;}}' x=$x | bedtools intersect -s -v -a $id.NOTCLIPPED -b stdin | wc -l` 
		else
			echo unknown feature: $f
			continue;
		fi
		R --no-save --args $c_feat_clipped $c_notfeat_clipped $c_feat_notclipped $c_notfeat_notclipped 1 < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher" >> $statf
		echo $id $f DONE
	done
done

rel=free #xaxis scale
if [ $SGE_TASK_ID == 10 ]; then
	while [ "`ls feature.stats* | wc -l`" -lt 10 ]; do
		echo waiting for other jobs
		sleep 30
	done

	#merge
	allstatf=all.feature.stats
	cat feature.stats* | grep -w "method" | sort | uniq > $allstatf
	cat feature.stats* | grep -w -v "method" | sort -k 11,11g -t $'\t' >> $allstatf

	#clean
	rm feature.stats* dataprep

	#plot
	omitfeatures="";
	for i in `seq 1 ${#features[@]}`; do
		f=${features[$((i-1))]};
		fjc=${featureJobCounts[$((i-1))]};
		if [ $fjc -gt 1 ]; then
			less $allstatf | grep -w -P "method|$f" | sed "s/$f *//g" > $f.feature.stats
			omitfeatures="$omitfeatures|$f"
		fi
		echo $f done
	done
	omitfeatures="`echo $omitfeatures | cut --bytes 2-`";
	less $allstatf | grep -v -w -P "$omitfeatures" > other.feature.stats

#	f=RegionExpression
#	if [ "`echo ${features[@]} | grep $f`" != "" ]; then
#		less $f.feature.stats | grep -w -P "tag|CLIPPED"  | grep -P "tag|Rep" | \
#			awk -F "\t" 'BEGIN{OFS="\t";}{tag="CTRL"; if ($3 !~ /CTRL/) { tag="KD"; } if (NR==1) { print $0;} else {print $1,tag,$3,$4;}}' | \
#			sed "s/CTRL'/'/"  | sed "s/EWSR1'/'/" > $f.CTRL_KD.feature.stats
#
#		head -n 1 $allsigstatf > $f.significance.CTRL_KD.feature.stats
#		R --no-save --args $f.CTRL_KD.feature.stats 1 < ~bilebi00/_CLIP/scripts/get_ttest_stats_for_lib_type.R | \
#			grep "^Two" | grep -v NaN | sort -t $'\t' -k 3,3g -k 7,7g >> $f.significance.CTRL_KD.feature.stats
#	fi

fi

date
echo DONE
