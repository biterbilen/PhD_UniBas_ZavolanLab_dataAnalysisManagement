#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_short@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-200
#$ -o LOG._Get_features.$TASK_ID

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
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
over=~bilebi00/_EWSR1/data/over/hg18ToHg19.over.chain.gz
subregion=intron
tag=EWSR1
cell=HEK293

date
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_posterior-7_xlinkEnrichment/features_$subregion
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_mutatedInCLIP/$tag/annotation_xlinkEnrichment_topN10000_exptagCTRL/features_$subregion$cell
mkdir -p $outdir; pushd $outdir;

idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk -F "\t" '$1==tag && ($7==cell || cell == "") { print $2; }' tag=$tag cell="$cell" | sort`)

statf=feature.stats.$SGE_TASK_ID
features=(RegionExpression Histone TfbsClustered RHM GC_pct G_pct Conservation Rloop Top1)
featureIs=(11 11 11 20 11 14 7 10 10)
featureJobCounts=(10 10 200 6 1 1 1 1 1)
#needs a hash
#features should be ordered by jobcount fix this requirement 

#feature specific files
#----------------------
TfbsClusteredf=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/wgEncodeRegTfbsClustered/wgEncodeRegTfbsClustered.bed.gz
tfs=(`zless $TfbsClusteredf | cut -f 4 | sort | uniq`)

#FIXME strange parametrization
RegionExpressiond=~bilebi00/_EWSR1/Analysis/Expression/Erfc_intron/
expressionLibs=(`ls $RegionExpressiond/*txt`)

Histonef=~bilebi00/DATA/hg19_ucsc_tracks/ENCODE/Histone_tracks.4download
histonemarks=(`less $Histonef | grep -P "H3K4me1|H3K27ac" | awk '{print $1"/"$2;}'`)

Rloopf=~bilebi00/_EWSR1/data/Ginno2012//molcel_4184_mmc2_all_GC_skewed_regions.txt

Top1f=~bilebi00/_EWSR1/data/Tuduri2009/GSM437526_shTop1.bed.gz

Conservationf=~bilebi00/_CLIP/Analysis/regions/hg19_GMAP_GENE${subregion}s_refseq_wPhastCons_opmedian_ws1000.gtf.gz
Conservationf=~bilebi00/_CLIP/Analysis/regions/hg19_GMAP_GENE${subregion}s_refseq_wPhastCons_opmean_ws1.gtf.gz

recombHotspotMotifs=(CC.CC.T..CC.C CCTCCCT CCCCACCCC G.GG..A.GG.GG AGGGAGG GGGGTGGGG)

#get psitive and negative sets
if [ $SGE_TASK_ID == 1 ]; then
	ln -sf ../$subregion$cell subregion
	ln -sf $Conservationf original.Conservation
	ln -sf $Rloopf original.Rloop
	ln -sf $Top1f original.Top1
	ln -sf $TfbsClusteredf original.TfbsClustered
	ln -sf $Histonef original.Histone
	ln -sf $RegionExpressiond original.RegionExpression
	for id in ${idsf[@]}; do
		idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
		if [ ! -e $id.CLIPPED ] || [ ! -e $id.NOTCLIPPED ]; then
			~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl ~/_CLIP/Analysis/regions/hg19_GMAP_GENE${subregion}s_refseq.gtf.gz subregion/common_subregion_genes | \
				bedtools intersect -s -u -a stdin -b subregion/$id.bed.w_common_subregion_genes > $id.CLIPPED
			~bilebi00/_EWSR1/scripts/grep_gene_list_from_gtf.pl ~/_CLIP/Analysis/regions/hg19_GMAP_GENE${subregion}s_refseq.gtf.gz subregion/common_subregion_genes | \
				bedtools intersect -s -v -a stdin -b subregion/$id.bed.w_common_subregion_genes > $id.NOTCLIPPED
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
			cat $files | sort | uniq -D | uniq > joint.$suff
		else
			echo joint.$suff exists
		fi
	done
	touch dataprep
fi

while [ ! -e dataprep ]; do
	echo waiting for dataprep;
	sleep 1m;
done

idsf=(${idsf[@]} joint)

echo -e "tag\ttype\tmu\tsd\tN" > $statf
#features
for id in ${idsf[@]}; do
	idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
	if [ $id == joint ]; then
		idn=$id
	fi
	for i in `seq 1 ${#features[@]}`; do
		f=${features[$((i-1))]};
		fii=${featureIs[$((i-1))]};
		for file in $id.CLIPPED $id.NOTCLIPPED; do 
			tag=`echo $file | sed "s/$id.//"`;
			if [ $SGE_TASK_ID -gt ${featureJobCounts[$((i-1))]} ]; then
				continue;
			fi
			if [ $f == RegionExpression ]; then
				if [ $SGE_TASK_ID -gt ${#expressionLibs[@]} ]; then
					echo Skip $id $f $fii $file $tag
					continue;
				fi
				n=0;
				x=${expressionLibs[$((SGE_TASK_ID-1))]}
				echo -en $x " "
				rid=`basename $x .txt`
				xx=`less $sid_prot_f | awk -F "\t" '$2==id{print $6}' id=$rid`;
				zless $x | cut -f 1-6 | awk 'BEGIN{OFS="\t"} { if (NR>1) {$2=$2-1; print}} ' | \
					bedtools intersect -wao -a $file -b stdin | \
					awk -F "\t" 'BEGIN{OFS="\t";}{print $0,$14;}' | cut -f 1-9,16- | \
					sort | bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 10,11 -o sum,sum | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ if (n==1) { $fii=($fii/($5-$4+1));} print $0;}' fii=$fii n=$n | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f $xx'" >> $statf
			elif [ $f == Histone ]; then #Histone sites are wo strand
				if [ $SGE_TASK_ID -gt ${#histonemarks[@]} ]; then
					echo Skip $id $f $fii $file $tag
					continue;
				fi
				n=1;
				dirn="`dirname $Histonef`";
				x=${histonemarks[$((SGE_TASK_ID-1))]}
				echo -en $x " "
				xx=`perl -e '$a=shift; $a=~/Histone([\w\d]*)\./; print $1,"\n"' $x`;
				zless $dirn/$x | cut -f 1-6 | \
					bedtools intersect -wao -a $file -b stdin | \
					awk -F "\t" 'BEGIN{OFS="\t";}{print $0,$16*$14;}' | cut -f 1-9,16- | \
					sort | bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 10,11 -o sum,sum | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ if (n==1) { $fii=($fii/($5-$4+1));} print $0;}' fii=$fii n=$n | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f $xx'" >> $statf
			elif [ $f == TfbsClustered ]; then #TfbsClustered sites are wo strand
				if [ $SGE_TASK_ID -gt ${#tfs[@]} ]; then
					echo Skip $id $f $fii $file $tag
					continue;
				fi
				n=1;
				x=${tfs[$((SGE_TASK_ID-1))]}
				echo -en $x " "
				zless original.$f | grep -w $x | cut -f 1-6 | \
					bedtools intersect -wao -a $file -b stdin | \
					awk -F "\t" 'BEGIN{OFS="\t";}{print $0,$16*$14;}' | cut -f 1-9,16- | \
					sort | bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 10,11 -o sum,sum | \
					awk -F "\t" 'BEGIN{OFS="\t";}{ if (n==1) { $fii=($fii/($5-$4+1));} print $0;}' fii=$fii n=$n | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f $x'" >> $statf
			elif [ $f == GC_pct ] || [ $f == G_pct ]; then
				n=0;
				if [ $f == G_pct ]; then n=1; fi
				less $file | bedtools nuc -s -C "-fi" $genomefa -bed stdin | \
					awk -F "\t" 'BEGIN{OFS="\t";}{if (NR>1) { if (n==1) { $fii=($fii/($5-$4+1));} print $0;}}' fii=$fii n=$n | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f'" >> $statf
			elif [ $f == Conservation ]; then
				zless original.$f | awk 'BEGIN{OFS="\t";}{print $1,$4,$5,"tag",0,$7,$13;}' | \
					bedtools intersect -s -a stdin -b $file | \
					bedtools groupby -g 4 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f'" >> $statf
			elif [ $f == Rloop ]; then #Rloop regions are wo strand
				zless original.$f | grep -v "^#" | awk 'BEGIN{c=0;OFS="\t";}{c=c+1; print $1, $2-1, $3, c"_"$4"skew", 0, "+"; }' > $f.hg18.bed
				~bilebi00/www/MARA/bin/liftOver $f.hg18.bed $over $f.bed $f.unmapped
				less $f.bed | \
					bedtools intersect -c -a $file -b stdin | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f'" >> $statf
			elif [ $f == Top1 ]; then #1. Top1 sites are wo strand and 2. introns are reduced to chr1 and chr6 
				zless original.$f > $f.hg18.bed
				~bilebi00/www/MARA/bin/liftOver $f.hg18.bed $over $f.bed $f.unmapped
				less $file | grep -P -w "chr1|chr6" | \
					bedtools intersect -c -a stdin -b $f.bed | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f'" >> $statf
			elif [ $f == RHM ]; then
				if [ $SGE_TASK_ID -gt ${#recombHotspotMotifs[@]} ]; then
					echo Skip $id $f $fii $file $tag
					continue;
				fi
				m=${recombHotspotMotifs[$((SGE_TASK_ID-1))]}
				echo -en $m " "
				less $file | bedtools nuc -seq -s "-fi" $genomefa -bed stdin | \
					~bilebi00/_CLIP/scripts/calc_motif_freqs_for_nucbed.pl $m $((fii-1)) | \
					bedtools groupby -g 2 -c $fii,$fii,$fii -o mean,stdev,count | cut -f 2- | \
					awk 'BEGIN{OFS="\t";}{print tag,type,$0;}' tag="'$idn $tag'" type="'$f $m'" >> $statf
			else
				echo unknown feature: $f
				continue;
			fi
			echo $id $f $fii $file $tag DONE
		done
	done
done

touch runfinished.$SGE_TASK_ID

rel=free #xaxis scale
if [ $SGE_TASK_ID == 200 ]; then
	while [ "`ls runfinished* | wc -l`" -lt 200 ]; do
		echo waiting for other jobs
		sleep 1m
	done

	#merge
	allstatf=all.feature.stats
	cat feature.stats* | sort -r | uniq > $allstatf

	#clean
	rm runfinished* feature.stats* dataprep

	#plot
	omitfeatures="";
	for i in `seq 1 ${#features[@]}`; do
		f=${features[$((i-1))]};
		fjc=${featureJobCounts[$((i-1))]};
		if [ $fjc -gt 1 ]; then
			less $allstatf | grep -w -P "^tag|$f" | sed "s/$f *//g" > $f.feature.stats
			~bilebi00/bin/R --no-save --args $f.feature.stats "$f" "" "" 0 $rel < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
			omitfeatures="$omitfeatures|$f"
		fi
		echo plot done for $f
	done
	omitfeatures="`echo $omitfeatures | cut --bytes 2-`";
	less $allstatf | grep -v -w -P "$omitfeatures" > other.feature.stats
	if [ "`less other.feature.stats | wc -l`" -gt 1 ]; then
		~bilebi00/bin/R --no-save --args other.feature.stats "Feature" "" "" 0 $rel < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
	fi

fi

date
echo DONE
