#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-20
#$ -o LOG._Get_bothtailed_xlinkEnrichment.$TASK_ID

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

#TODO set
project=DIS3L2
project=EWSR1
project=MethodValidation
project=AREx
project=ARE
project=DIS3L2_mut

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
fraction=1 #sampling fraction for backward compatibility
indir=~bilebi00/_CLIP/data/clipz_UniqueGenomicAlignments

idbs=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 == "xlinkBackground"{ print $2;}' project=$project`)
idb=${idbs[0]}; #takes the first one hard coded
ids=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 == "xlinkEnrichment"{ print $2;}' project=$project`)
annots=(`less $sid_prot_f | awk -F "\t" '$9==project && $8 == "xlinkEnrichment"{ print $10;}' project=$project`)
annot=${annots[0]}

#TODO get data if not done before
#rsync -t mirz@web08:~/BITER/UniqueGenomicAlignments/$id.bedplus.gz $indir/.
#exit;

if [ $SGE_TASK_ID -gt ${#ids[@]} ]; then
	exit;
fi

if [ `echo ${annots[@]} | sed 's/ /\n/g' | sort | uniq | wc -l` -gt 1 ]; then
	echo inconsistent annotation categories for project $project ids:${ids[@]};
	exit;
fi

annotStr=" -wP $annot"
if [ $annot == ALL ]; then
	annotStr=" -wvP bacterial|fungus|vector|Markers_Adaptors|viral"
fi

id=${ids[$((SGE_TASK_ID-1))]}
echo Doing $id $idb in $indir for $annotStr

themedir="/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/"
vs=135
snpT2C=~bilebi00/DATA/Ensembl/Annotation/snp${vs}Common.T2C.gtf.gz
db=hg19
paramfile=all.$id.params.txt
outfile=all.$id.bedplus
imagesfile=all.$id.pdf

tag=`less $sid_prot_f | awk '$2==id{ print $1}' id=$id`

outdir=Analysis/Project_$project/xlinkEnrichment/$tag/
mkdir -p $outdir; pushd $outdir;

if [ $SGE_TASK_ID == 1 ]; then
	#prepare background
	#1. extract tags w annot
	if [ ! -e ../$idb.annot.bedplus.gz ]; then
		zless $indir/$idb.bedplus.gz | \
			grep "^chr" | \
			grep $annotStr | sort -k 1,1 -k 2,2g | \
			gzip -c > ../$idb.annot.bedplus.gz
	fi
	#2. extract T2C positions
	if [ ! -e ../$idb.annot.MTC.bed.gz ]; then
		zless ../$idb.annot.bedplus.gz | \
			~bilebi00/_CLIP/scripts/decode_mutation_from_clipzFormat.pl MTC | \
			sort -k 1,1 -k 2,2g | bedtools merge -s -scores sum -d -1 -i stdin -nms | \
			sort -k 1,1 -k 2,2g | \
			gzip -c > ../$idb.annot.MTC.bed.gz
	fi
	touch ../prep.done
	echo prep done
fi


while [ ! -e ../prep.done ]; do
	echo waiting for ../prep.done
	sleep 1m
done

#prepare bedplus for id
if [ ! -e $id.bedplus.gz ]; then
	echo $id | perl -e '$_=<>; @t=split/_/; print join("\n",@t);' | \
		while read id; do 
			zless $indir/$id.bedplus.gz | grep "^chr"
		done | \
		sort -k 8,8 | bedtools groupby -g 8 -c 5 -o sum -full | \
		awk 'BEGIN{OFS="\t"}{$5=$NF; $NF=null; print }' | \
		sort -k 1,1 -k 2,2g | \
		gzip -c > $id.bedplus.gz 
else
	echo $id.bedplus.gz exists
	echo EXIT
	exit
fi

echo $id.bedplus.gz created

#extract tags w annot
zless $id.bedplus.gz | \
	grep $annotStr | \
	gzip -c > $id.annot.bedplus.gz

echo $id.annot.bedplus.gz created

#get some statistics
ets=`zless $id.bedplus.gz | awk 'BEGIN{t=0;}{ if ($1 ~ /^chr/) { t = t + $c; print t; } }' c=5 | tail -n 1`
echo $id.bedplus.gz Sample_Total_Size $ets | sed 's/ /\t/g' >> $paramfile
ess=`zless $id.annot.bedplus.gz | awk 'BEGIN{t=0;}{ if ($1 ~ /^chr/) { t = t + $c; print t; } }' c=5 | tail -n 1`
echo $id.bedplus.gz Effective_Sample_Size $ess | sed 's/ /\t/g' >> $paramfile

#extract positions
mut=MTC
zless $id.annot.bedplus.gz | \
	~bilebi00/_CLIP/scripts/decode_mutation_from_clipzFormat.pl $mut | \
	bedtools merge -s -scores sum -d -1 -i stdin -nms | \
	sort -k 1,1 -k 2,2g | \
	gzip -c > $id.annot.$mut.bed.gz

echo $id.annot.$mut.bed.gz created

#map intervals onto MTC intervals
#chrom start end name mut strand neigbour_mut copies all_copies mut_back copies_back SNPdb[FALSE/TRUE]
bedtools slop -b 10 -i $id.annot.MTC.bed.gz -g ~bilebi00/aux/human.$db.genome | \
	bedtools map -null 0 -s -o sum -a $id.annot.MTC.bed.gz -b stdin | \
	bedtools map -null 0 -s -o sum -a stdin -b $id.annot.bedplus.gz | \
	bedtools map -null 0 -s -o sum -a stdin -b $id.bedplus.gz | \
	bedtools map -null 0 -s -o sum -a stdin -b ../$idb.annot.MTC.bed.gz | \
	bedtools map -null 0 -s -o sum -a stdin -b ../$idb.annot.bedplus.gz | \
	bedtools map -null 0 -s -c 9 -o count_distinct -a stdin -b $snpT2C | \
	awk -F '\t' 'BEGIN{OFS="\t"}{if ($12 == 0) { $12 = "FALSE" } else { $12 = "TRUE"; } print; }' | \
	gzip -c > $id.bed12.gz

#calculate pbeta
echo "R --no-save --args $id.bed12.gz $id $themedir $paramfile $outfile < ~bilebi00/_CLIP/scripts/calc_pbeta.R > /dev/null"
R --no-save --args $id.bed12.gz $id $themedir $paramfile $outfile < ~bilebi00/_CLIP/scripts/calc_pbeta.R > /dev/null
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$imagesfile $id*.pdf

outfile=$outfile.gz

#get final statistics
epc=`zless $id.annot.MTC.bed.gz | wc -l`
echo $id.bedplus.gz Evaluated_Positions_Count $epc | sed 's/ /\t/g' >> $paramfile
xpc=`zless $outfile | grep "^chr" | awk '$7=="TRUE" && $8==0{ print }' | wc -l`
echo $id.bedplus.gz Crosslinked_Positions_Count $xpc | sed 's/ /\t/g' >> $paramfile

#rm -rf $id\.*

popd

echo DONE

