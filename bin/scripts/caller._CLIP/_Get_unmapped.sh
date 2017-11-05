#===============================================================================
#       USAGE: 
#       PURPOSE: 
#       REQUIREMENTS:  
#       AUTHOR:  Biter Bilen
#       COMPANY:  Biozentrum, University of Basel
#===============================================================================

#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-4
#$ -o LOG._Get_unmapped.$TASK_ID

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
project=ARE
project=EWSR1
project=DIS3L2
project=DIS3L2_mut; analysis_type=CCA; analysis_type=Trun;
project=DIS3L2_mut_reseq; analysis_type=CCA; analysis_type=Trun;

sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
indir=~bilebi00/_CLIP/data/clipz_Unmapped

ids=(`less $sid_prot_f | awk -F "\t" '$9==project && ($8=="tagBackground" || $8 == "tagEnrichment" || $8 == "xlinkEnrichment"){ print $2;}' project=$project`)
annots=(`less $sid_prot_f | awk -F "\t" '$9==project && ($8=="tagBackground" || $8 == "tagEnrichment" || $8 == "xlinkEnrichment"){ print $10;}' project=$project`)
annot=${annots[0]}

#TODO get data if not done before
#for id in ${ids[@]}; do echo $id; mkdir -p $indir; rsync -t mirz@web08:~/BITER/Unmapped/$id.tab.gz $indir/.; done
#for id in ${ids[@]}; do echo $id; mkdir -p $indir; rsync -t ~/BITER/Unmapped/$id.tab.gz $indir/.; done
#exit;

outdir=Analysis/Project_$project/rawdata
mkdir -p $outdir; pushd $outdir

#Trun
if [ $analysis_type == Trun ]; then
	motif=T
	id=${ids[$((SGE_TASK_ID-1))]}
	echo $id
	#>783_1_7 => id_readCopyRank_trimmedTcount
	zless $indir/$id.tab.gz | \
		perl -e '$m=shift; $k=1; while(<>){ next if ($_=~ /^id/); chomp; @t=split; $t[1] =~ s/($m*)$//;  $n=length($1)||0; for(1..$t[2]) { print ">${t[0]}_${_}_${n}\n","$t[1]\n" } $k++; }' $motif | \
		gzip -c >	$id.${motif}run.fasta.gz
	#>783_1_7 => id_readCopyRank_trimmedTcount
	zless $id.${motif}run.fasta.gz | \
		perl -e 'while(<>){ chomp; @t=split/_/; $_=<>; chomp; print join("\t", $_, $t[2] ),"\n";}' | sort -k 1,2 > $id.${motif}run.sorted.tab
	echo $id done
#CCA addition
elif [ $analysis_type == CCA ]; then
	motif=CCA
	echo "$motif ${motif}_over_unmapped ${motif}_over_tRNAmapped ${motif}_count tRNA_count" | sed 's/ /\t/g' > stats
	for f in $indir/*tab.gz; do 
		id=`basename $f .tab.gz`
		echo `zless $f | awk 'NR>1{ print tag"\t"$0}' tag=$id | bedtools groupby -g 1 -c 4 -o sum` `zless $f | grep -P "$motif\t" | awk 'NR>1{ print tag"\t"$0}' tag=$id | bedtools groupby -g 1 -c 4 -o sum` `zless ~/_CLIP/data/clipz_UniqueGenomicAlignments/$id.bedplus.gz | grep tRNA | awk '{ print tag"\t"$0}' tag=$id | bedtools groupby -g 1 -c 6 -o sum` | awk '{ print $1"\t"$4/$2,$4/$6,$4,$6;}'
	done >> stats

	for f in $indir/30?.tab.gz $indir/1???.tab.gz; do 
		id=`basename $f .tab.gz`
		zless $indir/$id.tab.gz | \
			perl -e '$k=1; while(<>){ next if ($_=~ /^id/); @t=split; $a=$t[1]; $t[1] =~ s/CCA|CCAC|CCACC|CCACCA//g; for(1..$t[2]) { print ">${t[0]}_$_\n$t[1]\n" } $k++; }' > $id.CCAp.fasta
		gzip -f $id.CCAp.fasta
	done
fi

exit
#garbage-------------


if [ $SGE_TASK_ID -gt ${#ids[@]} ]; then
	exit;
fi

#TODO write me

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
	#clean #previous run
	rm -rf $paramfile $outdir $imagesfile

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
	zless $indir/$id.bedplus.gz | grep "^chr" | \
		sort -k 1,1 -k 2,2g | \
		gzip -c > $id.bedplus.gz 
else
	echo $id.bedplus.gz exists
	echo exiting
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
#	echo chrom start end name mut strand neigbour_mut copies all_copies mut_MTA mut_MTG SNP | sed 's/ /\t/g' > $id.bed8
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

