#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-300
# $ -o LOG._Get_xlinkEnrichment_pairwise_distance.sh.HuR.$TASK_ID
#$ -o LOG._Get_xlinkEnrichment_pairwise_distance.sh.AGO2.$TASK_ID

#TODO run and compare me with basedOnNucleotide/HuR.xtalk
exit;
#merge
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=xtalk.pdf *xtalk*/*closest.distance.pdf
#-----------------
#XXX TODO set and log filename above
tag=AGO2
#tag=HuR
type=mRNA #TODO implement others types
#-----------------
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#---
export PATH=$HOME/bin:$PATH
#c=clusters
c=Nucleotide
#-----------------
outdir=Analysis/xlinkEnrichment/xtalk/DIS3L2_$tag
mkdir -p $outdir; pushd $outdir;

d1=pernuc_merged_XL.w_type_genename_$type
less ../../annotation/$tag/pernuc_merged_XL.w_type_genename | \
	awk '$7==type{ print; }' type=$type > $d1
d2=pernuc_merged_XL.w_type_genename_DIS3L2$type
less ../../annotation/DIS3L2/pernuc_merged_XL.w_type_genename | \
	awk '$7==type{ print; }' type=$type > $d2
#-----------------
if [ $SGE_TASK_ID == 1 ]; then
	#prep T positions of clipped genes
	less $d1 | cut -f 8 | sort | uniq | sed 's/";*//g' > xlinked.pernuc_merged_XL
	cp ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz .
	gzip -d hg19_GMAP_GENEexons_refseq.gtf.gz
	#XXX takes long
	grep -w -f xlinked.pernuc_merged_XL hg19_GMAP_GENEexons_refseq.gtf > xlinked.hg19_GMAP_GENEexons_refseq.gtf
	gzip xlinked.hg19_GMAP_GENEexons_refseq.gtf
	#---
	pref=xlinked.hg19_GMAP_GENEexons
	bedtools nuc "-fi" $genomefa -bed ${pref}_refseq.gtf.gz -s -seq | cut -f 1-9,19 | awk 'NR>1{ print; }' | \
		gzip -c > ${pref}_refseq.gtf.nucbed.gz
	nuc=T
	otag=_${pref}_refseq.gtf
	~bilebi00/_EWSR1/scripts/split_gtf_nucbed.pl ${pref}_refseq.gtf.nucbed.gz "$nuc" $otag 0
	gzip -f $nuc$otag;

	touch prepfinished

	#original distance distribution
	bedtools closest -s -d -D "a" -t first -a $d2 -b $d1 > closest.$tag

	less closest.$tag | awk '{if ($15<0) { $15=-$15;} print tag"\t"$15;}' tag=closest | \
		bedtools groupby -g 1 -c 2 -o freqdesc | \
		sed 's/,$//g' | sed 's/,/\n/g' | sed 's/closest\t//' | sed 's/:/\t/' | sort -k 1,1g | sed 's/\n0/0/g'  > closest.distance.$tag

	less closest.$tag | awk '{ if ($15 <0 ) { $15=-$15 }; print tag"\t"$15; }' tag=distance | \
		bedtools groupby -g 1 -c 2,2,2,2 -o mean,stdev,mode,count > closest.stats.$tag
	less closest.$tag | awk -F "\t" '{ if ($15 >= -d && $15 <= d) { print;} }' d=100000000 | \
		sort -k 5,5g | cut -f 7 | sort | uniq -c | sort -k 1,1gr > closest.genes.$tag
else
	while [ ! -e prepfinished ]; do
		echo Waiting for prep for 1min!
		sleep 1m
	done
fi

less $d1 |	bedtools shuffle -i stdin -incl T_xlinked.hg19_GMAP_GENEexons_refseq.gtf.gz -g $genomelen > $SGE_TASK_ID.int
bedtools closest -s -d -D "a" -t first -a $SGE_TASK_ID.int -b $d1 > $SGE_TASK_ID.closest
less $SGE_TASK_ID.closest | awk '{if ($15<0) { $15=-$15;} print tag"\t"$15;}' tag=closest | \
	bedtools groupby -g 1 -c 2 -o freqdesc | \
	sed 's/,$//g' | sed 's/,/\n/g' | sed 's/closest\t//' | sed 's/:/\t/' | sort -k 1,1g | sed 's/\n0/0/g'  > $SGE_TASK_ID.closest.distance
rm $SGE_TASK_ID.int $SGE_TASK_ID.closest

touch $SGE_TASK_ID.runfinished

if [ $SGE_TASK_ID == 1 ]; then
	while [ "`ls *.runfinished | wc -l`" -lt 300 ]; do
		echo Waiting 1 minute for the other tasks!
		sleep 1m
	done

	rm -rf *runfinished prepfinished

	~bilebi00/_DIS3L2/scripts/merge_distances_confInt.pl "*closest.distance" > closest.distance.random.$tag

	rm -rf *\.closest.distance 

	~bilebi00/_PAPD5/scripts/outerJoin.pl closest.distance.$tag closest.distance.random.$tag  1 1 '1-2,4-6' 0 | sort -k 1,1g > DIS3L2_exons2${tag}_exons.closest.distance
	R --no-save --args DIS3L2_exons2${tag}_exons.closest.distance < ~bilebi00/_DIS3L2/scripts/closest_confInt.R > /dev/null

	dis3="`less ../../annotation/DIS3L2/pernuc_merged_XL.w_type_genename | awk '{ if ($7 == type) { print;} }' type=$type | wc -l`"
	other="`less ../../annotation/$tag/pernuc_merged_XL.w_type_genename | awk '{ if ($7 == type) { print;} }' type=$type | wc -l`"
	all=13000
	echo "Number of genes in closest.genes :" `less closest.genes.$tag | wc -l` >> closest.stats.$tag;

	for dist in 0 10 100 200 250; do
		echo -en "distance:$dist"
		less closest.$tag | awk -F "\t" '{ if ($15 >= -d && $15 <= d) { print;} }' d=$dist | sort -k 5,5g | cut -f 7 | sort | uniq -c | sort -k 1,1gr > closest.genesd$dist.$tag
		echo -en "\tNumber of intersecting Dima genes:" `grep -w -f ../dima_target_genes closest.genesd$dist.$tag  | wc -l`
		echo -en "\tsignificance of intersection for d=$dist\t"
		both="`less closest.genesd$dist.$tag | wc -l`"
		~bilebi00/bin/R --no-save --args $both $dis3 $other $all < ~bilebi00/_EWSR1/scripts/fisher.ex.t.R | grep "^Fisher"
	done >> closest.stats.$tag

fi

popd
echo DONE
exit;
