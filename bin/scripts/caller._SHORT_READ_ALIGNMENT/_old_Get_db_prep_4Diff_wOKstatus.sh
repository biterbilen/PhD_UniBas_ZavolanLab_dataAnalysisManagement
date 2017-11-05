#!/bin/bash

#TODO run just one at a time
wd=`pwd`;

#set
project=DIS3L2
outdir=Stepanka
sa1=(mRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-0 mRNAseq-HEK-0       totRNAseq-HEK-0_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0 totRNAseq-HEK-0 mRNAseq-HEK-0 mRNAseq-HEK-0 totRNAseq-Hela-0)
sa2=(totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-0 totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2 totRNAseq-Hela-siDis3L2)

project=EWSR1
outdir=Shivendra
sa1=(mRNAseq_siCTRL-A_HeLa)
sa2=(mRNAseq_siEWSR1_HeLa)

project=_SAM68
outdir=__Scheiffele
sa1=(mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_KD_brainStem)
sa2=(mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KD_cerebellum)

dboutdir=~bilebi00/DB
weboutdir=~bilebi00/www/RNAseq/data
dbname=Diff$project
dbfile=${project}_data
soutdir=4DB
mkdir -p $outdir/$soutdir; pushd $outdir/$soutdir;

#prep tables
rm -rf *data
for d in CLIP_comparison_cuffdiff_trusted_superClusters CLIP_comparison_clipz_RNAseq_qn_superClusters; do
	for c in gene isoform splicing promoters tss_group; do 
		if [ $project == "DIS3L2" ]; then
		#	if [ $c != "gene" ] && [ $c != "isoform" ]; then
			if [ $c != "gene" ]; then
				continue;
			fi
#			if [ "`echo $d | grep clipz_RNAseq`" == "" ]; then
#				continue;
#			fi
		fi	 
		tag=CuffdiffTrusted$c;
		if [ "`echo $d | grep clipz`" != "" ]; then 
			if [ $c == "gene" ]; then
				tag=clipzQNgene; 
			else
				continue;
			fi
		fi
		dir=${d}_${c}
		echo $d $c $tag $dir
		for si in `seq 0 $((${#sa1[@]}-1))`; do
			s1=${sa1[$si]}
			s2=${sa2[$si]}
			cat ../$dir/*diff | cut -f 1-14 | awk 'BEGIN{OFS="\t";} { if ($5==s1 && $6==s2 && $7=="OK") { $7=tag; print project,$0; }}' project=$project tag=$tag s1=$s1 s2=$s2 >> $tag.data
		done
	done
done

exit;

#load data
pushd $dboutdir
for i in $wd/$outdir/$soutdir/*data; do
	tag=`basename $i .data`;
	echo $tag
	echo -e '.separator "\t"\ndelete from '$dbfile' where project="'$project'" and ctype="'$tag'";\n.import '$wd/$outdir/$soutdir'/'$tag'.data '$dbfile | sqlite3  $dbname
done
echo -e '.output '$dbfile'\n.separator "\t"\nselect * from '$dbfile';' | sqlite3 $dbname
popd

pushd $weboutdir;
ln -sf $dboutdir/$dbname .

echo DONE
