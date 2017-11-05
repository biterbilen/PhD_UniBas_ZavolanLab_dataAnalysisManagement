#!/bin/bash

#TODO run just one at a time
wd=`pwd`;

#set
project=DIS3L2
outdir=Stepanka
sa1=(mRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-0 mRNAseq-HEK-0       totRNAseq-HEK-0_2 totRNAseq-HEK-0_2 totRNAseq-HEK-0 totRNAseq-HEK-0 mRNAseq-HEK-0 mRNAseq-HEK-0 totRNAseq-Hela-0)
sa2=(totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-0 totRNAseq-HEK-0_2 totRNAseq-HEK-oeDis3L2_2 totRNAseq-HEK-siDis3L2_2 totRNAseq-HEK-oeDis3L2 totRNAseq-HEK-siDis3L2 mRNAseq-HEK-oeDis3L2 mRNAseq-HEK-siDis3L2 totRNAseq-Hela-siDis3L2)

project=SAM68
outdir=Scheiffele
sa1=(mRNAseq_Sam68_KO_brainStem mRNAseq_Sam68_KO_cerebellum mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_KO_brainStem)
sa2=(mRNAseq_Sam68_WT_brainStem mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_WT_cerebellum mRNAseq_Sam68_KO_cerebellum)

project=EMT
outdir=Yoana
sa1=(mRNAseq_NMuMG_EMT_0 mRNAseq_NMuMG_EMT_0 mRNAseq_NMuMG_EMT_0 mRNAseq_NMuMG_EMT_1 mRNAseq_NMuMG_EMT_10 mRNAseq_NMuMG_EMT_10)
sa2=(mRNAseq_NMuMG_EMT_1 mRNAseq_NMuMG_EMT_7 mRNAseq_NMuMG_EMT_10 mRNAseq_NMuMG_EMT_7 mRNAseq_NMuMG_EMT_1 mRNAseq_NMuMG_EMT_7)

project=EWSR1
outdir=Shivendra
sa1=(mRNAseq_EWSR1_siCTRL_1 mRNAseq_EWSR1_siCTRL_2)
sa2=(mRNAseq_EWSR1_siEWSR1_1 mRNAseq_EWSR1_siEWSR1_2)

dboutdir=~bilebi00/DB
weboutdir=~bilebi00/www/RNAseq/data
dbname=Diff
dbfile=${project}_data
soutdir=4DB
mkdir -p $outdir/$soutdir; pushd $outdir/$soutdir;

#prep tables
rm -rf *data
for c in gene isoform splicing promoters tss_group; do 
	tag=CuffdiffTrusted$c;
	for si in `seq 0 $((${#sa1[@]}-1))`; do
		s1=${sa1[$si]}
		s2=${sa2[$si]}
		cat ../Cuffdiff_trusted_RNAseq/$c*diff | awk 'BEGIN{OFS="\t";} { if ($5==s1 && $6==s2 && $7=="OK") { $7=tag; print project,$0; }}' project=$project tag=$tag s1=$s1 s2=$s2 >> $tag.data
	done
done

#load data
pushd $dboutdir
for i in $wd/$outdir/$soutdir/*data; do
	tag=`basename $i .data`;
	echo $tag
	if [ ! -e $dbname ]; then
		echo creating database $dbname
		~bilebi00/DB/create_diff.db.sh $dbname
	fi
	echo -e '.separator "\t"\ndelete from diff_data where project="'$project'" and ctype="'$tag'";\n.import '$wd/$outdir/$soutdir'/'$tag'.data diff_data\n' | sqlite3  $dbname
done
echo -e '.output '$dbfile'\n.separator "\t"\nselect * from diff_data;' | sqlite3 $dbname
popd

pushd $weboutdir;
ln -sf $dboutdir/$dbname .

echo DONE
