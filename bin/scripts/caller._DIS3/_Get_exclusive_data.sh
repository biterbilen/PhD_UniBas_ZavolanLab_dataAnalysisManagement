#!/bin/bash

outdir=Analysis/EXC
mkdir -p $outdir; pushd $outdir;

for i in DIS3L2 DIS3L DIS3; do 
	cutoff=1
	if [ "`echo $i | grep DIS3L2`" != "" ]; then
		cutoff=4;
	fi
	less ../$i.data | awk '$6>cut && $14 == "EXON" && $16 != "NA"{ print $16}' cut=$cutoff | sort | uniq > $i.genenames
done

less ~bilebi00/_DIS3/data/clipz_RNAseq/RNAseq.raw.qn.w_geneSymbol | cut -f 13 | grep -v gid | grep -v NA | sort | uniq > clipz.genenames

~bilebi00/bin/R --no-save --args DIS3L2.genenames DIS3L.genenames DIS3.genenames .genenames clipz.genenames < ~bilebi00/_DIS3/scripts/venn3.R > /dev/null

sed -i 's/universe/gid/g' venn3.txt

less venn3.txt | grep -P "\t[01]\t0\t0$" > DIS3L2.exclusive
less venn3.txt | grep -P "\t0\t[01]\t0$" > DIS3L.exclusive
less venn3.txt | grep -P "\t0\t0\t[01]$" > DIS3.exclusive

for i in DIS3L2 DIS3L DIS3; do 
	~bilebi00/CDS_CLIP_ElMMo/scripts/grep_f_w.pl $i.exclusive ../$i.data 0 15 > $i.data
done

echo DONE
