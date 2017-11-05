#!/bin/bash

bedfile=$1;
#bedfile=_40_41.bed;
bedfile=_04_14.bed;
tag=PAPD5_IGF

odir=ANNOT
mkdir -p $odir; cd $odir;

#for annot in geneid_genes refseq_genes rmsk rna_genes snomi_genes trna_genes; do
outs="";
for annot in refseq_genes rna_genes snomi_genes trna_genes rmsk; do
	out=${bedfile/.bed/}_${annot/.bed};
	echo $out

	less ~bilebi00/_PAPD5/Normalization/$bedfile | perl -e 'while(<>) { ($c,$b,$e,$id,$score,$s, @o)=split; print join ("\t",$id,$c,$s,$b,$e,$score,@o)."\n"; } ' > $out.inp
	less ~bilebi00/DATA/hg18_ucsc_tracks/$annot.bed | perl -e '$_ = <>; while(<>) { ($c,$b,$e,$id,$score,$s, @o)=split; print join ("\t",$id,$c,$s,$b,$e,$score,@o)."\n"; } ' >> $out.inp
	~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $out.inp > $out.cls 

	~bilebi00/_PAPD5/scripts/annot_from_cls.pl $tag $out.cls > $out.annot;
	outs="$outs $out.annot";
	#echo "Press a key to continue";
	#read;
	rm $out.inp $out.cls
done

allout=${bedfile/.bed/}.csv

echo -e "chromosome\tstrand\tbegin\tend\tz-score\tUCSC RefSeq Genes Track\tUCSC Non-coding RNA Genes Track\tUCSC snoRNA and miRNA Genes Track\tUCSC tRNA Genes Track\tUCSC Repeat Masker Track" > $allout;
paste $outs | awk 'OFS="\t" { print $3,$4,$5,$6,$7,$1,$8,$15,$22,$29;}' >> $allout;
rm $outs;


