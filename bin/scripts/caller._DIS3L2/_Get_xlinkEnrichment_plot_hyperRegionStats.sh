#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-1
#$ -o LOG._Get_xlinkEnrichment_plot_hyperRegionStats.$TASK_ID

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.mutrate.pdf */all.*.mutrate.pdf
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.freqstats.pdf */all.*.freqstats.pdf
export PATH=$HOME/bin:$PATH
#SGE_TASK_ID=4
#-----------------
#TODO set
skip=0;
#-----------------
sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/hyperRegionStats/$tag
mkdir -p $outdir; pushd $outdir;

#-----------------
if [ $skip == 0 ] ; then 
nuc=T
id1=${ids[0]}
id2=${ids[1]}

c=Nucleotide;
libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
Nucleotide_mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
Nucleotide_mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
Nucleotide_merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.clusters.$nuc.xlink_mutrate$Nucleotide_mu_ref_XL.bed
#Nucleotide_top_XL=../../basedOnNucleotide/$tag/top_xlinked_$tag.${id1}_$id2.clusters.$nuc.xlink_mutrate$Nucleotide_mu_ref_XL.bed

c=clusters;
libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
clusters_mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
clusters_mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
clusters_merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$clusters_mu_ref_XL.bed
#clusters_top_XL=../../basedOnclusters/$tag/top_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$clusters_mu_ref_XL.bed

c=introns_refseq; 
libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
introns_refseq_mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
introns_refseq_mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
introns_refseq_merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$introns_refseq_mu_ref_XL.bed
#introns_refseq_top_XL=../../basedOn$c/$tag/top_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$introns_refseq_mu_ref_XL.bed

c=exons_refseq; 
libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
exons_refseq_mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
exons_refseq_mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
exons_refseq_merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$exons_refseq_mu_ref_XL.bed
#exons_refseq_top_XL=../../basedOn$c/$tag/top_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$exons_refseq_mu_ref_XL.bed


#-----------------
posterior=0.000005
posterior=0.05
less $Nucleotide_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > pernuc_merged_XL
less $clusters_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > clusters_merged_XL
less $introns_refseq_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > introns_refseq_merged_XL
less $exons_refseq_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > exons_refseq_merged_XL

wc -l pernuc_merged_XL clusters_merged_XL introns_refseq_merged_XL exons_refseq_merged_XL

#refseq annotations
#-----------------
#-----------------
#distinct entities (genes)
#-----------------
bedtools slop -r 1000 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > clusters_merged_XL.genes
bedtools intersect -s -a clusters_merged_XL -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > clusters_merged_XL.genes_exon
bedtools intersect -s -a clusters_merged_XL -b ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > clusters_merged_XL.genes_intron
bedtools flank -r 0 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > clusters_merged_XL.genes_promoter
bedtools flank -r 1000 -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > clusters_merged_XL.genes_downstream3UTR
echo -e "\ncountOfRefseqGenes\ttype\tcountOfRefseqGenes_wmoreThanOneCluster"
echo "-----------------------------------------------------------------";
for i in clusters_merged_XL.genes*; do echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; done | sort -k 1,1gr


bedtools slop -r 1000 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > pernuc_merged_XL.genes
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > pernuc_merged_XL.genes_exon
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > pernuc_merged_XL.genes_intron
bedtools flank -r 0 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > pernuc_merged_XL.genes_promoter
bedtools flank -r 1000 -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | awk '$7~ /chr/ { print $6,$16; }' | sort | uniq -c | sort -k1,1gr > pernuc_merged_XL.genes_downstream3UTR
echo -e "\ncountOfRefseqGenes\ttype\tcountOfRefseqGenes_wmoreThanOnePernuc"
echo "---------------------------------------------------------------";
for i in pernuc_merged_XL.genes*; do echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; done | sort -k 1,1gr

#distinct entities (clusters and pernucs)
#-----------------
bedtools slop -r 1000 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > pernuc_merged_XL.entity_mapping_refseq
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > pernuc_merged_XL.entity_mapping_refseq_exon
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > pernuc_merged_XL.entity_mapping_refseq_intron
bedtools flank -r 0 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > pernuc_merged_XL.entity_mapping_refseq_promoter
bedtools flank -r 1000 -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > pernuc_merged_XL.entity_mapping_refseq_downstream3UTR
echo -e "\ncountOfpernuc\ttype";
echo "------------------------";
for i in pernuc_merged_XL pernuc_merged_XL.entity_mapping_refseq*; do echo `wc -l $i`; done | sort -k 1,1gr  

bedtools slop -r 1000 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > clusters_merged_XL.entity_mapping_refseq
bedtools intersect -s -a clusters_merged_XL -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > clusters_merged_XL.entity_mapping_refseq_exon
bedtools intersect -s -a clusters_merged_XL -b ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > clusters_merged_XL.entity_mapping_refseq_intron
bedtools flank -r 0 -l 1000 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > clusters_merged_XL.entity_mapping_refseq_promoter
bedtools flank -r 1000 -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > clusters_merged_XL.entity_mapping_refseq_downstream3UTR
echo -e "\ncountOfclusters\ttype";
echo "------------------------";
for i in clusters_merged_XL clusters_merged_XL.entity_mapping_refseq*; do echo `wc -l $i`; done | sort -k 1,1gr  

#ensembl annotations -wo protein_coding genes for the not annotated clusters 
#-----------------
#-----------------
zless ~/DATA/Ensembl/Annotation/ensGene.gtf.gz | grep -v -P "protein_coding|processed_transcript|pseudogene" > ensGene.noncoding.gtf
#tRNAs are not represented in Ensembl good enough so that added here
zless ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs_genename.gtf.gz >> ensGene.noncoding.gtf

gzip ensGene.noncoding.gtf
#-----------------
less clusters_merged_XL.entity_mapping_refseq | cut -f 1-6 | bedtools subtract -s -a clusters_merged_XL -b stdin | uniq > clusters_merged_XL_entity_not_mapping_refseq
less pernuc_merged_XL.entity_mapping_refseq | bedtools subtract -s -a pernuc_merged_XL -b stdin | uniq > pernuc_merged_XL_entity_not_mapping_refseq

#-----------------
bedtools slop -pct -r 0.1 -l 0.1 -s -i ensGene.noncoding.gtf.gz -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL_entity_not_mapping_refseq -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$8,$20; }}' | uniq > clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding
bedtools slop -pct -r 0.1 -l 0.1 -s -i ensGene.noncoding.gtf.gz -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL_entity_not_mapping_refseq -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$8,$20; }}' | uniq > pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding

less clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | awk '{ print $7; }' | sort | uniq | while read i; do 
	less clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | grep -w $i | cut -f 1-6,8 > \
		clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding_$i
done
less pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | awk '{ print $7; }' | sort | uniq | while read i; do 
	less pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | grep -w $i | cut -f 1-6,8 > \
		pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding_$i
done
echo -e "\ncountOfclusters\ttype_wnoncoding";
echo "----------------------------------";
for i in clusters_merged_XL clusters_merged_XL.entity_mapping_refseq* clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding*; do echo `wc -l $i`; done | sort -k 1,1gr  

echo -e "\ncountOfpernucs\ttype_wnoncoding";
echo "----------------------------------";
for i in pernuc_merged_XL pernuc_merged_XL.entity_mapping_refseq* pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding*; do echo `wc -l $i`; done | sort -k 1,1gr  

less clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | awk '{ print $7; }' | sort | uniq | while read i; do 
	less clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding_$i | cut -f 7 | sort | uniq -c | \
		sort -k 1,1gr > clusters_merged_XL_entity_not_mapping_refseq.noncoding_$i 
done
less pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | awk '{ print $7; }' | sort | uniq | while read i; do 
	less pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding_$i | cut -f 7 | sort | uniq -c | \
		sort -k 1,1gr > pernuc_merged_XL_entity_not_mapping_refseq.noncoding_$i 
done
echo -e "\ncountOfnoncodingGenes\ttype\tcountOfnoncodingGenes_wmoreThanOneClusters"
echo "---------------------------------------------------------------";
for i in clusters_merged_XL_entity_not_mapping_refseq.noncoding*; do echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; done | sort -k 1,1gr

echo -e "\ncountOfnoncodingGenes\ttype\tcountOfnoncodingGenes_wmoreThanOnePernuc"
echo "---------------------------------------------------------------";
for i in pernuc_merged_XL_entity_not_mapping_refseq.noncoding*; do echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; done | sort -k 1,1gr

#mapping antisense refseqs
#-----------------
#-----------------
less clusters_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | cut -f 1-6 | \
	bedtools subtract -s -a clusters_merged_XL_entity_not_mapping_refseq -b stdin | uniq > clusters_merged_XL_entity_not_mapping_refseq_or_noncoding
less pernuc_merged_XL_entity_not_mapping_refseq.entity_mapping_noncoding | cut -f 1-6 | \
	bedtools subtract -s -a pernuc_merged_XL_entity_not_mapping_refseq -b stdin | uniq > pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding

zless ../../regions/hg19_GMAP_GENE_refseq.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{ if ($7 == "-") { $7="+"; } else { $7="-"; } print $0;  }' | \
	bedtools intersect -s -a clusters_merged_XL_entity_not_mapping_refseq_or_noncoding -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq
less clusters_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq | cut -f 7 | sort | uniq -c | sort -k 1,1gr  > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding.antisense_refseq_genes

zless ../../regions/hg19_GMAP_GENE_refseq.gtf.gz | awk -F "\t" 'BEGIN{OFS="\t"}{ if ($7 == "-") { $7="+"; } else { $7="-"; } print $0;  }' | \
	bedtools intersect -s -a pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq
less pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq | cut -f 7 | sort | uniq -c | sort -k 1,1gr > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding.antisense_refseq_genes

#mapping 4000 flanking regions of refseq
#-----------------
#-----------------
less clusters_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq | cut -f 1-6 | \
	bedtools subtract -s -a clusters_merged_XL_entity_not_mapping_refseq_or_noncoding -b stdin | uniq > clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense
less pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding.entity_mapping_antisense_refseq | cut -f 1-6 | \
	bedtools subtract -s -a pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding -b stdin | uniq > pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense


for ext in 2000 3000 4000; do 
bedtools slop -r $ext -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendeddownstream3UTR
less clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendeddownstream3UTR | cut -f 7 | sort | uniq -c | sort -k 1,1gr  > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.${ext}extendeddownstream3UTR_genes

bedtools slop -r $ext -l 0 -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendeddownstream3UTR
less pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendeddownstream3UTR | cut -f 7 | sort | uniq -c | sort -k 1,1gr  > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.${ext}extendeddownstream3UTR_genes
done


for ext in 2000 3000 4000; do 
bedtools slop -r 0 -l $ext -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendedpromoter
less clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendedpromoter | cut -f 7 | sort | uniq -c | sort -k 1,1gr  > \
	clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.${ext}extendedpromoter_genes

bedtools slop -r 0 -l $ext -s -i ../../regions/hg19_GMAP_GENE_refseq.gtf.gz  -g ~bilebi00/aux/human.hg19.genome | \
	bedtools intersect -s -a pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense -b stdin -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$16; }}' | sort | uniq > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendedpromoter
less pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.entity_mapping_refseq_${ext}extendedpromoter | cut -f 7 | sort | uniq -c | sort -k 1,1gr  > \
	pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense.${ext}extendedpromoter_genes
done

echo -e "\nnumberOfclusters_mappingRegionType"
echo "----------------------------------"
wc -l clusters_merged_XL clusters_merged_XL_entity_not_mapping_refseq clusters_merged_XL_entity_not_mapping_refseq_or_noncoding clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense clusters_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense*entity* 

echo -e "\nnumberOfpernucs_mappingRegionType"
echo "----------------------------------"
wc -l pernuc_merged_XL pernuc_merged_XL_entity_not_mapping_refseq pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense  pernuc_merged_XL_entity_not_mapping_refseq_or_noncoding_or_antisense*entity*

echo -e "\ncountOfRefseqGenes\ttype\tcountOfRefseqGenes_wmoreThanOneCluster"
echo "-----------------------------------------------------------------";
for i in clusters*antisense_refseq_genes clusters*000extendeddownstream3UTR_genes clusters*000extendedpromoter_genes; do 
	echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; 
done | sort -k 1,1gr

echo -e "\ncountOfRefseqGenes\ttype\tcountOfRefseqGenes_wmoreThanOnePernuc"
echo "-----------------------------------------------------------------";
for i in pernuc*antisense_refseq_genes pernuc*000extendeddownstream3UTR_genes pernuc*000extendedpromoter_genes; do 
	echo `wc -l $i` `less $i | grep -w -v " 1 " | wc -l `; 
done | sort -k 1,1gr


#-----UCSC-----------
#-----------------
for i in clusters_merged_XL `ls clusters_merged_XL* | grep  "\.entity"`; do 
	less $i | cut -f 1-6 > _$i.ucsc.bed; 
done

for i in pernuc_merged_XL `ls pernuc_merged_XL* | grep  "\.entity"`; do 
	less $i | cut -f 1-6 > _$i.ucsc.bed; 
done


#percent of genic regiosn covered
#-----------------
#-----------------
echo -e "\nwhat  percent of the gene exons and introns are covered by clusters? "
echo "---------------------------------------------------------------------------"
rand=$RANDOM;
zless ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz | awk 'BEGIN{OFS="\t";}{print $10, $5-$4+1; }' | bedtools groupby -g 1 -c 2 -o sum | sed -e 's/";*//g' > $rand.all;
less exons_refseq_merged_XL | awk 'BEGIN{OFS="\t";} {print $4, $3-$2+1;}' > $rand.exons
~/_PAPD5/scripts/leftJoin.pl $rand.exons $rand.all 1 1 '' NULL 1 | awk '{r=$2/$4; print type"\t"$1"\t"r;}' type=exonCoverage | bedtools groupby -g 1 -c 3,3 -o mean,stdev 

zless ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz | awk 'BEGIN{OFS="\t";}{print $10, $5-$4+1; }' | bedtools groupby -g 1 -c 2 -o sum | sed -e 's/";*//g' > $rand.all;
less introns_refseq_merged_XL | awk 'BEGIN{OFS="\t";} {print $4, $3-$2+1;}' > $rand.introns
~/_PAPD5/scripts/leftJoin.pl $rand.introns $rand.all 1 1 '' NULL 1 | awk '{r=$2/$4; print type"\t"$1"\t"r;}' type=intronCoverage | bedtools groupby -g 1 -c 3,3 -o mean,stdev 
rm $rand*

#cluster length statistics
echo -e "\ncluster length statistics";
echo "------------------------------";

for i in clusters_merged_XL clusters_merged_XL.entity*; do 
	less $i | awk 'BEGIN{OFS="\t";}{print type, $3-$2+1; }' type=$i | bedtools groupby -g 1 -c 2,2 -o mean,stdev
done



#where along the refseq gene they map?
#-----------------
#-----------------
echo -e "\nwhere along the refseq gene pernucOnIntrons map"
echo "------------------------------------------------"
less pernuc_merged_XL.entity_mapping_refseq_intron | cut -f 1-6  | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

echo -e "\nwhere along the refseq gene pernucOnExons map"
echo "------------------------------------------------"
less pernuc_merged_XL.entity_mapping_refseq_exon | cut -f 1-6  | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

echo -e "\nwhere along the refseq gene clusterOnIntrons map"
echo "------------------------------------------------"
less clusters_merged_XL.entity_mapping_refseq_intron | cut -f 1-6  | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

echo -e "\nwhere along the refseq gene clusterOnExons map"
echo "------------------------------------------------"
less clusters_merged_XL.entity_mapping_refseq_exon | cut -f 1-6  | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

echo -e "\nwhere along the refseq gene refseq_exon map"
echo "------------------------------------------------"
less exons_refseq_merged_XL | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

echo -e "\nwhere along the refseq gene refseq_intron map"
echo "------------------------------------------------"
less introns_refseq_merged_XL | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | \
	awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

for N in 1 5 200; do
echo -e "\nbins for at most N=$N refseq_exon crosslinked refseq genes"
echo "----------------------------------------------------------------"
less exons_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do 
	grep -w "$i" exons_refseq_merged_XL ; 
done | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done

for N in 1 5 200; do
echo -e "\nbins for at most N=$N refseq_intron crosslinked refseq genes" 
echo "----------------------------------------------------------------"
less introns_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do 
	grep -w "$i" introns_refseq_merged_XL ;
done |  bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done

for N in 1 5 200; do
echo -e "\nwhere are the refseq_introns of genes that has at most N=$N crosslinked refseq_exons"
echo "------------------------------------------------------------------------------------------"
less exons_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do 
	grep -w "$i" introns_refseq_merged_XL ; 
done | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done

for N in 1 5 200; do
echo -e "\nwhere are the refseq_exons of genes that has at most N=$N crosslinked refseq_introns"
echo "------------------------------------------------------------------------------------------"
less introns_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do 
	grep -w "$i" exons_refseq_merged_XL ;
done |  bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done

for BIN in 9 8; do
for N in 1 5 200; do
echo -e "\nwhere do BIN>=$BIN N<=$N crosslinked refseq_exons have binding sites in refseq_introns"
echo "-------------------------------------------------------------------------------------------------"
less exons_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do grep -w "$i" exons_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin>=BIN || bin<=-BIN) { print $4 } }' BIN=$BIN | while read i; do 
	grep -w "$i" introns_refseq_merged_XL ;
done | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done
done

for BIN in 0 1 2; do
for N in 1 5 200; do
echo -e "\nwhere do BIN<=$BIN N<=$N crosslinked refseq_introns have binding sites in refseq_exons"
echo "--------------------------------------------------------------------------------------------------"
less introns_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1<=N{print $2;}' N=$N | while read i; do grep -w "$i" exons_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin<=BIN && bin>=-BIN) { print $4 } }' BIN=$BIN | while read i; do 
	grep -w "$i" exons_refseq_merged_XL ;
done | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'
done
done

for BIN in 0 1 2; do
for N in `seq 1 10`; do
echo -e "\nGenes with BIN<=$BIN N==$N refseq_intron examples"  
echo "--------------------------------------------"
less introns_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1==N{print $2;}' N=$N | while read i; do grep -w "$i" introns_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin<=BIN && bin>=-BIN) { print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6; } }' BIN=$BIN | head
done
done

for BIN in 9 8; do
for N in `seq 1 10`; do
echo -e "\nGenes with BIN>=$BIN N==$N refseq_exon examples"  
echo "--------------------------------------------"
less exons_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1==N{print $2;}' N=$N | while read i; do grep -w "$i" exons_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin>=BIN || bin<=-BIN) { print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6; } }' BIN=$BIN | head
done
done


echo -e "\nwhere does the longest intron of a gene locate -does pbeta get the longest intron?"
echo "---------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz | awk 'BEGIN{OFS="\t";} {print $1,$4-1,$5,$10,$6,$7, $5-$4+1;}' | \
	bedtools groupby -g 4 -c 7 -o max -full | cut -f 1-6 | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); print $1,$2,$3,$4,$5,$6,bin;  }' |   cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3;} print $2, $1, $1/t*100}' | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print $0,c}'

for r in promoter downstream3UTR intron exon; do
#for r in exon; do
	t=clusters_merged_XL.genes_; 
	echo -e "\nexamples of how many $r's in $t* exist for the top $t$r" 
	echo -e "--------------------------------------------------------------"
	fs=`ls $t* | grep -v $r` 
	less $t$r | awk '{ print $3; }'| while read i; do 
		echo `grep $i $t$r $fs`; 
	done | sed "s/$t//g"  | head -n 30
done



echo -e "\n#anchor exons and introns wrt geneEND to check gradient of clusters along the region; where the cluster midpoint lie wrt to the start of the region (intron|exon)"
#anchor regions wrt to gene start and check the gradient of clusters
#for region in introns exons ; do
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "binCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
for binlen in 2 5; do
#if [ $region == introns ]; then binlen=$((binlen*4)); fi
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "+") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+$3+1)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11; } d=x-b; if (bin<0) { d=-d;} bin=int(binlen*(d+1)/l);  print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterMid_diff_regionStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneEnd mnoac=$N r=$region;
done
done
done
done




echo -e "\n#anchor exons and introns wrt geneSTART to check gradient of clusters along the region; where the cluster midpoint lie wrt to the start of the region (intron|exon)"
#anchor regions wrt to gene start and check the gradient of clusters
#for region in introns exons ; do
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "binCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
for binlen in 2 5; do
#if [ $region == introns ]; then binlen=$((binlen*4)); fi
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "-") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=($2+$3+1)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11; } bin=int(binlen*(x-b)/l); if (bin<0) { bin=-bin;} print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterMid_diff_regionStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneStart mnoac=$N r=$region;
done
done
done
done







echo -e "\n#anchor exons and introns wrt geneSTART to check if the cluster ends extend the boundary of the region; IS EWSR1 binding mRNA/pre-mRNA"
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "\nbinCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
binlen=1
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "-") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$3; l=$11-$10+1; b=$10; if ($6=="-") { b=$11; x=$2+1;} d=x-b; if (bin<0) { d=-d;} bin=int(binlen*(d+1)/l); print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterEnd_diff_regionStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneStart mnoac=$N r=$region
done
done
done

echo -e "\n#anchor exons and introns wrt geneSTART to check if the cluster starts extend the boundary of the region; IS EWSR1 binding mRNA/pre-mRNA"
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "\nbinCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
binlen=1
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "-") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$2+1; l=$11-$10+1; e=$11; if ($6=="-") { e=$10; x=$3;} d=x-e; if (bin<0) { d=-d;} bin=int(binlen*(d+1)/l); print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterStart_diff_regionEnd__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneStart mnoac=$N r=$region
done
done
done


echo -e "\n#anchor exons and introns wrt geneEND to check if the cluster ends extend the boundary of the region; IS EWSR1 binding mRNA/pre-mRNA"
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "\nbinCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
binlen=1
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "+") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$3; l=$11-$10+1; b=$10; if ($6=="-") { b=$11; x=$2+1;} d=x-b; if (bin<0) { d=-d;} bin=int(binlen*(d+1)/l); print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterEnd_diff_regionStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneEND mnoac=$N r=$region
done
done
done

echo -e "\n#anchor exons and introns wrt geneEND to check if the cluster starts extend the boundary of the region; IS EWSR1 binding mRNA/pre-mRNA"
for region in introns exons; do
for order in 1 2 3 4 5; do
echo -e "\nbinCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
binlen=1
for N in 1 2 3 5 10 20 100 200; do
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "+") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$2+1; l=$11-$10+1; e=$11; if ($6=="-") { e=$10; x=$3;} d=x-e; if (bin<0) { d=-d;} bin=int(binlen*(d+1)/l); print $1,$2,$3,$4,$5,$6,bin; }' binlen=$binlen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterStart_diff_regionEnd__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,$0,c}' bc=$binlen ro=$order rot=geneEND mnoac=$N r=$region
done
done
done

fi #end skip

echo DONE
exit;


















#XXX
echo -e "\n#anchor exons and introns wrt geneSTART to check gradient of clusters along the region; where the cluster extends wrt to the start of the region (intron|exon)"
for region in exons introns; do
for order in 1 2 3 4 5; do
echo -e "binCount\tregion\tmaxNumberOfClustersAllowed\tbin\tnumberOfRegions\tpercent\tcumulativePercent"
for binlen in 2 5; do
for N in 1 2 3 5 10 20 100 200; do
if [ $region == introns ]; then minreglen=$((minreglen*10)); binlen=$((binlen*2)); fi
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "-") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$3; l=$11-$10+1; b=$10; if ($6=="-") { b=$11; x=$2+1;} bin=int(binlen*(x-b)/l); if (bin<0) { bin=-bin;} if (l>=mrl) { print $1,$2,$3,$4,$5,$6,bin; } }' binlen=$binlen mrl=$minreglen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__clusterEnd_diff_regionStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,mrl,$0,c}' bc=$binlen ro=$order rot=geneStart mnoac=$N mrl=$minreglen r=$region;
done
done
done
done

echo -e "\n#anchor exons wrt to gene END for PA processing related analysis; how much is the cluster (start position) is in wrt the end of the exon"
for region in exons; do
for order in 1 2 3 5 10; do
echo -e "binCount\tregion\tmaxNumberOfClustersAllowed\tbin\tminRegionLen\tnumberOfRegions\tpercent\tcumulativePercent"
for N in 1 2 5 50; do
for binlen in 2 5; do
#for minreglen in 200 500 1000; do
if [ $region == introns ]; then minreglen=$((minreglen*10)); binlen=$((binlen*2)); fi
echo "#------------------------------------------------------------------------------------------------------------------"
zless ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz | bedtools groupby -g 9 -c 1,7,6,4,5 -o distinct,distinct,distinct,collapse,collapse | \
	sed 's/";//g' | sed 's/"//g' | perl -e '$order=shift; $oi=$order-1; while(<>){@t=split/\s/; @begs=split/,/,$t[7]; @ends=split/,/,$t[8]; if ($t[5] eq "+") { @tmp=reverse @begs; @begs=@tmp; @tmp=reverse @ends; @ends=@tmp;} if ($oi < scalar @begs) { print join("\t", $t[4], $begs[$oi], $ends[$oi], $t[1], $t[6], $t[5]),"\n"; }}' $order | \
 	bedtools intersect -s -a clusters_merged_XL -b stdin -wb | \
	cut -f 1-6,10 | sort -k 7,7 | \
	bedtools groupby -g 7 -c 1,6,2,3,4,5 -o distinct,distinct,collapse,collapse,collapse,collapse | \
	perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs <= $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2]), "\n"; }} }' $N | \
	#perl -e '$N=shift; while(<>){ @t=split/\s/; @begs=split/,/,$t[3]; @ends=split/,/,$t[4]; @names=split/,/,$t[5],@scrs=split/,/,$t[6]; if (scalar @begs == $N) { for my $i (0..$#begs) { print join("\t", $t[1],$begs[$i],$ends[$i],$names[$i],$scrs[$i],$t[2], $t[0] ), "\n"; }} }' $N | \
	sort -k 5,5g | \
	bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE${region}_refseq.gtf.gz -wb | sort -k 5,5g | \
	awk 'BEGIN{OFS="\t";} { x=$2+1; l=$11-$10+1; e=$11; if ($6=="-") { e=$10; c=$3;} bin=int(binlen*(e-x)/l); if (bin<0) { bin=-bin;} if (mrl>=l) { print $1,$2,$3,$4,$5,$6,bin; } }' binlen=$binlen mrl=$minreglen | \
	cut -f 7 | sort | uniq -c | awk 'BEGIN{c=0;}{c=c+$1; print $1,$2,c;}' | tac | \
	awk 'BEGIN{OFS="\t"}{if (NR==1) { t=$3; if (t<minNumberOfRegions) { exit; } } print $2, $1, $1/t*100}' minNumberOfRegions=50 | sort -k 3,3gr | awk 'BEGIN{OFS="\t"; c=0;}{c=c+$3; print bc,"binOf__regionEnd_diff_clusterStart__where_regionIs__"r"_"ro"_wrt_"rot,mnoac,mrl,$0,c}' bc=$binlen ro=$order rot=geneEnd mnoac=$N mrl=$minreglen r=$region;
#done
done
done
done
done


echo DONE
exit;


for BIN in 0 1 2; do
for N in `seq 1 10`; do
echo -e "\nBIN<=$BIN N==$N refseq_exon examples"  
echo "--------------------------------------------"
less exons_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1==N{print $2;}' N=$N | while read i; do grep -w "$i" exons_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin<=BIN && bin>=-BIN) { print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6; } }' BIN=$BIN | head
done
done

for BIN in 9 8; do
for N in `seq 1 10`; do
echo -e "\nBIN>=$BIN N==$N refseq_intron examples"  
echo "--------------------------------------------"
less introns_refseq_merged_XL | cut -f 4 | sort | uniq -c | sort -k 1,1gr | awk '$1==N{print $2;}' N=$N | while read i; do grep -w "$i" introns_refseq_merged_XL ; done |   bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENE_refseq.gtf.gz -wb | sort -k 5,5g |   awk 'BEGIN{OFS="\t";} { x=($2+1+$3)/2; l=$11-$10+1; b=$10; if ($6=="-") { b=$11;} bin=int((x-b)/l*10); if (bin>=BIN || bin<=-BIN) { print $1":"$2"-"$3,$1,$2,$3,$4,$5,$6; } }' BIN=$BIN | head
done
done


echo DONE
exit;

#---------------
#stats for all clusters of EWSR1
bedtools intersect -s -v -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf | wc -l
#1123781
bedtools intersect -s -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf | wc -l
#1040149
#Half of them are outside the introns

zless ~/DATA/hg19_ucsc_tracks/refGenePromoterStart.gtf.gz | bedtools slop -r 0 -l 5000 -s -i stdin -g ~/_DOWNLOAD/_FOR_ALIGNMENT/BEDTools-Version-2.15.0/genomes/human.hg19.genome |  bedtools intersect -s -v -a EWSR1_clusters.gtf.gz -b stdin | wc -l
#2049175
less ~/DATA/hg19_ucsc_tracks/refGenePromoterStart.gtf.gz | bedtools slop -r 0 -l 5000 -s -i stdin -g ~/_DOWNLOAD/_FOR_ALIGNMENT/BEDTools-Version-2.15.0/genomes/human.hg19.genome |  bedtools intersect -s -a EWSR1_clusters.gtf.gz -b stdin | wc -l
#177341
#180KB promoter associated EWSR1 clusters exist

bedtools intersect -s -v -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEexons.refseq.gtf | wc -l
#1917142
bedtools intersect -s -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEexons.refseq.gtf | wc -l
#246638
#200KB exon associated EWSR1 clusters exist

#--------------------
#stats for trusted clusters and introns
#XXX
less trusted_288_289.clusters.xlinkT_mutrate0.05.txt | awk 'BEGIN{OFS="\t";}{ if (NR>1) {  print  $1, $2, $3, $5, $6, $3; } }' > trusted_288_289.clusters.bed 
less trusted_288_289.introns.xlinkT_mutrate0.05.txt | awk 'BEGIN{OFS="\t";}{ if (NR>1) {  print  $1, $2, $3, $5, $6, $3; } }' > trusted_288_289.introns.bed 
less trusted_288_289.cls.bed | bedtools intersect -s -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf 

