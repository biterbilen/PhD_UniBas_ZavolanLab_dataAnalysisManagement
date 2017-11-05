#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-8
#$ -o LOG._Get_xlinkEnrichment_annotation.$TASK_ID

export PATH=$HOME/bin:$PATH
#-----------------
#XXX set
sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/annotation/$tag
mkdir -p $outdir; pushd $outdir;

#-----------------
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

#-----------------
posterior=0.000005
posterior=0.05
less $Nucleotide_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > pernuc_merged_XL
#less $clusters_merged_XL | awk '$5<posterior{ print }' posterior=$posterior > clusters_merged_XL

#prep ensembGenes and UCSC tRNAs
#tRNAs are not represented in Ensembl good enough so that added here
zless ~bilebi00/DATA/Ensembl/Annotation/ensGene.gtf.gz > other.gtf;
zless ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs_genename.gtf.gz | sed 's/hg19_tRNAs/tRNA/' >> other.gtf

bedtools intersect -s -a pernuc_merged_XL -b other.gtf -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,$8,$20; }}' | sort -T ~bilebi00/scratch/ | uniq > pernuc_merged_XL.w_type_genename
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,type,$16; }}' type=mRNA | sort -T ~bilebi00/scratch/ | uniq >> pernuc_merged_XL.w_type_genename
bedtools intersect -s -a pernuc_merged_XL -b ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz -wao | \
	awk 'BEGIN{OFS="\t";} { if($7~ /chr/) { print $1,$2,$3,$4,$5,$6,type,$16; }}' type=intron | sort -T ~bilebi00/scratch/ | uniq >> pernuc_merged_XL.w_type_genename

less pernuc_merged_XL.w_type_genename | awk '{ print $7; }' | sort -T ~bilebi00/scratch/ | uniq | while read i; do
	less pernuc_merged_XL.w_type_genename | awk '{ if ($7 == type) { print; } }' type=$i | cut -f 8 | sort -T ~bilebi00/scratch/ | uniq -c | \
		sort -T ~bilebi00/scratch/ -k 1,1gr > pernuc_merged_XL.gene_$i
done


echo -e "type\tgeneCount" > gene.counts
echo -e mRNA"\t""`zless ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz | cut -f 9 | sort | uniq | wc -l `" >> gene.counts
echo -e intron"\t""`zless ../../regions/hg19_GMAP_GENEintrons_refseq.gtf.gz | cut -f 9 | sort | uniq | wc -l`" >> gene.counts
zless other.gtf | cut -f 2,9 | sort | uniq | cut -f 1 | sort | uniq -c | awk 'BEGIN{OFS="\t";}{ print $2,$1}' >> gene.counts

echo -e "type\txlinkedGeneCount" > xlinkedGene.counts
less pernuc_merged_XL.w_type_genename | cut -f 7-8 | sort -T ~bilebi00/scratch/ | uniq | bedtools groupby -g 1 -c 2 -o count | sort -T ~bilebi00/scratch/ -k 2,2gr >> xlinkedGene.counts

total=`less pernuc_merged_XL | wc -l`;
echo -e "type\txlinkCount\ttotalXlinkCount" > xlink.counts
less pernuc_merged_XL.w_type_genename | sort -T ~bilebi00/scratch/ -k 7,7 | bedtools groupby -g 7 -c 7 -o count | sort -T ~bilebi00/scratch/ -k 2,2gr | awk '{ print $1"\t"$2"\t"total;}' total=$total >> xlink.counts

R --no-save --args xlink.counts xlinkedGene.counts gene.counts $tag < ~bilebi00/_DIS3L2/scripts/annot2D.R > /dev/null

echo DONE

