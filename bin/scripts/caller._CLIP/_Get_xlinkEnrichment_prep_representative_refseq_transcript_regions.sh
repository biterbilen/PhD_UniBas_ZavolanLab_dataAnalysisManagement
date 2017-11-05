#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_xlinkEnrichment_prep_representative_refseq_transcript_regions.sh

export PATH=$HOME/bin:$PATH

#TODO set
outdir=Analysis/clipz_RNAseq/GeneExpression_GMExp_scaleNorminner_Untreated/representative_transcript_regions
mkdir -p $outdir; pushd $outdir;

less ../RNAseq.w_geneSymbol | grep "NM_" | cut -f 1 > exp_trx
index=0
wheader=0
extension=5000

#get only NM_ tagged transcripts which are not gene fusions
for i in mRNA EXON INTRON CDS 3UTR 5UTR; do
	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl ~bilebi00/_EWSR1/data/hg19_GMAP_$i.gtf.gz exp_trx $index $wheader transcript_id | \
		grep "NM_" | awk -F "\t" '$9 !~ /\-/ { print }' | \
		gzip -c > hg19_GMAP_$i.refseq.gtf.gz
done

tag=UGS
zless hg19_GMAP_mRNA.refseq.gtf.gz | \
	awk -F "\t" 'BEGIN{OFS="\t";}{ $3=tag; if ($7=="+") {$5=$4-1;} else {$4=$5+1;} print; }' tag=$tag | \
	bedtools slop -s -l $extension -r 0 -i stdin -g ~/aux/human.hg19.genome | \
	gzip -c > hg19_GMAP_$tag.refseq.gtf.gz

tag=DGE
zless hg19_GMAP_mRNA.refseq.gtf.gz | \
	awk -F "\t" 'BEGIN{OFS="\t";}{ $3=tag; if ($7=="-") {$5=$4-1;} else {$4=$5+1;} print; }' tag=$tag | \
	bedtools slop -s -r $extension -l 0 -i stdin -g ~/aux/human.hg19.genome | \
	gzip -c > hg19_GMAP_$tag.refseq.gtf.gz

#overlapping with gene loci
tag=DGS
zless hg19_GMAP_mRNA.refseq.gtf.gz | \
	awk -F "\t" 'BEGIN{OFS="\t";}{ $3=tag; if ($7=="+") {$5=$4;} else {$4=$5;} print; }' tag=$tag | \
	bedtools slop -s -r $extension -l 0 -i stdin -g ~/aux/human.hg19.genome | \
	gzip -c > hg19_GMAP_$tag.refseq.gtf.gz

tag=UGE
zless hg19_GMAP_mRNA.refseq.gtf.gz | \
	awk -F "\t" 'BEGIN{OFS="\t";}{ $3=tag; if ($7=="-") {$5=$4;} else {$4=$5;} print; }' tag=$tag | \
	bedtools slop -s -l $extension -r 0 -i stdin -g ~/aux/human.hg19.genome | \
	gzip -c > hg19_GMAP_$tag.refseq.gtf.gz



exit

#this script finds non-overlapping gene loci; 
#if any gene overlaps with the locus of another the number of such genes are indicated in the score field of the gtf file
#
#set feature as the gene name in mRNA.gtf file
#use -nms of mergeBed and convert regions to gtf using the number of genes in the cluster as the score of each entry
zless hg19_GMAP_mRNAs.refseq.gtf.gz  | \
	perl -e 'while(<>){$_=~/gene_id \"([^"]*)/; @t=split/\t/; $t[2]=$1; print join ("\t", @t); }' | \
	bedtools merge -s -nms -i stdin | \
	perl -e 'while(<>){ @t=split; @nms=split/;/,$t[3]; %h=(); map { $h{$_}=1; } @nms; $gene_cluster=join(";", sort keys %h); print join ("\t", $t[0], "GMAP", "gene", $t[1]+1, $t[2],scalar keys %h, $t[4], ".", "gene_id \"$gene_cluster\"; transcript_id \"$gene_cluster\";"  ), "\n"; }' | \
	gzip -c > hg19_GMAP_GENE_refseq.gtf.gz


#get non-exon regions as INTRONic regions from GENE.gtf file
bedtools subtract -s -a hg19_GMAP_GENE_refseq.gtf.gz -b hg19_GMAP_EXONS.refseq.gtf.gz | sed -e 's/\tgene\t/\tintron\t/g' | gzip -c > hg19_GMAP_GENEintrons_refseq.gtf.gz

#get non-intron regions as EXONic regions
bedtools subtract -s -a hg19_GMAP_GENE_refseq.gtf.gz -b hg19_GMAP_GENEintrons_refseq.gtf.gz | sed -e 's/\tgene\t/\texon\t/g' | gzip -c > hg19_GMAP_GENEexons_refseq.gtf.gz

#pernuc regions 
if [ $pernuc == 1 ]; then
	db=hg19
	genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
	nucs="A C G T"
	for tag in hg19_GMAP_GENEexons hg19_GMAP_GENEintrons; do
		bedtools nuc "-fi" $genomefa -bed ${tag}_refseq.gtf.gz -s -seq | cut -f 1-9,19 | awk 'NR>1{ print; }' | \
			gzip -c > ${tag}_refseq.gtf.nucbed.gz
		#----------------
		otag=_${tag}_refseq.gtf
		~bilebi00/_EWSR1/scripts/split_gtf_nucbed.pl ${tag}_refseq.gtf.nucbed.gz "$nucs" $otag 0
		for nuc in $nucs; do
			gzip -f $nuc$otag;
		done
	done
fi

#rm -rf hg19_GMAP_mRNAs.refseq.gtf hg19_GMAP_EXONS.refseq.gtf
rm -rf hg19_GMAP_EXONS.refseq.gtf
echo DONE
