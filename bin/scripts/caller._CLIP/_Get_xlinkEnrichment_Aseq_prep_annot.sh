#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -o LOG._Get_xlinkEnrichment_Aseq_prep_annot.sh

export PATH=$HOME/bin:$PATH

#TODO set
project=DIS3L2
project=EWSR1
project=ARE

#FIXME 
#for those NR_ genes exons are not included in the region based data; but intron,ugs, and dge
#check folder Project_EWSR1/xlinkEnrichment/_ISSUES/NR_exon_problem as an example
exit


RefGene_exp_genes=../RNAseq/GM/RNAseq.w_geneSymbol
RefGene_exp_genes=../RNAseq/expressed.id.gid
RefGene_allexp_genes=../RNAseq/refgene.gid.id.bioType.status.locusLength.regionLength
RefGene_gene_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_gene.gtf.gz
RefGene_exon_regions=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_exon.gtf.gz
flanklen=2500
ENCODE_gene_info_f=~/DATA/hg19_ucsc_tracks/hg19_UCSC_wgEncodeGencodeBasicV12_tab.gz
ENCODE_tracks=(snRNA rRNA) 
ENDODE_genes_pref=~bilebi00/DATA/hg19_ucsc_tracks/hg19_UCSC_wgEncodeGencodeBasicV12_
RefGene_regions_pref=~bilebi00/DATA/hg19_ucsc_tracks/Processed_RefGene/refGene_
ofile=hg19.gtf.gz

outdir=Analysis/Project_$project/annot_regions
mkdir -p $outdir; pushd $outdir;


if [ $project == EWSR1 ] || [ $project == ARE ]; then
	get_exp_genes_subregions=1
	subregions="5UTR CDS 3UTR intron"
	get_ENCODE=0
	get_mRNA=0
	get_tRNA=0
	get_miRNA=0
	get_snoRNA=0
	get_all_exons=0
	if [ $project == ARE ]; then
		get_all_exons=1
	fi
elif [ $project == DIS3L2 ]; then
	get_exp_genes_subregions=0
	get_ENCODE=1
	get_mRNA=1
	get_tRNA=1
	get_miRNA=1
	get_snoRNA=1
	get_all_exons=0
else
	echo project $project parameters are not defined
	exit
fi

#link files
ln -sf $RefGene_exp_genes .
ln -sf $RefGene_allexp_genes .

if [ $get_all_exons == 1 ]; then
	~/_PAPD5/scripts/leftJoin.pl $RefGene_exp_genes $RefGene_allexp_genes 2 1 '' NULL | awk '{ print $4"\t"$3}' > allexpressed.id.gid
	for i in exon; do
		if [ -e $RefGene_regions_pref$i.gtf.gz ]; then
			~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_regions_pref$i.gtf.gz allexpressed.id.gid 0 0 transcript_id | \
				perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
				awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
				gzip -c	> ALL_hg19_UCSC_$i.gtf.gz
			zless ALL_hg19_UCSC_$i.gtf.gz | \
				awk 'BEGIN{OFS="\t"}{ if ( $7 == "+") { print $12,$1,$7,$4,$5,$10} }' | sed 's/";*//g' | \
			 	sort -k 2,2 -k 6,6 -k 1,1 -k 4,4g > ALL_hg19_UCSC_$i.andreasFormat
			zless ALL_hg19_UCSC_$i.gtf.gz | \
				awk 'BEGIN{OFS="\t"}{ if ( $7 == "-") { print $12,$1,$7,$4,$5,$10} }' | sed 's/";*//g' | \
			 	sort -k 2,2 -k 6,6 -k 1,1 -k 4,4gr >> ALL_hg19_UCSC_$i.andreasFormat
			less ALL_hg19_UCSC_$i.andreasFormat | \
				sort -k 2,2 -k 6,6 -k 1,1 | bedtools groupby -g 2,6 -c 4,5,1 -o min,max,distinct -full | \
				awk 'BEGIN{OFS="\t"}{ print $6,$2,$3,$7,$8,$6,$9}' > ALL_hg19_UCSC_geneinfo.andreasFormat
		fi
	done
fi

#expressed RefGene gene subregions
if [ $get_exp_genes_subregions == 1 ]; then
	for i in $subregions; do
		if [ -e $RefGene_regions_pref$i.gtf.gz ]; then
			~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_regions_pref$i.gtf.gz $RefGene_exp_genes 0 0 transcript_id | \
				perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
				awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
				gzip -c	> hg19_UCSC_$i.gtf.gz
		fi
	done
	i=UGS
	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_gene_regions $RefGene_exp_genes 0 0 transcript_id | \
		bedtools flank -s -l $flanklen -r 0 -i stdin -g ~/aux/human.hg19.genome | \
		perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n"; }' | \
		awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
		gzip -c	> hg19_UCSC_$i.gtf.gz
	echo $i
	i=DGE
	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_gene_regions $RefGene_exp_genes 0 0 transcript_id | \
		bedtools flank -s -l 0 -r $flanklen -i stdin -g ~/aux/human.hg19.genome | \
		perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n"; }' | \
		awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
		gzip -c	> hg19_UCSC_$i.gtf.gz
	echo $i
#	i=DGS
#	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl hg19_UCSC_UGS.gtf.gz $RefGene_exp_genes 0 0 transcript_id | \
#		bedtools flank -s -l 0 -r $flanklen -i stdin -g ~/aux/human.hg19.genome | \
#		awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
#		gzip -c	> hg19_UCSC_$i.gtf.gz
#	echo $i
#	i=UGE
#	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl hg19_UCSC_DGE.gtf.gz $RefGene_exp_genes 0 0 transcript_id | \
#		bedtools flank -s -l $flanklen -r 0 -i stdin -g ~/aux/human.hg19.genome | \
#		awk -F "\t" 'BEGIN{OFS="\t"}{$3=toupper(tag); print}' tag=$i | \
#		gzip -c	> hg19_UCSC_$i.gtf.gz
fi

if [ $get_mRNA == 1 ]; then
	#expressed RefGene gene loci
	for i in mRNA; do
		~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl $RefGene_exon_regions $RefGene_exp_genes 0 0 transcript_id | \
			perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
			gzip -c	> hg19_UCSC_$i.gtf.gz
	done
fi

if [ $get_ENCODE == 1 ]; then
	#encode non-coding tracks
	for i in ${ENCODE_tracks[@]}; do
		zless $ENDODE_genes_pref$i.gtf.gz | \
			awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=$i | \
			gzip -c > _hg19_UCSC_$i.gtf.gz
		~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl _hg19_UCSC_$i.gtf.gz $ENCODE_gene_info_f 0 1 transcript_id 1 | \
			perl -e 'while(<>){ chomp; @t=split/\t/; $t[8] .= " gene_name \"$t[10]\";"; print join("\t",@t),"\n"; }' | \
			cut -f 1-9 | gzip -c > hg19_UCSC_$i.gtf.gz
		rm _hg19_UCSC_$i.gtf.gz
	done
fi

if [ $get_tRNA == 1 ]; then
	i=tRNA
	zless ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs.gtf.gz | \
		awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=$i | \
		perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
		gzip -c > hg19_UCSC_$i.gtf.gz
fi

if [ $get_miRNA == 1 ]; then
	i=miRNA
	zless ~bilebi00/DATA/hg19_ucsc_tracks/snomir.gtf.gz | \
		grep "hsa-" | \
		awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=$i | \
		perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
		gzip -c > hg19_UCSC_$i.gtf.gz
fi

if [ $get_snoRNA == 1 ]; then
	i=snoRNA
	zless ~bilebi00/DATA/hg19_ucsc_tracks/snomir.gtf.gz | \
		grep -v "hsa-" | \
		awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=$i | \
		perl -e 'while(<>){ chomp; $_ =~ /gene_id "([^"]+)";/; $_ .= " gene_name \"$1\";"; print $_,"\n";  }' | \
		gzip -c > hg19_UCSC_$i.gtf.gz
fi

#merge all tracks
for i in hg19_UCSC*gz; do
	zless $i;
done | gzip -c > $ofile

echo DONE
exit;


#snomir in UCSC could be more curated than encode non-coding tracks
#REMARK snomir and tRNA
ln -sf ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs.gtf.gz hg19_UCSC_tRNA.gtf.gz
#split or tag
zless ~bilebi00/DATA/hg19_ucsc_tracks/snomir.gtf.gz hg19_UCSC_snomir.gtf.gz | \
	grep "hsa-" | awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=miRNA | gzip -c > hg19_UCSC_miRNA.gtf.gz
zless ~bilebi00/DATA/hg19_ucsc_tracks/snomir.gtf.gz hg19_UCSC_snomir.gtf.gz | \
	grep -v "hsa-" | awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=snoRNA | gzip -c > hg19_UCSC_snoRNA.gtf.gz
zless ~bilebi00/DATA/hg19_ucsc_tracks/tRNAs.gtf.gz hg19_UCSC_tRNA.gtf.gz | \
	awk -F "\t" 'BEGIN{OFS="\t"}{$3=tag; print}' tag=tRNA | gzip -c > hg19_UCSC_tRNA.gtf.gz


echo DONE
