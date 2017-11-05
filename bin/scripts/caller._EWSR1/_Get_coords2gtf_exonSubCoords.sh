#!/bin/bash 

hg19_TR=~bilebi00/_EWSR1/data/hg19_TR.info.all
#hg19
#----
#~bilebi00/_SHORT_READ_ALIGNMENT/scripts/clsformat2gtf.pl "~bilebi00/_EWSR1/data/hg19/coordinates/chr*" $hg19_TR GMAP mRNA > ~bilebi00/_EWSR1/data/hg19_GMAP_mRNA.gtf 2> ~bilebi00/_EWSR1/data/hg19_GMAP_mRNA.gtf.log
#gzip -f ~bilebi00/_EWSR1/data/hg19_GMAP_mRNA.gtf
#~bilebi00/_SHORT_READ_ALIGNMENT/scripts/clsformat2gtf.pl "~bilebi00/_EWSR1/data/hg19/coordinates/INTRON/chr*" $hg19_TR GMAP intron > ~bilebi00/_EWSR1/data/hg19_GMAP_INTROS.gtf 2> ~bilebi00/_EWSR1/data/hg19_GMAP_INTRON.gtf.log
#gzip -f ~bilebi00/_EWSR1/data/hg19_GMAP_INTRON.gtf
#~bilebi00/_SHORT_READ_ALIGNMENT/scripts/clsformat2gtf.pl "~bilebi00/_EWSR1/data/hg19/coordinates/EXON/chr*" $hg19_TR GMAP exon > ~bilebi00/_EWSR1/data/hg19_GMAP_EXON.gtf 2> ~bilebi00/_EWSR1/data/hg19_GMAP_EXON.gtf.log
#gzip -f ~bilebi00/_EWSR1/data/hg19_GMAP_EXON.gtf

~bilebi00/_CLIP/scripts/TRinfo2bed4.pl $hg19_TR > subregions

for r in 5UTR CDS 3UTR; do
	grep -w $r subregions > subregions.$r

	~bilebi00/_CLIP/scripts/grep_append_geneTrx_list_from_gtf.pl ~bilebi00/_EWSR1/data/hg19_GMAP_EXON.gtf.gz subregions.$r 0 0 transcript_id 1 \
		> exons_region_appended.$r

	less exons_region_appended.$r | awk -F "\t" 'BEGIN{OFS="\t"}{ if($7==s){ $3=r; print }}' r=$r s="+" | sort -k 9,9 -k 1,1 -k 4,5g | \
		~bilebi00/_CLIP/scripts/convertTrxRegionCooordinates2Genome.pl > ~bilebi00/_EWSR1/data/hg19_GMAP_$r.gtf
	less exons_region_appended.$r | awk -F "\t" 'BEGIN{OFS="\t";}{ if ($7==s){ $3=r; print }}' r=$r s="-" | sort -k 9,9 -k 1,1 -k 4,5gr | \
		~bilebi00/_CLIP/scripts/convertTrxRegionCooordinates2Genome.pl >> ~bilebi00/_EWSR1/data/hg19_GMAP_$r.gtf

	gzip -f ~bilebi00/_EWSR1/data/hg19_GMAP_$r.gtf
done      

rm -rf subregions* exons_region_appended*
