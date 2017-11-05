#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@high_mem_node
# $ -q fs_long@@qc_nehalem
# $ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-1
#$ -o LOG._Get_xlinkEnrichment_stats.sh

#merge runs
cat *stat | sort -r | uniq > all.stat
R --no-save --args all.stat < ~bilebi00/_EWSR1/scripts/CI.R

#generate trusted taking the worst score from the replicates
c=clusters
mut=0.18
~bilebi00/_PAPD5/scripts/outerJoin.pl 288.$c.xlinkT_mutrate$mut.txt 289.$c.xlinkT_mutrate$mut.txt \
  1,2,3,4,5,6,7,8,9 1,2,3,4,5,6,7,8,9 '1-9,12,10-11,13-14,24-25,27-28' NA > _288_289.$c.xlinkT_mutrate$mut.txt
~bilebi00/bin/R --no-save --args _288_289.$c.xlinkT_mutrate$mut.txt < ~bilebi00/_EWSR1/scripts/get_truested_pbeta.R > /dev/null

#~bilebi00/_PAPD5/scripts/outerJoin.pl 271.introns.xlinkT_mutrate0.05.txt 272.introns.xlinkT_mutrate0.05.txt \
# 1,2,3,4,5,6,7,8,9 1,2,3,4,5,6,7,8,9 '1-9,12,10-11,13-14,24-25,27-28' NA > _271_272.introns.xlinkT_mutrate0.05.txt
#~bilebi00/bin/R --no-save --args _271_272.introns.xlinkT_mutrate0.05.txt < ~bilebi00/_EWSR1/scripts/get_truested_pbeta.R > /dev/null
#
#~bilebi00/_PAPD5/scripts/outerJoin.pl trusted_271* trusted_288* 1,2,3,4,5 1,2,3,4,5 '1-6,12' NA | grep -v -w "NA$" | sort -k 7,7g | grep -w -v NA | wc -l 
#XXX None of the trusted CLIP sites are among the trusted RNAseq sites
#TODO are the trusted RNAseq sample SNPs?

#---------------
#stats for all clusters of EWSR1
bedtools intersect -s -v -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf | wc -l
#1123781
bedtools intersect -s -a EWSR1_clusters.gtf.gz -b ../RNAseq_exonIntronCoverage/hg19_GMAP_GENEintrons.refseq.gtf | wc -l
#1040149
#Half of them are outside the introns

less ~/DATA/hg19_ucsc_tracks/refGenePromoterStart.gtf.gz | bedtools slop -r 0 -l 5000 -s -i stdin -g ~/_DOWNLOAD/_FOR_ALIGNMENT/BEDTools-Version-2.15.0/genomes/human.hg19.genome |  bedtools intersect -s -v -a EWSR1_clusters.gtf.gz -b stdin | wc -l
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

