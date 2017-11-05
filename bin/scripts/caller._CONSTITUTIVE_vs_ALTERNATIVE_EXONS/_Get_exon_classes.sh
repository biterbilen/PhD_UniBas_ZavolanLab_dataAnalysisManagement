#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o log_file

export PERL5LIB=~bilebi00/lib:\$PERL5LIB;

outdir=output_directory

mkdir -p $outdir; cd $outdir;

N=5
clsProg=~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel
p=0.8
outF=exons.GMAP.ge${N}mRNAs.cut$p.gz
accoutF=gene_accessions.GMAP.ge${N}mRNAs.gz
clsinF=mrna_exon_cls.inp
geneAccF=~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/gene_accessions
mrnaP=~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/coordinates/chr
exonP=~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/data/coordinates/EXON/chr
#mouse TODO generalize
geneAccF=~bilebi00/_SAM68/data/mm_auxData/gene_accessions
mrnaP=~bilebi00/_SAM68/data/mm_auxData/coordinates/chr
exonP=~bilebi00/_SAM68/data/coordinates/EXON/chr

#Cluster GMAP mapped exons of geneids having more than 5 mRNAs

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/cluster_exons_of_N_mRNA_having_genes.pl $N $clsinF $geneAccF $mrnaP $exonP $clsProg | gzip -c > $accoutF

~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/classify_exons.pl $clsinF.cls $p $N | gzip -c > $outF

rm $clsinF;
echo Exons are in $outF
echo Transcript gene ids are on $accoutF

echo Total Number of CONSTITUTIVE Exons: `zless $outF | grep CONSTITUTIVE_EXON | wc -l`;
echo Total Number of ALTERNATIVE Exons: `zless $outF | grep ALTERNATIVE_EXON | wc -l`;
echo Total Number of OTHER Exons: `zless $outF | grep OTHER_EXON | wc -l`;
echo Total Number of SINGLE Exons: `zless $outF | grep SINGLE_EXON | wc -l`;
echo Total Number of FIRST Exons: `zless $outF | grep FIRST_EXON | wc -l`;
echo Total Number of LAST Exons: `zless $outF | grep LAST_EXON | wc -l`;

