#!/bin/bash
#$ -S /bin/bash
#$ -q fs_short@@qc_nehalem
#$ -P project_zavolan
#$ -j y
#$ -cwd
#$ -o LOG._Get_expressed_strict_exons.sh

#mouse
outdir=EXONS_ExpressedStrictBoundary_mm9
indir=EXONS_mm9
mkdir -p $outdir; cd $outdir;
inp=exons

zless ~/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/$indir/exons.GMAP.ge5mRNAs.cut0.8.gz | grep -P "$chr\t$strand" >> $inp
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $inp.cls
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/unify_exon_classes.pl $inp.cls EXON "^SS" > $inp.strict
~bilebi00/CDS_CLIP_ElMMo/scripts/add_geneSymbol.pl "^chr" $inp.strict ~bilebi00/_SAM68/data/mm_auxData/coordinates/EXON/ ~bilebi00/_SAM68/data/mm_TR.info.all > $inp.strict_w_genesymbols
#TODO generate expressed genes for mouse
echo "TODO" > $inp.strict_w_genesymbols.expressed_genes
#innerJoin.pl $inp.strict_w_genesymbols ~/_KNOCKDOWNS/exp.raw.w_gene_symbol 7 1 '1-7' > $inp.strict_w_genesymbols.expressed_genes

exit;

#human #TODO genralize mouse and human
outdir=EXONS_ExpressedStrictBoundary
indir=EXONS
mkdir -p $outdir; cd $outdir;

inp=exons
zless ~/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/$indir/exons.GMAP.ge5mRNAs.cut0.8.gz | grep -P "$chr\t$strand" >> $inp
~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $inp > $inp.cls
~bilebi00/CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/unify_exon_classes.pl $inp.cls EXON "^SS" > $inp.strict
~bilebi00/CDS_CLIP_ElMMo/scripts/add_geneSymbol.pl "^chr" $inp.strict > $inp.strict_w_genesymbols
innerJoin.pl $inp.strict_w_genesymbols ~/_KNOCKDOWNS/exp.raw.w_gene_symbol 7 1 '1-7' > $inp.strict_w_genesymbols.expressed_genes
