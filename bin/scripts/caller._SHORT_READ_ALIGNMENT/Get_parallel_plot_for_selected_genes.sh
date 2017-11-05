#!/bin/bash

#TODO set
moutdir=Stepanka;
outdir=$moutdir/Parallel_plots
mkdir -p $outdir; pushd $outdir

#genes=(DIS3L2 `less ../CLIP_comparison_cuffdiff_trusted_superClusters_splicing/totRNAseq-Hela-0VStotRNAseq-Hela-siDis3L2_splicing.diff | grep -w yes | sort -k 12,12g | cut -f 3 | head -n 5`)
#ftag=CuffdifTtrustedsplicing_significantTop10_HelaWTvKD
#genes=(DIS3L2 `less ../CLIP_comparison_cuffdiff_trusted_superClusters_promoters/totRNAseq-Hela-0VStotRNAseq-Hela-siDis3L2_promoters.diff | grep -w yes | sort -k 12,12g | cut -f 3 | head -n 5`)
#ftag=CuffdifTtrustedpromoters_significantTop10_HelaWTvKD
genes=(DIS3L2 HSPA1B ACTG1 PMAIP1 MAZ BIRC5 RNF44 SCD EEF2);
ftag=other_clip_targets_in_Hela

genes=(DIS3L2 DIS3L DIS3)
ftag=DIS3s

genes=(DIS3L2)
ftag=DIS3L2

genes=(DIS3L2 EXO1 EXOC4 EXOC8 REXO1L1 EXOSC9 EXOG REXO1 REXO4 REXO2 EXOC5 EXOSC10 EXOC6 EXOSC2 EXOSC5 EXOSC4 EXOC3L EXOC6B EXOC2 EXOC7)
ftag=EXOs

genes=(DIS3L2 ANKRD1 MID1IP1 WHAMM ECSIT BRPF1 PPP1R15A DNAJC3 DUSP4)
ftag=top_anticorrelated_in_clipz

genes=(DIS3L2 ACVR1C PNMAL1 PCDH10 GBX2 HMGCLL1 MAP2 GHR SERPINB9)
ftag=top_correlated_in_clipz

#patsuffix=gene.raw.qn.w_geneSymbol #only gene black and white plots are not suitable for isoform and tss_group
patsuffix=gene.*.w_geneSymbol #only gene black and white plots are not suitable for isoform and tss_group
indir=Expression_matrix;

for f in ../$indir/*$patsuffix; do
	outf=$ftag.`basename $f $patsuffix`
	echo $outf
	head -n 1 $f > $outf
	for g in ${genes[@]}; do
		grep -w $g $f;
	done >> $outf
	~bilebi00/bin/R --no-save --args $outf < ~bilebi00/_SHORT_READ_ALIGNMENT/scripts/parallel_exp_plot_w_scaled_and_real.R > /dev/null
done


echo DONE
