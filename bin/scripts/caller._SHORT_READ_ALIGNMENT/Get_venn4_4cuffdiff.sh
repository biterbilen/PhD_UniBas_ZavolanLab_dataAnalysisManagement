#!/bin/bash

outdir=_Yoana/Cuffdiff_combinatorial
mkdir -p $outdir; pushd $outdir;

for i in gene_exp isoform_exp tss_group_exp splicing promoters; do
	echo $i
	less ../Cuffdiff_trusted_RNAseq/$i.diff | grep -P "mRNAseq_NMuMG_EMT_0\tmRNAseq_NMuMG_EMT_7" | grep -w yes | cut -f 3 | perl -e 'while(<>){chomp; (@a)= split/,/,$_; print join("\n",@a), "\n"; }' | grep -v NA | sort | uniq > $i.genenames 
done

less ../Expression_matrix/CuffdiffTrustedgene.raw.w_geneSymbol | awk '{ print $NF; }' | grep -v gid | perl -e 'while(<>){chomp; (@a)= split/,/,$_; print join("\n",@a), "\n"; }' | grep -v NA | sort | uniq > cufflinks.genenames

~bilebi00/bin/R --no-save --args gene_exp.genenames isoform_exp.genenames tss_group_exp.genenames .genenames cufflinks.genenames _expression < ~bilebi00/_DIS3/scripts/venn3.R > /dev/null
~bilebi00/bin/R --no-save --args promoters.genenames isoform_exp.genenames tss_group_exp.genenames .genenames cufflinks.genenames _promoters < ~bilebi00/_DIS3/scripts/venn3.R > /dev/null
~bilebi00/bin/R --no-save --args splicing.genenames isoform_exp.genenames tss_group_exp.genenames .genenames cufflinks.genenames _splicing < ~bilebi00/_DIS3/scripts/venn3.R > /dev/null

echo DONE
