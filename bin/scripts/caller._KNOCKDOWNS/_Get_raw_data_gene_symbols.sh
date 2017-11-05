#!/bin/bash

file=Normalization/293wt-RNAseq_HEK293_Nanofectin_RNAseq_HEK293_siAUF1_RNAseq_HEK293_siGFP_RNAseq_HEK293_siTIA1_RNAseq_HEK293_sihnRNPC_RNAseq_mRNASeq_No_4SU_No_XL_rep_A_Clip13_mRNASeq_No_4SU_No_XL_rep_B_Clip13_si_GFP_mRNASEQ_si_HuR_mRNASEQ.raw

#TODO extremely slow; write the file input version
awk 'NR>1 { print $1; }' $file | while read i; do ./scripts/get_gene_symbol_given_rep_trx.pl $i; done > gene_symbols.raw 

#add header to the first line
sed -i '1 i id gene_symbol' gene_symbols.raw 

innerJoin.pl gene_symbols.raw $file 1 1 '2-13' | grep -v "^RP" | cut -f 2-12 > wo_ribosomal_prot_mrnas.raw

