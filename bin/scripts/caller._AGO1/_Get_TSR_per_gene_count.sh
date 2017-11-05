#!/bin/bash

outdir=TSR

mkdir -p $outdir; cd $outdir;

region_extension=-500;

less ~bilebi00/_AGO1/data/human_TSR_v1.gff | perl -e '$pr=shift; $_=<>;$_=<>;$_=<>; while(<>){@t=split; print join("\t",$t[9],$t[0],$t[6],$t[3]+1,$t[4]+1,$pr,$t[12])."\n" }' TSR > human_TSR_v1.gff.cls_format

~/CDS_CLIP_ElMMo/scripts/add_geneSymbol.pl "^TSR" human_TSR_v1.gff.cls_format ~bilebi00/CDS_CLIP_ElMMo/data/coordinates/EXON/ ~bilebi00/CDS_CLIP_ElMMo/data/hg_TR.info.all ~bilebi00/SRC/cluster_for_moreMM/generate_oriented_clusters_formutmodel $region_extension > human_TSR_v1.gff.cls_format.w_gene_symbol 

~bilebi00/GEORGES_POLYA/scripts/annotate_sites_with_eventCount.pl human_TSR_v1.gff.cls_format.w_gene_symbol 7 > human_TSR_v1.gff.cls_format.w_gene_symbol.w_event_count

rm human_TSR_v1.gff.cls_format.w_gene_symbol human_TSR_v1.gff.cls_format
