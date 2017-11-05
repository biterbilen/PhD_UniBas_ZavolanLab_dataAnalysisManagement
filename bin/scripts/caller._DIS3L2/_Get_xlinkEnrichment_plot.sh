#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-1
#$ -o LOG._Get_xlinkEnrichment_plot.Nucleotide.sh.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_plot.clusters.sh.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_plot.exons.sh.$TASK_ID
# $ -o LOG._Get_xlinkEnrichment_plot.introns.sh.$TASK_ID

#for i in LOG._Get_xlinkEnrichment.sh.289.pernuc.*; do echo $i `less $i | grep -w waiting -v | wc -l`; done | grep -v -w 7
#-----------------
export PATH=$HOME/bin:$PATH
#SGE_TASK_ID=2
#-----------------
#TODO set
#c=introns_refseq
#c=exons_refseq
c=clusters
#c=Nucleotide
#-----------------
sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/basedOn$c/$tag
mkdir -p $outdir; pushd $outdir;

#-----------------
if [ $c == "Nucleotide" ]; then
	c=clusters;
fi
#-----------------
if [ "`ls *freqstat 2> /dev/null `" != "" ]; then
	if [ ! -e all.$tag.freqstats ]; then
		echo -e "tag\ttype\tmu\tsd\tN" > all.$tag.freqstats
	fi
	cat $tag*freqstat >> all.$tag.freqstats
	rm -rf $tag*freqstat
	echo `date` freqstats merged
	~bilebi00/bin/R --no-save --args all.$tag.freqstats < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
	echo `date` freqstats plotted
fi

#-----------------
libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
mu_ref_XL=`grep -w -P "$libpat" all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
mu_ref_woXL=`grep -w -P "$libpat" all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `

#-----------------
if [ ! -e all.$tag.mutrate.pdf ]; then
	nuc=T
	for id in ${ids[@]}; do
			~bilebi00/bin/R --no-save --args $tag.$id.$c.${nuc}.xlink_count.gtf.gz $mu_ref_XL < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
			~bilebi00/bin/R --no-save --args $tag.$id.$c.${nuc}.xlink_count.gtf.gz $mu_ref_woXL < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
	done
	
	for id in ${idsXL[@]}; do
			~bilebi00/bin/R --no-save --args $tag.$id.$c.${nuc}.xlink_count.gtf.gz $mu_ref_woXL < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
	done
	
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.$tag.mutrate.pdf $tag.*_mutrate*pdf; 
	echo `date` Merged mutrate plots
fi

nuc=T

id1=${ids[0]}
id2=${ids[1]}
for mu_ref in $mu_ref_XL $mu_ref_woXL; do
	echo merging mu_ref=$mu_ref for $id1 and $id2
	f1=$tag.$id1.$c.$nuc.xlink_mutrate$mu_ref.bed
	f2=$tag.$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	#foreground xlinked
	merged=full_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	top=top_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl min $f1 $f2 $merged $top
	#background not xlinked
	merged=full_notxlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	top=top_notxlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl max $f1 $f2 $merged $top
done
echo `date` all_xlinked and top_xlinked files produced for CLIP

id1=${idsXL[0]}
id2=${idsXL[1]}
for mu_ref in $mu_ref_woXL; do
	echo merging mu_ref=$mu_ref for $id1 and $id2
	f1=$tag.$id1.$c.$nuc.xlink_mutrate$mu_ref.bed
	f2=$tag.$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	#foreground xlinked
	merged=full_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	top=top_xlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl min $f1 $f2 $merged $top
	#background not xlinked
	merged=full_notxlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	top=top_notxlinked_$tag.${id1}_$id2.$c.$nuc.xlink_mutrate$mu_ref.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl max $f1 $f2 $merged $top
done
echo `date` all_xlinked and top_xlinked files produced for XL mRNAseq

echo DONE
exit;
