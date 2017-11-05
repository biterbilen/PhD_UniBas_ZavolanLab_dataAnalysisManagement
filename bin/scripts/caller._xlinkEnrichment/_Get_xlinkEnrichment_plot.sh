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
annot=mRNA
#TODO set
#c=introns_refseq
#c=exons_refseq
#c=clusters
#c=Nucleotide
c=mutatedInCLIP
#-----------------
sid_prot_f=~bilebi00/_xlinkEnrichment/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_${annot}_$c/$tag
mkdir -p $outdir; pushd $outdir;

#merge freqstats and plot
#-----------------
if [ "`ls *freqstat 2> /dev/null `" != "" ]; then
	if [ ! -e all.$tag.freqstats ]; then
		echo -e "tag\ttype\tmu\tsd\tN" > all.$tag.freqstats
	fi
	cat $tag*freqstat >> all.$tag.freqstats
	rm -rf $tag*freqstat
	echo Freqstats merged
	~bilebi00/bin/R --no-save --args all.$tag.freqstats "Mutation Rate" "" < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
	echo Freqstats plotted
fi

#calculate pbeta if not done before
#-----------------
nuc=T
id1=${ids[0]}
id2=${ids[1]}
v="$id1.$id1|$id1.$id2|$id2.$id1|$id2.$id2"
rm -rf $tag.pbeta.params
#get the minimum and maximum mut_rate from RNAseq data and estimate
#P[mut_count < copies_count(mu_ref)]
#|||
	#P[mut_count > copies_count(1-mu_ref)]
for id in $id1 $id2; do
	mu_ref_max=`less all.$tag.freqstats | grep xlink | grep -P -v $v | grep -P "\.$id\t" | sort -k 3,3gr | head -n 1 | awk '{printf("%.6f", $3);}'`
	mu_ref_min=`less all.$tag.freqstats | grep xlink | grep -P -v $v | grep -P "\.$id\t" | sort -k 3,3g  | head -n 1 | awk '{printf("%.6f", $3);}'`
	lib_max=`less all.$tag.freqstats | grep xlink | grep -P -v $v | grep -P "\.$id\t" | sort -k 3,3gr | head -n 1 | awk '{print $1, $2;}'`
	lib_min=`less all.$tag.freqstats | grep xlink | grep -P -v $v | grep -P "\.$id\t" | sort -k 3,3g  | head -n 1 | awk '{print $1, $2;}'`

	#save the parameters in a file
	echo "$lib_min mu_ref_min $mu_ref_min" >> $tag.pbeta.params
	echo "$lib_max mu_ref_max $mu_ref_max" >> $tag.pbeta.params

	mkdir -p xlinkEnrichment copies_count mut_count
	~bilebi00/bin/R --no-save --args $tag.$id.$id.${nuc}.xlink_count.gtf.gz $mu_ref_min min < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
	~bilebi00/bin/R --no-save --args $tag.$id.$id.${nuc}.xlink_count.gtf.gz $mu_ref_max max < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
done
echo Calculated posterior probabilities that the mutrate in the Binomially distributed sample is less than a baseline mutation rate

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.$tag.mutrate.pdf $tag.*_mutrate*pdf; 
rm -rf $tag.*_mutrate*pdf;
echo Merged mutrate plots

#merge the replicates for min and max mut_rate
#-----------------
for mu_ref in min max; do
	echo merging mu_ref_type=$mu_ref for $id1 and $id2
	odir=xlinkEnrichment
	f1=$odir/$tag.$id1.$id1.$nuc.xlink_mutrate$mu_ref.bed
	f2=$odir/$tag.$id2.$id2.$nuc.xlink_mutrate$mu_ref.bed
	merged=$odir/full_xlinked_$tag.${id1}_$id2.$nuc.xlink_mutrate$mu_ref.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl min $f1 $f2 $merged
	echo full_xlinked produced for CLIP
done

#merge replicates for mut_count and copies_count
for odir in mut_count copies_count; do
	f1=$odir/$tag.$id1.$id1.$nuc.xlink.bed
	f2=$odir/$tag.$id2.$id2.$nuc.xlink.bed
	if [ ! -e $f1 ]; then
		mv $odir/$tag.$id1.$id1.$nuc.xlink_mutrate$mu_ref.bed $f1
		mv $odir/$tag.$id2.$id2.$nuc.xlink_mutrate$mu_ref.bed $f2
		rm $odir/*mutrate*bed
	fi
	merged=$odir/full_xlinked_$tag.${id1}_$id2.$nuc.xlink.bed
	~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl max $f1 $f2 $merged
	echo full_$odir produced for CLIP
done

echo DONE
exit;
