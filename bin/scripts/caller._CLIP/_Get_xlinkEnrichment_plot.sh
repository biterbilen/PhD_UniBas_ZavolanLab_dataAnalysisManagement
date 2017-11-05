#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@qc_nehalem
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 3-5
#$ -o LOG._Get_xlinkEnrichment_plot.sh.$TASK_ID

handler() {
	echo "Job $SGE_TASK_ID receives signal : $1";
}

trap "handler Usr1" SIGUSR1
trap "handler Usr2" SIGUSR2
trap "handler Term" SIGTERM
trap "handler Quit" SIGQUIT
trap "handler  Int" SIGINT
trap "handler  Bus" SIGBUS
trap "handler  Seg" SIGSEGV

export PATH=$HOME/bin:$PATH
#TODO set
#c=Nucleotide
c=mutatedInCLIP
#-----------------
sid_prot_f=~bilebi00/_CLIP/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
idsf=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | awk '$1==tag { print $2; }' tag=$tag | sort`)
idsb=(`less $sid_prot_f | grep -w mRNAseq | grep -w xlinkBackground | cut -f 2 | sort`)
ids=(${idsf[@]} ${idsb[@]})

#-----------------
outdir=Analysis/xlinkEnrichment/basedOnNucleotide_$c/$tag
mkdir -p $outdir; pushd $outdir;

rm -rf runfinished*;

#merge freqstats
#-----------------
if [ "`ls *freqstat 2> /dev/null `" != "" ]; then
	if [ ! -e all.$tag.freqstats ]; then
		echo -e "tag\ttype\tmu\tsd\tN" > all.$tag.freqstats
	fi
	cat $tag*freqstat >> all.$tag.freqstats
	rm -rf $tag*freqstat
	echo Freqstats merged
fi
#plot freqstats
f=all.$tag.freqstats; 
if [ `less $f | wc -l` -gt 30 ]; then #plot separately if there are too many entries 
	less $f | cut -f 2 | sort | uniq | grep -v type | while read i; do 
	grep -P "type|$i" $f > $f.$i; 
	~bilebi00/bin/R --no-save --args $f.$i "Mutation Rate" "" "name.X.mutatedPositions" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null; 
	rm -rf $f.$i
done 
gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=$f.pdf $f\.*\.pdf
rm -rf $f\.*\.pdf
else
	~bilebi00/bin/R --no-save --args $f "Mutation Rate" "" "name.X.mutatedPositions" 1 < ~bilebi00/_EWSR1/scripts/CI.R > /dev/null
fi
echo Freqstats plotted


#calculate pbeta if not done before
#-----------------
nuc=T
rm -rf $tag.pbeta.params
#get the minimum and maximum mut_rate from RNAseq data and estimate
#P[mut_count < copies_count(mu_ref)]
#|||
	#P[mut_count > copies_count(1-mu_ref)]
for id in ${idsf[@]}; do
	p="";
	for idb in ${idsb[@]}; do
		p="$p|$idb.$id\t"
	done
	p=`echo $p | cut --bytes=2-`;
	mu_ref_max=`less all.$tag.freqstats | grep xlink | grep -P $p | sort -k 3,3gr | head -n 1 | awk '{printf("%.6f", $3);}'`
	mu_ref_min=`less all.$tag.freqstats | grep xlink | grep -P $p | sort -k 3,3g  | head -n 1 | awk '{printf("%.6f", $3);}'`
	lib_max=`less all.$tag.freqstats | grep xlink | grep -P $p | sort -k 3,3gr | head -n 1 | awk '{print $1, $2;}'`
	lib_min=`less all.$tag.freqstats | grep xlink | grep -P $p | sort -k 3,3g  | head -n 1 | awk '{print $1, $2;}'`

	#save the parameters in a file
	echo "$lib_min mu_ref_min $mu_ref_min" >> $tag.pbeta.params
	echo "$lib_max mu_ref_max $mu_ref_max" >> $tag.pbeta.params

	mkdir -p xlinkEnrichment tagCount xlinkCount
	R --no-save --args $tag.$id.$id.$nuc.xlink_count.gtf.gz $mu_ref_min min < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
	R --no-save --args $tag.$id.$id.$nuc.xlink_count.gtf.gz $mu_ref_max max < ~bilebi00/_EWSR1/scripts/plot_pbeta_fromGtfPlus.R > /dev/null 
done
echo Calculated xlinkEnrichment Scores

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.$tag.mu_ref.pdf $tag.*_mu_ref*pdf; 
#rm -rf $tag.*_mu_ref*pdf;
#echo Merged mu_ref plots

id1=${idsf[0]}
id2=${idsf[1]}
#patch for samples wo a replicate
if [ "$id2" == "" ]; then
	id2=$id1;
fi
#merge the replicates for xlinkEnrichment min and max mut_rate
#-----------------
odir=xlinkEnrichment
pushd $odir
for mu_ref in min max; do

	revfiles=()
	idns=()
	f1=$tag.$id1.$id1.$nuc.xlink_mu_ref$mu_ref.bed
	f2=$tag.$id2.$id2.$nuc.xlink_mu_ref$mu_ref.bed
	
	#revcum
	~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f1 4 0 0 1 > $f1.reverse
	~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f2 4 0 0 1 > $f2.reverse
	revfiles=(${revfiles[@]} $f1.reverse $f2.reverse)
	idn1="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id1`";
	idn2="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id2`";
	idns=("${idns[@]}" "'$idn1'" "'$idn2'")
	R --no-save --args "$tag.$id1.$id2.$mu_ref.revcum" "-$odir" 1 ${#idns[@]} ${revfiles[@]} "${idns[@]}" < ~bilebi00/_CLIP/scripts/rev_cum.R > /dev/null

	#merge
	if [ $id1 -ne $id2 ]; then
		merged=full_xlinked_$tag.${id1}_${id2}.$nuc.xlink_mu_ref$mu_ref.bed
		echo merging mu_ref_type=$mu_ref for $f1 and $f2
		~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl min $f1 $f2 $merged
	fi
	#merge others too
	for id in ${idsf[@]}; do
		if [ $id == $id1 ] || [ $id == $id2 ]; then
			continue;
		fi
		f1=$merged
		f2=$tag.$id.$id.$nuc.xlink_mu_ref$mu_ref.bed
		merged=full_xlinked_$tag.$nuc.xlink_mu_ref$mu_ref.bed
		echo merging mu_ref_type=$mu_ref for $f1 and $f2
		~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl min $f1 $f2 $merged.tmp
		mv $merged.tmp $merged

		#revcum
		~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f2 4 0 0 1 > $f2.reverse
		revfiles=(${revfiles[@]} $f2.reverse)
		idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
		idns=("${idns[@]}" "'$idn'")
	done

	#singleton processing
	for id in ${idsf[@]}; do
		bedtools closest -s -io -t first -d -a $tag.$id.$id.$nuc.xlink_mu_ref$mu_ref.bed -b $tag.$id.$id.$nuc.xlink_mu_ref$mu_ref.bed | \
			awk '$13<10{print;}' | cut -f 1-6 \
			> $tag.$id.$id.$nuc.xlink_mu_ref$mu_ref.wo_singleton.bed
	done

	echo full_xlinked produced for CLIP
	R --no-save --args $tag.$mu_ref.revcum "-$odir" 1 ${#idns[@]} ${revfiles[@]} "${idns[@]}" < ~bilebi00/_CLIP/scripts/rev_cum.R > /dev/null
done
popd

#merge replicates for xlinkCount and tagCount 
mu_ref=max
for odir in xlinkCount tagCount; do
	pushd $odir

	revfiles=()
	idns=()
	ff1=$tag.$id1.$id1.$nuc.xlink_mu_ref$mu_ref.bed
	ff2=$tag.$id2.$id2.$nuc.xlink_mu_ref$mu_ref.bed
	f1=$tag.$id1.$id1.$nuc.xlink.bed
	f2=$tag.$id2.$id2.$nuc.xlink.bed
	
	if [ $id1 != $id2 ]; then
		merged=full_xlinked_$tag.${id1}_$id2.$nuc.xlink.bed
		echo merging odir=$odir for $ff1 and $ff2
		~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl max $ff1 $ff2 $merged
		mv $ff2 $f2
	fi
	mv $ff1 $f1

	#revcum
	~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f1 4 0 0 0 > $f1.reverse
	~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f2 4 0 0 0 > $f2.reverse
	revfiles=(${revfiles[@]} $f1.reverse $f2.reverse)
	idn1="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id1`";
	idn2="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id2`";
	idns=("${idns[@]}" "'$idn1'" "'$idn2'")
	R --no-save --args "$tag.$id1.$id2.revcum" "$odir" 0 ${#idns[@]} ${revfiles[@]} "${idns[@]}" < ~bilebi00/_CLIP/scripts/rev_cum.R > /dev/null
	#
	#merge others
	for id in ${idsf[@]}; do
		if [ $id == $id1 ] || [ $id == $id2 ]; then
			continue;
		fi
		f1=$merged
		ff2=$tag.$id.$id.$nuc.xlink_mu_ref$mu_ref.bed
		f2=$tag.$id.$id.$nuc.xlink.bed
		merged=full_xlinked_$tag.$nuc.xlink.bed
		echo merging odir=$odir for $f1 and $ff2
		~bilebi00/_EWSR1/scripts/summarizeScores_from2BedFiles.pl max $f1 $ff2 $merged.tmp
		#
		mv $merged.tmp $merged
		mv $ff2 $f2
		#revcum
		~bilebi00/_CONSTITUTIVE_vs_ALTERNATIVE_EXONS/scripts/reverse_cum_dist_of_density.pl $f2 4 0 0 0 > $f2.reverse
		revfiles=(${revfiles[@]} $f2.reverse)
		idn="`less $sid_prot_f | awk -F "\t" '$2==id{print $6; }' id=$id`";
		idns=("${idns[@]}" "'$idn'")
	done
	rm -rf *xlink_mu_ref*;

	#singleton processing
	for id in ${idsf[@]}; do
		bedtools closest -s -io -t first -d -a $tag.$id.$id.$nuc.xlink.bed -b $tag.$id.$id.$nuc.xlink.bed | \
			awk '$13<10{print;}' | cut -f 1-6 \
			> $tag.$id.$id.$nuc.xlink.wo_singleton.bed
	done

	echo full_$odir produced for CLIP
	R --no-save --args $tag.revcum "$odir" 0 ${#idns[@]} ${revfiles[@]} "${idns[@]}" < ~bilebi00/_CLIP/scripts/rev_cum.R > /dev/null
	popd
done

for odir in xlinkEnrichment xlinkCount tagCount; do
	pushd $odir;
	rm -rf *gz
	gzip *bed;
	popd
done

echo DONE
exit;
