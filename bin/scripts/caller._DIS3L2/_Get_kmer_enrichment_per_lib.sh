#!/bin/bash
#$ -S /bin/bash
#$ -q fs_long@@qc_nehalem
#$ -P project_zavolan
#$ -l sjpn=1
#$ -j y
#$ -cwd
#$ -t 1-3
#$ -o LOG._Get_kmer_enrichment_per_lib.sh$TASK_ID

#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.mutrate.pdf */all.*.mutrate.pdf
#gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=all.freqstats.pdf */all.*.freqstats.pdf
export PATH=$HOME/bin:$PATH
#SGE_TASK_ID=4
#-----------------
#TODO set
#c=clusters; N=2000; f="1 / 100";
c=Nucleotide; n=30; N=2000; f="1 /100"; #top 1 percent
b=bottom
b=neighbour
type=exon
#--
db=hg19
genomefa=~bilebi00/DATA/Bowtie_DB_INDEX/$db.fa
genomelen=~bilebi00/aux/human.$db.genome
#-----------------
sid_prot_f=~bilebi00/_DIS3L2/data/protein_sampleid_list
tags=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | cut -f 1 | sort | uniq`)
tag=${tags[$((SGE_TASK_ID-1))]}
ids=(`less $sid_prot_f | grep -w -P "MCLIP|PARCLIP" | grep -w $tag | cut -f 2 | sort`)
idsXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep XL | cut -f 2 | sort`)
idswoXL=(`less $sid_prot_f | grep -w mRNAseq | grep -w background | grep -v XL | cut -f 2 | sort`)

#-----------------
outdir=Analysis/xlinkEnrichment/kmerEnrichment_perlib_$type${c}_${b}_n${n}/$tag
mkdir -p $outdir; pushd $outdir;

#-----------------
nuc=T
id1=${ids[0]}
id2=${ids[1]}

libpat=`echo ${idsXL[@]} | sed 's/ /|/g'`
mu_ref_XL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
libpat=`echo ${idswoXL[@]} | sed 's/ /|/g'`
mu_ref_woXL=`grep -w -P "$libpat" ../../basedOn$c/$tag/all.$tag.freqstats | grep xlink | bedtools groupby -g 2 -c 3 -o mean | awk '{printf("%.6f", $2);}' `
#merged_XL=../../basedOn$c/$tag/full_xlinked_$tag.${id1}_$id2.*.$nuc.xlink_mutrate$mu_ref_XL.bed
merged_XL1=../../basedOn$c/$tag/$tag.$id1.*.$nuc.xlink_mutrate$mu_ref_XL.bed
merged_XL2=../../basedOn$c/$tag/$tag.$id2.*.$nuc.xlink_mutrate$mu_ref_XL.bed

if [ $type == "exon" ]; then
	less $merged_XL1 | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz > merged_XL1
	less $merged_XL2 | bedtools intersect -s -a stdin -b ../../regions/hg19_GMAP_GENEexons_refseq.gtf.gz > merged_XL2
		merged_XL1=merged_XL1
		merged_XL2=merged_XL2
fi



#T=`less $merged_XL | wc -l` 
#N=`echo "$T * $f" | bc`;

for i in 1 2; do
	merged_XL=$merged_XL1;
	if [ $i == 2 ]; then merged_XL=$merged_XL2; fi
	echo Doing $merged_XL

	if [ $c == "clusters" ]; then
		less $merged_XL | head -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > foreground$i.nucbed
		if [ $b == "bottom" ]; then
			less $merged_XL | tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background$i.nucbed
		elif [ $b == "neighbour" ]; then
			cluster_len=`less $merged_XL | awk 'BEGIN{OFS="\t";} {print tag,$3-$2+1}' tag=$tag | bedtools groupby -g 1 -c 2 -o mean | awk '{ print int($2)+1;}'`;
			less $merged_XL | bedtools flank -b $cluster_len -i stdin -g $genomelen | bedtools intersect -v -s -a stdin -b $merged_XL -f 0.1 | head -n $N | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background$i.nucbed
		fi
	elif [ $c == "Nucleotide" ]; then
		cluster_len=$((n*2+1))
		less $merged_XL | bedtools slop -b $n -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$cluster_len > all$i.CCRs
		less all$i.CCRs | head -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > foreground$i.nucbed
		if [ $b == "bottom" ]; then
			less all$i.CCRs |	tail -n $N | bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background$i.nucbed
		elif [ $b == "neighbour" ]; then
			less all$i.CCRs | bedtools flank -b $cluster_len -i stdin -g $genomelen | awk '$3-$2 == len { print; }' len=$cluster_len | \
				bedtools intersect -v -s -a stdin -b all$i.CCRs -f 0.1 | head -n $N | \
				bedtools nuc -s -seq "-fi" $genomefa -bed stdin > background$i.nucbed
		fi
		rm all$i.CCRs
	fi

	echo nucbed finished
	for kmer in 1 4 6; do
		~bilebi00/_EWSR1/scripts/count_kmers.pl foreground$i.nucbed $kmer bed count $tag ${kmer}mer > ${kmer}mers$i 2> ${kmer}mers$i.onames
		~bilebi00/_EWSR1/scripts/count_kmers.pl background$i.nucbed $kmer bed count $tag ${kmer}mer > ${kmer}mers$i.background 2> /dev/null
		R --no-save --args ${kmer}mers$i ${kmer}mers$i.onames ${kmer}mers$i.background < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
		R --no-save --args ${kmer}mers$i ${kmer}mers$i.onames < ~bilebi00/_EWSR1/scripts/kmer_enrichment_w_background_mu_ref.R > /dev/null
		for m in ${kmer}mers$i.*pbeta; do
			less $m | sort -k 5,5g -k 2,2gr > ${m}_
			mv ${m}_ $m;
		done

		echo $kmer finished
	done

done

topn=8
for mfix in backgroundFile_mu_ref uniform_mu_ref; do
	for kmer in 1 4 6; do
		R --no-save --args ${kmer}mers1.$mfix.pbeta ${kmer}mers2.$mfix.pbeta ${kmer}mer.$mfix $topn < ~bilebi00/_DIS3L2/scripts/kmer_rep.R > /dev/null
	done
done

echo DONE
